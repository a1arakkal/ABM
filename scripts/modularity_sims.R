
library(janitor)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(Matrix)

# Load in interactions by time list --------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")

int_and_neighbors_by_t_true <- int_and_neighbors_by_t

# Load in LSHM results --------------------------------------------

load(file = "cns_week_1_2_res.RData")

clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise")))

load(file = "fb_net_LSM_res.RData")
clusters_fb_network <- JANE:::summary.JANE(cns_fits)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise",
                          "clusters_fb_network")))

# Source helper functions ------------------------------------------------------

source("https://raw.githubusercontent.com/a1arakkal/ABM/refs/heads/master/scripts/helper_functions.R")

# Function to adjust modularity ------------------------------------------------

network <- c
assumed_true_communities <- clusters_accounting_for_noise
prob_within <- .5

adjust_modularity <- function(network,
                              assumed_true_communities,
                              prob_within){
  
  for (t in names(network)){
    
    temp <- lapply(names(network[[t]]$neighbors),
                   FUN = function(x){
                     ints <- network[[t]]$neighbors[[x]]
                     ints[as.numeric(names(ints)) <= as.numeric(x)] <- 0
                     clust_id <- assumed_true_communities[x]
                     
                     if(length(ints) > 0 && sum(ints)>0){
                       
                       if(!is.na(clust_id)){
                         
                         within <- intersect(names(ints), names(which(assumed_true_communities == clust_id)))
                         between <- setdiff(names(ints), names(which(assumed_true_communities == clust_id)))
                         keep_within <- within[runif(length(within)) < prob_within]
                         keep_between <- between[runif(length(between)) < 1-prob_within]
                         p1 <- ints[c(keep_within, keep_between)]
                         
                       } else {
                         
                         between <- names(ints)
                         keep_between <- between[runif(length(between)) < 1-prob_within]
                         p1 <- ints[keep_between]
                         
                       }
                       
                       return(p1)
                       
                     } else {
                       return(numeric(0))
                     }
                     
                   })
    
    names(temp) <- names(network[[t]]$neighbors)
    
    temp2 <- lapply( names(temp), function(x) {
      temp1 <- unlist(lapply(temp, function(y){if(!is.na(y[x])){y[[x]]}}))
      combined <- c(temp[[x]],temp1)
      out <- tapply(combined, names(combined), max)
      out1 <- as.vector(out)       
      names(out1) <- names(out)     
      return(out1[order(as.numeric(names(out1)))])
    })
    
    names(temp2) <- names(temp)
    
    network[[t]]$neighbors<- lapply(temp2,
                                    FUN = function(x){x})
    
    network[[t]]$n_total_t <- sapply(temp2,
                                     FUN = function(x){sum(x)})
    
  }
  
  # all.equal(int_and_neighbors_by_t_true, network)
  
  modularity <- sapply(names(network), 
                       function(x){
                         
                         A <- matrix(0, nrow = total_actors, ncol= total_actors)
                         colnames(A) <- 1:total_actors
                         rownames(A) <- 1:total_actors
                         dat <- network[[x]]$neighbors
                         
                         for(i in names(dat)){
                           A[i, names(dat[[i]])] <- unname(dat[[i]])
                         }
                         # image((as.matrix(A>0)*1)[names(sort(assumed_true_communities)), names(sort(assumed_true_communities))])
                         
                         sum_within <- 0
                         for (i in unique(assumed_true_communities)){
                           sum_within <- sum_within + sum(A[names(which(assumed_true_communities == i)),
                                                            names(which(assumed_true_communities == i))])
                         }
                         
                         return(sum_within/sum(A))
                       })
  
  names(modularity) <- names(network)
  
  return(list(modularity = modularity,
              network = network))
  
}


test <- adjust_modularity(network = int_and_neighbors_by_t_true,
                          assumed_true_communities = clusters_accounting_for_noise,
                          prob_within = .5)

# Parameters ------------------------------------------------------------------

# prob of transmisson among infected
p_infected <- 0.0025

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- c("asymptomatic" = 10L,
                        "pre_symptomatic" = 5L, # if 0 not pre_symptomatic period
                        "symptomatic" = 5L)

# prob of being asymptomatic given infected
p_asym <- 0.8

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector of length 1 will isolate infected and their 1-hop neighbors
quarantine_days_seq  <- 5 # ignored if cluster is NULL (i.e., no intervention)

# N time-steps
timesteps <- names(int_and_neighbors_by_t)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks

# Minimum degree for seed options if 0 will randomly select from 1-692
min_degree_t1 <- 0

# How far to look back for digital contact tracing if used as intervention, inclusive of current day
digital_contact_tracing_look_back <- 4

# Sensitivity and specificity of digital contact tracing if used as intervention (allows for false positives and false negatives)
DCT_sensitivity <- 1
DCT_specificity <- 1

n_trial <- 1000
cores <- 100

set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_for_R0 <- parallel::mclapply(1:n_trial,
                                     FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                                      mean_exposure_days = 1000, # makes it so seed is the only infective for duration of ABM
                                                                      mean_infected_days = mean_infected_days,
                                                                      actor_labels = actor_labels,
                                                                      min_degree_t1 = min_degree_t1,
                                                                      timesteps = timesteps,
                                                                      n_repeat = n_repeat,
                                                                      p_asym = p_asym,
                                                                      clusters = NULL,
                                                                      quarantine_days = quarantine_days,
                                                                      digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                      DCT_sensitivity = DCT_sensitivity,
                                                                      DCT_specificity = DCT_specificity)},
                                     mc.preschedule = TRUE,
                                     mc.cores = cores)


