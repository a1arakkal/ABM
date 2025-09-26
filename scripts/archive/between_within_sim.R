# Load in interactions by time list --------------------------------------------

setwd("/Users/atlan/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")

# Load in LSHM results --------------------------------------------

load(file = "cns_week_1_2_res.RData")

clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t",
                          "total_actors",
                          "clusters_accounting_for_noise",
                          "clusters_dichotomize_g_1")))

int_and_neighbors_by_t_true <- int_and_neighbors_by_t

# Source helper functions ------------------------------------------------------

source("https://raw.githubusercontent.com/a1arakkal/ABM/refs/heads/master/scripts/helper_functions.R")

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
p_asym <- 0 # if 0 all infections are symptomatic

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector of length 1 will isolate infected and their 1-hop neighbors
quarantine_days <- 5 # ignored if cluster is NULL (i.e., no intervention)

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
DCT_specificity_seq <- 1


n_trial <- 1e3
cores <- 100L


# Only allow for interaction within clusters
int_and_neighbors_by_t_only_within <- int_and_neighbors_by_t_true
for (t in names(int_and_neighbors_by_t_only_within)){
  temp <- lapply(names(int_and_neighbors_by_t_only_within[[t]]$neighbors),
                                                  FUN = function(x){
                                                    ints <- int_and_neighbors_by_t_only_within[[t]]$neighbors[[x]]
                                                    clust_id <- clusters_accounting_for_noise[x]
                                                    if(!is.na(clust_id)){
                                                      temp1 <- intersect(names(ints), names(which(clusters_accounting_for_noise == clust_id)))
                                                      p1 <- ints[temp1]
                                                      p2 <- sum(p1) 
                                                      
                                                      return(list(p1 = p1,
                                                                  p2 = p2))
                                                    } else {
                                                      return(list(p1 = numeric(0),
                                                                  p2 = 0))
                                                    }
                                                  })
  names(temp) <- names(int_and_neighbors_by_t_only_within[[t]]$neighbors)
  
  int_and_neighbors_by_t_only_within[[t]]$neighbors<- lapply(temp,
                                                            FUN = function(x){x$p1})
  
  int_and_neighbors_by_t_only_within[[t]]$n_total_t <- sapply(temp,
                                          FUN = function(x){x$p2})
                                          
}

# Only allow for interaction between clusters
int_and_neighbors_by_t_only_between <- int_and_neighbors_by_t_true
for (t in names(int_and_neighbors_by_t_only_between)){
  temp <- lapply(names(int_and_neighbors_by_t_only_between[[t]]$neighbors),
                 FUN = function(x){
                   ints <- int_and_neighbors_by_t_only_between[[t]]$neighbors[[x]]
                   clust_id <- clusters_accounting_for_noise[x]
                   if(!is.na(clust_id)){
                     
                     temp1 <- intersect(names(ints), names(which(clusters_accounting_for_noise == clust_id)))
                     ints[temp1] <- 0
                     p1 <-  ints[ints>0]
                     p2 <- sum(p1) 
                     
                    } else {
                      
                      p1 <-  ints
                      p2 <- sum(p1) 
                      
                    }
                     return(list(p1 = p1,
                                 p2 = p2))
                  
                 })
  names(temp) <- names(int_and_neighbors_by_t_only_between[[t]]$neighbors)
  
  int_and_neighbors_by_t_only_between[[t]]$neighbors<- lapply(temp,
                                                              FUN = function(x){x$p1})
  
  int_and_neighbors_by_t_only_between[[t]]$n_total_t <- sapply(temp,
                                                              FUN = function(x){x$p2})
  
}


nice_fun <- function(network, cluster_labels){
  
  assign("int_and_neighbors_by_t", network, envir = .GlobalEnv) # will add as argument to run_single_ABM
  ## Run ABM with clustering accounting for noise
  set.seed(1234, kind = "L'Ecuyer-CMRG")
  run_ABM <- parallel::mclapply(1:n_trial,
                                FUN = function(x){
                                  
                                  if(cluster_labels == "lshm"){
                                    clusters <- clusters_accounting_for_noise
                                  } else if(cluster_labels == "random_perm") {
                                    clusters_random_lshm <- clusters_accounting_for_noise
                                    permute <- sample(names(clusters_random_lshm), size = length(clusters_random_lshm), replace = FALSE)
                                    names(clusters_random_lshm) <- permute
                                    clusters <- clusters_random_lshm
                                  } else if (cluster_labels == "random"){
                                    non_isolates <- as.character(1:total_actors)[ as.character(1:total_actors) %in% names(clusters_accounting_for_noise) ]
                                    clusters_random_lshm <- rmultinom(length(non_isolates), 1, prob = prop.table(table(clusters_accounting_for_noise)))
                                    clusters_random_lshm <- apply(clusters_random_lshm , 2, which.max)
                                    names(clusters_random_lshm) <- non_isolates
                                    clusters <- clusters_random_lshm
                                  } else {
                                    clusters <-NULL
                                  }
                                  
                                  run_single_ABM(p_infected = p_infected,
                                                 mean_exposure_days = mean_exposure_days,
                                                 mean_infected_days = mean_infected_days,
                                                 actor_labels = actor_labels,
                                                 min_degree_t1 = min_degree_t1,
                                                 timesteps = timesteps,
                                                 n_repeat = n_repeat,
                                                 p_asym = p_asym,
                                                 clusters = clusters,
                                                 quarantine_days = quarantine_days,
                                                 digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                 DCT_sensitivity = DCT_sensitivity,
                                                 DCT_specificity = DCT_specificity)},
                                mc.preschedule = TRUE,
                                mc.cores = cores)
  
  # remove the global variable to clean up
  rm(int_and_neighbors_by_t, envir = .GlobalEnv)
  
  res <- list(
    attack_rate = quantile(do.call(c, lapply(run_ABM, function(x) x$attack_rate)), probs = c(0.25, 0.5, 0.75)),
    quarantined_rate = quantile(do.call(c, lapply(run_ABM, function(x) x$quarantined_rate)), probs = c(0.25, 0.5, 0.75)),
    total_infections = quantile(do.call(c, lapply(run_ABM, function(x) sum(x$incidence))), probs = c(0.25, 0.5, 0.75)),
    total_contacts = quantile(do.call(c, lapply(run_ABM, function(x) sum(x$average_interactions_by_time * 692))), probs = c(0.25, 0.5, 0.75)),
    person_days_quarantined = quantile(do.call(c, lapply(run_ABM, function(x) sum(x$quarantined))), probs = c(0.25, 0.5, 0.75))
  )
  
  # cat("Attack Rate:", quantile((do.call("c", lapply(run_ABM, function(x){(x$attack_rate)} ))), probs = c(0.25, 0.5, 0.75)), "\n")
  # cat("Quarantine Rate:", quantile((do.call("c", lapply(run_ABM, function(x){(x$quarantined_rate)} ))), probs = c(0.25, 0.5, 0.75)), "\n")
  # cat("Total infections:", quantile((do.call("c", lapply(run_ABM, function(x){sum(x$incidence)} ))), probs = c(0.25, 0.5, 0.75)), "\n")
  # cat("Total contacts:", quantile((do.call("c", lapply(run_ABM, function(x){sum(x$average_interactions_by_time * 692)} ))), probs = c(0.25, 0.5, 0.75)), "\n")
  # cat("Person-days quarantined:", quantile((do.call("c", lapply(run_ABM, function(x){sum(x$quarantined)} ))), probs = c(0.25, 0.5, 0.75)) , "\n")
  return(res)
}


## Run ABM with LSHM clustering accounting for noise on observed network
p1 <- nice_fun(int_and_neighbors_by_t_true, "lshm")

## Run ABM with LSHM clustering accounting for noise on network with only within cluster interaction
p2 <- nice_fun(int_and_neighbors_by_t_only_within, "lshm")

## Run ABM with LSHM clustering accounting for noise on network with only between cluster interaction
p3 <- nice_fun(int_and_neighbors_by_t_only_between, "lshm")

## Run ABM with random clustering accounting for noise on observed network
p4 <- nice_fun(int_and_neighbors_by_t_true, "random_perm")

## Run ABM with random clustering accounting for noise on network with only within cluster interaction
p5 <- nice_fun(int_and_neighbors_by_t_only_within, "random_perm")

## Run ABM with random clustering accounting for noise on network with only between cluster interaction
p6 <- nice_fun(int_and_neighbors_by_t_only_between, "random_perm")






