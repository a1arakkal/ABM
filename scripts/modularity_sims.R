
# Load in interactions by time list --------------------------------------------

setwd("/Users/atlan/dissertation/real_data_application/paper3/")

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

modularity <- function(network,
                       total_actors,
                       assumed_true_communities){
  
  out <-  sapply(names(network), 
                 function(x){
                   
                   A <- matrix(0, nrow = total_actors, ncol= total_actors)
                   colnames(A) <- 1:total_actors
                   rownames(A) <- 1:total_actors
                   dat <- network[[x]]$neighbors
                   
                   for(i in names(dat)){
                     A[i, names(dat[[i]])] <- 1
                   }
                   
                   sum_within <- 0
                   for (i in unique(assumed_true_communities)){
                     sum_within <- sum_within + sum(A[names(which(assumed_true_communities == i)),
                                                      names(which(assumed_true_communities == i))])
                   }
                   
                   sum_between <- sum(A) - sum_within
                   
                   return(c(sum_within, sum_between, sum_within/(sum(A))))
                 })
  
}

modularity_clusters_accounting_for_noise <- modularity(int_and_neighbors_by_t_true,
                                                       total_actors,
                                                       clusters_accounting_for_noise)

modularity_clusters_fb_network <- modularity(int_and_neighbors_by_t_true,
                                             total_actors,
                                             clusters_fb_network)
  
# Function to adjust modularity ------------------------------------------------

adjust_modularity <- function(network,
                              total_actors,
                              assumed_true_communities,
                              prob_within,
                              tol = 1e-2){
  
  new_network <- vector("list", length = length(network))
  names(new_network) <- names(network)
  modularity <- numeric(length = length(network))
  names(modularity) <- names(network)
    
  for (t in names(network)){
    
    # extract within and between edges -----------------------------------------
    A <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(A) <- 1:total_actors
    rownames(A) <- 1:total_actors
    dat <- network[[t]]$neighbors
    
    for(i in names(dat)){
      A[i, names(dat[[i]])] <- unname(dat[[i]])
    }
   
    # isSymmetric(A)
    indices_w_ints <- which(A>0, arr.ind = TRUE)
    indices_w_ints <- cbind(indices_w_ints, 0)
    colnames(indices_w_ints) <- c("i", "j", "within")
    indices_w_ints <- indices_w_ints[indices_w_ints[, "j"] > indices_w_ints[, "i"], ]
    indices_w_ints <- cbind(indices_w_ints, A[indices_w_ints[, 1:2]])
    colnames(indices_w_ints) <- c("i", "j", "within", "value")
    
    within_loc <- (as.character(indices_w_ints[,"i"]) %in% names(assumed_true_communities)) & 
      (as.character(indices_w_ints[,"j"]) %in% names(assumed_true_communities)) &
      (assumed_true_communities[as.character(indices_w_ints[,"i"])] ==  assumed_true_communities[as.character(indices_w_ints[,"j"])])
    within_loc[is.na(within_loc)] <- FALSE
    indices_w_ints[within_loc, "within"] <- 1
    
    obs_p_within <- prop.table(table(indices_w_ints[, "within"]))["1"]
    within_edges  <- indices_w_ints[indices_w_ints[, "within"] == 1, ]
    between_edges <-  indices_w_ints[indices_w_ints[, "within"] == 0, ]
    
    if(!is.null(prob_within)){
      
    if(prob_within == 1){
      
      out <- within_edges
      
    } else if (prob_within == 0){
      
      out <- between_edges
      
    } else if ( abs(prob_within - obs_p_within) < tol){
      
      out <- indices_w_ints
      
    } else {
      
      if (obs_p_within < prob_within) {
        
        # need to remove some between edges
        within <- within_edges
        N_between <- round((nrow(within) * (1 - prob_within)) / prob_within)
        between <- between_edges[order(runif(nrow(between_edges))), ]
        out <- rbind(within, between[1:N_between, ])
        
      } else {
        
        # need to remove some within edges
        between <- between_edges
        N_within <- round((prob_within * nrow(between)) / (1 - prob_within))
        within <- within_edges[order(runif(nrow(within_edges))), ]
        out <- rbind(within[1:N_within, ], between)
        
      }
      
    } 
      
    } else {
      
      out <- indices_w_ints
      
    }
             
    B <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(B) <- 1:total_actors
    rownames(B) <- 1:total_actors
    
    B[out[, 1:2]] <- out[, "value"]
    B[out[, 2:1]] <- out[, "value"]
    # isSymmetric(B)
    n_total_t <- rowSums(B)
    check_isolates <- n_total_t == 0
    
    sum_within <- 0
    for (i in unique(assumed_true_communities)){
      sum_within <- sum_within + sum((B>0)[names(which(assumed_true_communities == i)),
                                       names(which(assumed_true_communities == i))])
    }
    modularity[t] <- sum_within/sum(B>0)
    # modularity[t] <- prop.table(table(out[, "within"]))[["1"]]

    temp2 <- lapply(rownames(B)[!check_isolates], function(y){
      z <- B[y, ]
      z <- z[z!=0]
      return(z)
    })
    
    names(temp2) <- rownames(B)[!check_isolates]
    
    new_network[[t]] <- list(neighbors = temp2,
                             n_total_t = n_total_t[!check_isolates])

  }
  
  return(list(network = new_network,
              modularity = modularity))
  
}

# test
# plot_fun <- function(network,
#                      total_actors,
#                      assumed_true_communities,
#                      prob_within,
#                      tol = 1e-2,
#                      t){
#   
#   if (is.null(prob_within)){
#     
#     A <- matrix(0, nrow = total_actors, ncol= total_actors)
#     colnames(A) <- 1:total_actors
#     rownames(A) <- 1:total_actors
#     dat <- network[[t]]$neighbors
#     
#     for(i in names(dat)){
#       A[i, names(dat[[i]])] <- unname(dat[[i]])
#     }
#     
#   } else {
#     
#     temp <- adjust_modularity(network = network,
#                               total_actors = total_actors,
#                               assumed_true_communities = assumed_true_communities,
#                               prob_within = prob_within,
#                               tol = tol-2)
#     
#     A <- matrix(0, nrow = total_actors, ncol= total_actors)
#     colnames(A) <- 1:total_actors
#     rownames(A) <- 1:total_actors
#     dat <- temp$network[[t]]$neighbors
#     
#     for(i in names(dat)){
#       A[i, names(dat[[i]])] <- unname(dat[[i]])
#     }
#     
#   }
#   
#   image(((A>0)*1)[as.numeric(names(sort(assumed_true_communities))), as.numeric(names(sort(assumed_true_communities)))])
#   
# }
# 
# plot_fun(network = int_and_neighbors_by_t_true,
#          total_actors = total_actors,
#          assumed_true_communities = clusters_accounting_for_noise,
#          prob_within = NULL,
#          t = "19")
# 
# modularity_clusters_accounting_for_noise <- modularity(adjust_modularity(network = int_and_neighbors_by_t_true,
#                                                                          total_actors = total_actors,
#                                                                          assumed_true_communities = clusters_accounting_for_noise,
#                                                                          prob_within = 0.8,
#                                                                          tol = 1e-2)$network,
#                                                        total_actors,
#                                                        clusters_accounting_for_noise)
# 
# modularity_clusters_fb_network <- modularity(adjust_modularity(network = int_and_neighbors_by_t_true,
#                                                                total_actors = total_actors,
#                                                                assumed_true_communities = clusters_fb_network,
#                                                                prob_within = 1,
#                                                                tol = 1e-2)$network,
#                                              total_actors,
#                                              clusters_fb_network)

# Parameters ------------------------------------------------------------------

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- c("asymptomatic" = 10L,
                        "pre_symptomatic" = 5L, # if 0 not pre_symptomatic period
                        "symptomatic" = 5L)

# prob of being asymptomatic given infected
p_asym_seq <- round(c(0, 0.2, 0.4, 0.6, 0.8), 1) # if 0 all infections are symptomatic

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector of length 1 will isolate infected and their 1-hop neighbors
quarantine_days  <- 5 # ignored if cluster is NULL (i.e., no intervention)

# N time-steps
timesteps <- names(int_and_neighbors_by_t_true)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks

# Minimum degree for seed options if 0 will randomly select from 1-692
min_degree_t1 <- 0

# How far to look back for digital contact tracing if used as intervention, inclusive of current day
digital_contact_tracing_look_back <- 4

# Sensitivity and specificity of digital contact tracing if used as intervention (allows for false positives and false negatives)
DCT_sensitivity_seq  <- round(c(0, seq(.6, 1, by = .1)), 1)
DCT_specificity <- 1

# Number of trials and number of cores
n_trial <- 1e3
cores <- 100L

set.seed(1234, kind = "L'Ecuyer-CMRG")
seeds <- sample.int(1e6, size = n_trial, replace = FALSE)

for(DCT_sensitivity in DCT_sensitivity_seq){
  for(p_asym in p_asym_seq){
    for (c_labels in c("clusters_accounting_for_noise", "clusters_fb_network")){
      for (mod in round(seq(0, 1, 0.1), 1)){
        
        network <- int_and_neighbors_by_t_true
        
        if(c_labels == "clusters_accounting_for_noise"){
          assumed_true_communities <- clusters_accounting_for_noise
        } else {
          assumed_true_communities <- clusters_fb_network
        }
        
        R0_supply <- 5.3
        
        p_infected <- 0.0025
        n_while_loop <- 0
        change_p <- 0.0005
        min_R0 <- R0_supply*0.95
        max_R0 <- R0_supply*1.05
        R0 <- 0
        
        net_data <- adjust_modularity(network = network,
                                      total_actors = total_actors,
                                      assumed_true_communities = assumed_true_communities,
                                      prob_within = mod)
        
        while(! (R0 >= min_R0 & R0 <= max_R0) ){
          
          if(n_while_loop>100){
            break
          }
          
          run_ABM_for_R0 <- parallel::mclapply(1:n_trial,
                                               FUN = function(x){
                                                 
                                                 set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                 
                                                 run_single_ABM(p_infected = p_infected,
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
                                                                DCT_specificity = DCT_specificity,
                                                                int_and_neighbors_by_t = net_data$network)},
                                               mc.preschedule = TRUE,
                                               mc.cores = cores)
          
          R0 <- mean(sapply(run_ABM_for_R0, function(x){x$R0}))
          
          if (R0 > max_R0)  {
            p_infected <- p_infected - change_p
          }
          
          if (R0 < min_R0)  {
            p_infected <- p_infected + change_p
          }
          
          n_while_loop <- n_while_loop + 1
          
        }
        
        rm(net_data)
        
        ## Run ABM with clustering accounting for noise
        run_ABM_no_intervention <- parallel::mclapply(1:n_trial,
                                      FUN = function(x){
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        net_data <- adjust_modularity(network = network,
                                                                      total_actors = total_actors,
                                                                      assumed_true_communities = assumed_true_communities,
                                                                      prob_within = mod)
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        run_single_ABM(p_infected = p_infected,
                                                       mean_exposure_days = mean_exposure_days,
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
                                                       DCT_specificity = DCT_specificity,
                                                       int_and_neighbors_by_t = net_data$network)
                                      },
                                      mc.preschedule = TRUE,
                                      mc.cores = cores)
        
        ## Run ABM with clustering accounting for noise
        run_ABM_clustering <- parallel::mclapply(1:n_trial,
                                      FUN = function(x){
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        net_data <- adjust_modularity(network = network,
                                                                      total_actors = total_actors,
                                                                      assumed_true_communities = assumed_true_communities,
                                                                      prob_within = mod)
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        run_single_ABM(p_infected = p_infected,
                                                       mean_exposure_days = mean_exposure_days,
                                                       mean_infected_days = mean_infected_days,
                                                       actor_labels = actor_labels,
                                                       min_degree_t1 = min_degree_t1,
                                                       timesteps = timesteps,
                                                       n_repeat = n_repeat,
                                                       p_asym = p_asym,
                                                       clusters = assumed_true_communities,
                                                       quarantine_days = quarantine_days,
                                                       digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                       DCT_sensitivity = DCT_sensitivity,
                                                       DCT_specificity = DCT_specificity,
                                                       int_and_neighbors_by_t = net_data$network)
                                      },
                                      mc.preschedule = TRUE,
                                      mc.cores = cores)
        
        # apply(do.call("rbind", run_ABM), 2, quantile, probs = c(0.25,0.5,0.75))
        
        ## Run ABM with random clusters similar size to clustering accounting for noise (allow different assignment for each run)
        run_ABM_random_clustering <- parallel::mclapply(1:n_trial,
                                             FUN = function(x){
                                               
                                               set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                               
                                               net_data <- adjust_modularity(network = network,
                                                                             total_actors = total_actors,
                                                                             assumed_true_communities = assumed_true_communities,
                                                                             prob_within = mod)
                                               
                                               #random clusters same size as lshm clusters
                                               clusters_random <- assumed_true_communities
                                               permute <- sample(names(clusters_random), size = length(clusters_random), replace = FALSE)
                                               names(clusters_random) <- permute
                                               
                                               set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                               
                                               run_single_ABM(p_infected = p_infected,
                                                              mean_exposure_days = mean_exposure_days,
                                                              mean_infected_days = mean_infected_days,
                                                              actor_labels = actor_labels,
                                                              min_degree_t1 = min_degree_t1,
                                                              timesteps = timesteps,
                                                              n_repeat = n_repeat,
                                                              p_asym = p_asym,
                                                              clusters = clusters_random,
                                                              quarantine_days = quarantine_days,
                                                              digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                              DCT_sensitivity = DCT_sensitivity,
                                                              DCT_specificity = DCT_specificity,
                                                              int_and_neighbors_by_t = net_data$network)},
                                             
                                             mc.preschedule = TRUE,
                                             mc.cores = cores)
        
        ## Run ABM with isolating individual actors not cluster
        run_ABM_isolate_individuals <- parallel::mclapply(1:n_trial,
                                      FUN = function(x){
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        net_data <- adjust_modularity(network = network,
                                                                      total_actors = total_actors,
                                                                      assumed_true_communities = assumed_true_communities,
                                                                      prob_within = mod)
                                        
                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                        
                                        run_single_ABM(p_infected = p_infected,
                                                       mean_exposure_days = mean_exposure_days,
                                                       mean_infected_days = mean_infected_days,
                                                       actor_labels = actor_labels,
                                                       min_degree_t1 = min_degree_t1,
                                                       timesteps = timesteps,
                                                       n_repeat = n_repeat,
                                                       p_asym = p_asym,
                                                       clusters = 1, # if vector of length 1 will use digital contact tracing approach
                                                       quarantine_days = quarantine_days,
                                                       digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                       DCT_sensitivity = DCT_sensitivity,
                                                       DCT_specificity = DCT_specificity,
                                                       int_and_neighbors_by_t = net_data$network)
                                      },
                                      mc.preschedule = TRUE,
                                      mc.cores = cores)
        
        
        # test <- run_ABM
        # p1 <- list(
        #   attack_rate = quantile(do.call("c", lapply(test, function(x) x$attack_rate)), probs = c(0.25, 0.5, 0.75)),
        #   quarantined_rate = quantile(do.call("c", lapply(test, function(x) x$quarantined_rate)), probs = c(0.25, 0.5, 0.75)),
        #   total_infections = quantile(do.call("c", lapply(test, function(x) sum(x$incidence))), probs = c(0.25, 0.5, 0.75)),
        #   total_contacts = quantile(do.call("c", lapply(test, function(x) sum(x$average_interactions_by_time * 692))), probs = c(0.25, 0.5, 0.75)),
        #   person_days_quarantined = quantile(do.call("c", lapply(test, function(x) sum(x$quarantined))), probs = c(0.25, 0.5, 0.75))
        # )
        # 
        # test <- run_ABM_random
        # p2 <- list(
        #   attack_rate = quantile(do.call("c", lapply(test, function(x) x$attack_rate)), probs = c(0.25, 0.5, 0.75)),
        #   quarantined_rate = quantile(do.call("c", lapply(test, function(x) x$quarantined_rate)), probs = c(0.25, 0.5, 0.75)),
        #   total_infections = quantile(do.call("c", lapply(test, function(x) sum(x$incidence))), probs = c(0.25, 0.5, 0.75)),
        #   total_contacts = quantile(do.call("c", lapply(test, function(x) sum(x$average_interactions_by_time * 692))), probs = c(0.25, 0.5, 0.75)),
        #   person_days_quarantined = quantile(do.call("c", lapply(test, function(x) sum(x$quarantined))), probs = c(0.25, 0.5, 0.75))
        # )
        
        res <- list(R0 = R0,
                    run_ABM_no_intervention = run_ABM_no_intervention,
                    run_ABM_clustering = run_ABM_clustering,
                    run_ABM_random_clustering = run_ABM_random_clustering,
                    run_ABM_isolate_individuals = run_ABM_isolate_individuals,
                    cluster = c_labels,
                    modularity = mod,
                    p_asym = p_asym,
                    DCT_sensitivity = DCT_sensitivity,
                    DCT_specificity = DCT_specificity)
        
        saveRDS(res,
                file = paste0("res_modularity/p_asym_", p_asym,"_modularity_", mod, "_", c_labels, "_DCT_sensitivity_", DCT_sensitivity, ".RDS"))
        
        
      }
    }
  }
}


