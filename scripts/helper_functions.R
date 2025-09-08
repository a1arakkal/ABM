# Helper function --------------------------------------------------------------

# This function is to initialze the environment to keep track of the actors's states
# through each step of the ABM

initialize_actors <- function(actor_names){
  
  actor_states <- list(state = 0L, # 0L susceptible, 1L exposed, 2L infected, 3L recoverd
                       first_gen_infection = FALSE, # Infected by seed
                       n_contacts = 0, # track cumulative number of contacts 
                       ever_infected = FALSE, # Indicator if the agent ever became infected at any point 
                       duration_exposed = 0L, # If exposed tracks days of latency (non-infectious) period 
                       duration_infected = c(0L, 0L), # If infected tracks days of infectious period which is split into 
                                                      # two cases: 1) SYMPTOMATIC infections will have a pre-symtomatic period 
                                                      # follow by a symtomatic period or 2) ASYMPTOMATIC infections just have one
                                                      # period. The first element will be duration of the pre-symtomatic/asymtomatic period
                                                      # and the second element is the total length of the infectious period.
                                                      # So for ASYMPTOMATIC infections element 1 will equal element 2
                                                      # and for SYMPTOMATIC infections the symptomatic period is just element 2 minus element 1.
                       duration_quarantined = 0L, # if quarantined tracks days in quarantine
                       ever_quarantined = FALSE)  # Indicator if the agent was ever quarantined at any point 
  
  actors <- lapply(1:length(actor_names), function(x){actor_states})
  names(actors) <- as.character(actor_names)
  actors <- list2env(actors)
  
  return(actors)
  
}

# Function for single run of the ABM
run_single_ABM <- function(p_infected, mean_exposure_days, 
                           mean_infected_days, actor_labels,
                           timesteps, n_repeat, min_degree_t1,
                           p_asym, clusters, quarantine_days,
                           digital_contact_tracing_look_back){
  
  # initialize actors
  actors <- initialize_actors(actor_labels)
  
  # Randomly select patient 0
  # seed_index <- as.character(sample.int(n = length(actors), size = 1))
  # This approach will randomly sample from all actors 1-692. This may cause 
  # issues as there will be actors who will not interact at time point 1 causing
  # infection to die out quickly. Specifically about 52.5% dont have interactions
  # at 1st time point.
  # (length(int_and_neighbors_by_t[[1]]$n_total)/692)*100

  # we could instead select those with at least x degree at first time point to 
  # ensure propagation
  # sum(int_and_neighbors_by_t[[1]]$n_total >= min_degree_t1)
  
  if(min_degree_t1 == 0){
    seed_options <- as.character(1:length(actors))
  } else {
    seed_options <- names(which(int_and_neighbors_by_t[[1]]$n_total >= min_degree_t1))
  }

  seed_index <- seed_options[sample.int(n = length(seed_options), size = 1)]
  
  actors[[seed_index]]$state <- 2L 
  actors[[seed_index]]$ever_infected <- TRUE
  symptomatic_infection <- runif(1) > p_asym
  
  if(symptomatic_infection){
    
    durations <- c(ifelse(mean_infected_days[["pre_symptomatic"]] == 0, 0,
                          extraDistr::rtpois(n = 1, lambda = mean_infected_days[["pre_symptomatic"]], a = 0)),
                   extraDistr::rtpois(n = 1, lambda = mean_infected_days[["symptomatic"]], a = 0))
    
    actors[[seed_index]]$duration_infected <- c(durations[1], sum(durations))
    
  } else{
    
    actors[[seed_index]]$duration_infected <- rep(extraDistr::rtpois(n = 1, lambda = mean_infected_days[["asymptomatic"]], a = 0), 2)
    
  }
 
  # Pre-allocate results storage
  n_timesteps <- length(timesteps)
  res <- matrix(NA, nrow = 4, ncol = n_timesteps * n_repeat) # we will just repeat the observed interactions 2 week 4 times 
  rownames(res) <- c("Susceptible", "Exposed", "Infective", "Recovered")
  incidence <- numeric(n_timesteps * n_repeat) # keep track of new infections by day
  average_interaction <- numeric(n_timesteps * n_repeat) 
  quarantined_t <- numeric(n_timesteps * n_repeat)
  
  # initialize counter
  k <- 1
  
  # make copy of int_and_neighbors_by_t will be of length n_timesteps * n_repeat
  copy_int_and_neighbors_by_t <- rep(int_and_neighbors_by_t, n_repeat)
  
  # For loop number of time to repeat weeks
  for(R in 1:n_repeat){
    
    # For loop number for cycling through days
    for (i in 1:n_timesteps){
      
      # data at time
      t <- timesteps[i]
      int_and_neighbors_t <- int_and_neighbors_by_t[[t]]
      
      # States at time t
      state_t <- sapply(actors, function(x){x$state})
      
      # Who is suspectible at time t
      sup_t <- names(which(state_t == 0))
      
      # Who is exposed at time t
      exp_t <- names(which(state_t == 1))
      
      # Who is infective at time t
      inf_t <- names(which(state_t == 2))
      
      # If everyone recovered break
      if(length(c(sup_t,exp_t,inf_t)) == 0){
        break
      }
      
      # Who is quarantined at time t
      qua_t <- names(which(sapply(actors, function(x){x$duration_quarantined > 0})))
      quarantined_t[k] <- length(qua_t)
      
      # remove interations for those quarantined
      if(length(qua_t) > 0){
        int_and_neighbors_t$neighbors <- int_and_neighbors_t$neighbors[setdiff(names(int_and_neighbors_t$neighbors), qua_t)] 
        # test <- lapply(int_and_neighbors_t$neighbors,
        #                                         FUN = function(x){
        #                                           x[setdiff(names(x), qua_t)]
        #                                         })
        
        if(length(int_and_neighbors_t$neighbors) > 0){
          temp <- lapply(int_and_neighbors_t$neighbors,
                         FUN = function(x){
                           p1 <- x[setdiff(names(x), qua_t)] # remove the interactions with quarantined actors
                           p2 <- sum(p1) # total interactions excluding quarantined actors
                           return(list(p1 = p1,
                                       p2 = p2))
                         })
          
          int_and_neighbors_t$neighbors <- lapply(temp,
                                                  FUN = function(x){x$p1})
          
          int_and_neighbors_t$n_total_t <- sapply(temp,
                                                  FUN = function(x){x$p2})
        }
        
        copy_int_and_neighbors_by_t[[k]]$neighbors <- int_and_neighbors_t$neighbors
        copy_int_and_neighbors_by_t[[k]]$n_total_t <- int_and_neighbors_t$n_total_t
        
      }
      
      # Store counts of SEIR
      res[, k] <- c(length(sup_t), length(exp_t), length(inf_t), length(actors) - length(c(sup_t, exp_t, inf_t)))
      
      if(length(int_and_neighbors_t$neighbors) > 0){
        # Count number of interactions by actor
        n_interactions <- int_and_neighbors_t$n_total_t
        average_interaction[k] <- sum(n_interactions)/length(actors) # not just mean(n_interactions) as we want to include the 0's
        
        # Update cumulative number of interactions per actor
        invisible(lapply(names(n_interactions),
                         function(x){actors[[x]]$n_contacts <- actors[[x]]$n_contacts + n_interactions[[x]]}))
      }
      
      # Transmission
      if( length(sup_t) > 0 & length(inf_t) > 0 ){
        
        # Check if infective actors have interactions at time t
        inft_had_interaction_t <- inf_t %in% names(int_and_neighbors_t$neighbors)
        
        if(any(inft_had_interaction_t)){
          
          interaction_w_infect <- int_and_neighbors_t$neighbors[inf_t[inft_had_interaction_t]]
          
          interaction_w_infect <- lapply(interaction_w_infect,
                                         function(x){x[names(x) %in% sup_t]})
          
          if(length(unlist(interaction_w_infect)) > 0){
            
            total_ints <- tapply(unlist(unname(interaction_w_infect)), names(unlist(unname(interaction_w_infect))), sum)
            prob_trans <- 1-(1-p_infected)^total_ints
            transmission <- runif(length(prob_trans)) < prob_trans
            
            # Compute number of new infections 
            incidence[k] <- sum(transmission)
            
            if(any(transmission)){
              
              first_gen_interaction <- interaction_w_infect[[seed_index]]
              
              invisible(lapply(names(which(transmission)),
                               
                               function(x){
                                 
                                 actors[[x]]$ever_infected <- TRUE
                                 
                                 if( !is.null(names(first_gen_interaction)) && (x %in% names(first_gen_interaction)) ){
                                   
                                   probs <- sapply(interaction_w_infect,
                                                   FUN = function(y){
                                                     if(x %in% names(y)){
                                                       temp <- y[[x]]
                                                       ( (1 - (1-p_infected)^temp) * (1-p_infected)^(total_ints[[x]]-temp) )
                                                     } else {
                                                       0
                                                     }
                                                   })
                                   
                                   probs <- probs/sum(probs)
                                   
                                   if(any(is.nan(probs))) {
                                     probs <- rep(1/length(probs), length(probs))
                                   }
                                   
                                   actors[[x]]$first_gen_infection <- runif(1) < probs[[seed_index]]
                                   
                                 }
                                 
                                 if(is.null(mean_exposure_days)){
                                   
                                   actors[[x]]$state <- 2L
                                   symptomatic_inf <- runif(1) > p_asym
                                   
                                   if(symptomatic_inf){
                                     
                                     durations <- c(ifelse(mean_infected_days[["pre_symptomatic"]] == 0, 0,
                                                           extraDistr::rtpois(n = 1, lambda = mean_infected_days[["pre_symptomatic"]], a = 0)),
                                                    extraDistr::rtpois(n = 1, lambda = mean_infected_days[["symptomatic"]], a = 0))
                                     
                                     actors[[x]]$duration_infected <- c(durations[1], sum(durations))
                                     
                                   } else{
                                     
                                     actors[[x]]$duration_infected <- rep(extraDistr::rtpois(n = 1, lambda = mean_infected_days[["asymptomatic"]], a = 0), 2)
                                     
                                   }
                                   
                                 } else {
                                   
                                   actors[[x]]$state <- 1L
                                   actors[[x]]$duration_exposed <- extraDistr::rtpois(n = 1, lambda = mean_exposure_days, a = 0)
                                   
                                 }
                                 
                               }))
              
            }
            
          }
          
        }
        
      }
      
      # update days exposed
      if(length(exp_t) > 0){
        
        invisible(lapply(exp_t,
                         function(x){actors[[x]]$duration_exposed <- actors[[x]]$duration_exposed - 1L}))
        
        dur_exposed_t <- sapply(exp_t,
                                function(x){actors[[x]]$duration_exposed})
        
        if(any(dur_exposed_t == 0)){
          
          invisible(lapply(exp_t[which(dur_exposed_t == 0)],
                           function(x){
                             
                             actors[[x]]$state <- 2L # infective
                             symptomatic_inf <- runif(1) > p_asym
                             
                             if(symptomatic_inf){
                               
                               durations <- c(ifelse(mean_infected_days[["pre_symptomatic"]] == 0, 0,
                                                     extraDistr::rtpois(n = 1, lambda = mean_infected_days[["pre_symptomatic"]], a = 0)),
                                              extraDistr::rtpois(n = 1, lambda = mean_infected_days[["symptomatic"]], a = 0))
                               
                               actors[[x]]$duration_infected <- c(durations[1], sum(durations))
                               
                             } else{
                               
                               actors[[x]]$duration_infected <- rep(extraDistr::rtpois(n = 1, lambda = mean_infected_days[["asymptomatic"]], a = 0), 2)
                               
                             }
                           }))
          
        }  
        
      }
      
      # update days infective
      if(length(inf_t) > 0){
        
        invisible(lapply(inf_t,
                         function(x){actors[[x]]$duration_infected <- actors[[x]]$duration_infected - 1L
                                     actors[[x]]$duration_infected[1] <- max(actors[[x]]$duration_infected[1], 0L)}))
        
        dur_infect_t <- sapply(inf_t,
                               function(x){actors[[x]]$duration_infected[2]})
        
        if(any(dur_infect_t == 0)){
          
          invisible(lapply(inf_t[which(dur_infect_t == 0)],
                           function(x){
                             actors[[x]]$state <- 3L # recovered
                           }))
          
        }  
        
      }
      
      # update days quarantined
      if(length(qua_t) > 0){
        
        invisible(lapply(qua_t,
                         function(x){actors[[x]]$duration_quarantined <- actors[[x]]$duration_quarantined - 1L}))
    
      }
      
      # Implement intervention
      potential_qua_t <- setdiff(inf_t, qua_t) #check if infected is already quarantined
      
      if(!is.null(clusters) && length(potential_qua_t) > 0){ 
        
        # Find those who have >0 symptomatic infectious days and done with asymtomatic period
        dur_symptom_infect_t <- sapply(potential_qua_t,
                                       function(x){ (abs(diff(actors[[x]]$duration_infected)) > 0) && actors[[x]]$duration_infected[1] == 0})
        
        symptom_infect_t <- potential_qua_t[unname(dur_symptom_infect_t)]
        
        ####### need to add check if infected is already quarantined ########
          
        if(length(symptom_infect_t) > 0){
          
          if(length(clusters) > 1){
          
          # Find those in the same cluster as symtomatic actors
          
          ## Find those that are in a cluster
          symptom_infect_t_in_cluster <- symptom_infect_t[symptom_infect_t %in% names(clusters)]
          
          ## Find those that not in a cluster as they are isolates removed before running LSHM
          symptom_infect_t_not_in_cluster <- setdiff(symptom_infect_t, symptom_infect_t_in_cluster)
          
          ## Get the cluster lables for those in clusters
          cluster_infected <- unique(unname(clusters[symptom_infect_t_in_cluster]))
          
          ## Get all the actors in those clusters
          quarantine_in_cluster <- names(clusters[clusters %in% cluster_infected])
          
          ## Quarantine actors include both those in an infected cluster and isolates
          quarantine <- unique(c(quarantine_in_cluster,
                                 symptom_infect_t_not_in_cluster)) # deals with isolates that are removed from LSHM implementation
          
          } else {
            
            # Find 1 hop neighbors of symptom_infect_t during [max(t-digital_contact_tracing_look_back, 1), t]
            index_current <- k
            one_hop <- unique(unlist(lapply(copy_int_and_neighbors_by_t[max(index_current-digital_contact_tracing_look_back, 1):index_current],
                                            FUN = function(x){unique(names(unlist(unname(x$neighbors[symptom_infect_t]))))}),
                                     use.names = FALSE))
            
            ## Quarantine actors include those infected and their one-hop neighbors
            quarantine <- unique(c(one_hop,
                                   symptom_infect_t)) 
            
          }
          
          # Exclude those already in quarantine
          quarantine <- setdiff(quarantine, qua_t)
          
          if (length(quarantine) > 0){
            
            invisible(lapply(quarantine, function(x){actors[[x]]$duration_quarantined <- quarantine_days
                                                     actors[[x]]$ever_quarantined <- TRUE}))
            
          }
          
        }  
      }
      
      # Increment
      k <- k + 1
      
    }
    
  }
  
  return(list(res = res,
              incidence = incidence,
              quarantined = quarantined_t,
              average_interactions_by_time = average_interaction,
              average_cumulative_interactions_per_actor = mean(sapply(names(actors),
                                                            function(x){actors[[x]]$n_contacts})),
              attack_rate = mean(sapply(names(actors),
                                        function(x){actors[[x]]$ever_infected})),
              quarantined_rate = mean(sapply(names(actors),
                                        function(x){actors[[x]]$ever_quarantined})),
              R0 = sum(sapply(names(actors),
                              function(x){actors[[x]]$first_gen_infection}), na.rm = T)))
  
}
