# prob of transmisson among infected
p_infected <- 0.01

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- 5L

# Initilize actors
actors <- initialize_actors(user_values$name)

# collapse 5 min timestamps for intereactions into daily level
interactions_temp <- cns %>% 
  distinct(number_timestamp) %>%
  mutate(new_timestamp = (row_number() - 1) %/% 288 + 1) # each timestamp is 5 mins so 288 timestamps day, subtract 1 as we want groups to be closed on right interval
# table(cut(1:8064, breaks = c(0, seq(288, 288*4*7, 288)) ))

# tibble with total actor pair interactions by day
interactions <- interactions_temp %>% 
  inner_join(cns, ., by = "number_timestamp") %>% 
  rename(timestamp = new_timestamp) %>% 
  count(timestamp, user_a, user_b)

interactions <- as.data.frame(interactions)

# N time-steps
n_timesteps <- nrow(distinct(interactions, timestamp)) 

set.seed(123)
# randomly seed infection in actor
seed_index <- as.character(sample.int(n = length(actors), size = 1))
actors[[seed_index]]$state <- 2L # infected
actors[[seed_index]]$duration_infected <- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
actors[[seed_index]]$ever_infected <- T

# Store results
res <- matrix(NA, nrow = 4, ncol = n_timesteps * 2) # we will just repeat the observed interactions
incidence <- numeric(n_timesteps * 2) # keep track of new infections by day

k <- 1
for(R in 1:2){
  
  for (t in 1:n_timesteps){
    
    # States at time t
    state_t <- sapply(as.character(1:length(actors)), function(x){actors[[x]]$state})
    
    # Who is suspectible at time t
    sup_t <- which(state_t == 0)
    
    # Who is exposed at time t
    exp_t <- which(state_t == 1)
    
    # Who is infected at time t
    inf_t <- which(state_t == 2)
    
    # If everyone recovered break
    if(length(c(sup_t,exp_t,inf_t)) == 0){
      break
    }
    
    # Store counts of SEIR
    res[, k] <- c(length(sup_t), length(exp_t), length(inf_t),  length(actors)- length(c(sup_t,exp_t,inf_t)))

    # Count number of interactions by actor
    interactions_t <- interactions[interactions$timestamp == t, , drop = F]
    
    n_interactions <- tapply(X = c(interactions_t$n, interactions_t$n), 
                             INDEX = c(interactions_t$user_a, interactions_t$user_b),
                             FUN = sum)
    
    invisible(lapply(names(n_interactions),
                     function(x){actors[[x]]$n_contacts <- actors[[x]]$n_contacts + n_interactions[[x]]}))
    
    # Transmission
    if( length(sup_t) > 0 & length(inf_t) > 0 ){
      
      interaction_w_infect <- ( (interactions_t$user_a %in% inf_t) & (interactions_t$user_b %in% sup_t) ) |  ( (interactions_t$user_b %in% inf_t) & (interactions_t$user_a %in% sup_t) )
      
      if(sum(interaction_w_infect)>1){
        
        temp <- interactions_t[interaction_w_infect, , drop = F]
        n_interactions_w_infected <- lapply(inf_t,
                                            FUN = function(x){
                                              temp2 <- temp[temp$user_a %in% x |  temp$user_b %in% x, , drop = F]
                                              temp2 <- tapply(X = c(temp2$n, temp2$n), 
                                                     INDEX = c(temp2$user_a, temp2$user_b),
                                                     FUN = sum)
                                              return(temp2[setdiff(names(temp2), as.character(x))])})
        
        total_ints <- tapply(unlist(unname(n_interactions_w_infected)), names(unlist(unname(n_interactions_w_infected))), sum)
        prob_trans <- 1-(1-p_infected)^total_ints
        transmission <- runif(length(prob_trans)) < prob_trans
        
        # Compute number of new infections THINK!!!!!!!!!!
        if(k == 1){
          incidence[k] <- 1
        } else {
          incidence[k] <- sum(transmission)
        }
      
        if(any(transmission)){
          
          first_gen_interaction <- n_interactions_w_infected[[seed_index]]
          
          invisible(lapply(names(which(transmission)),
                           
                           function(x){
                             
                             actors[[x]]$ever_infected <- T
                             
                             if( !is.null(names(first_gen_interaction)) && (x %in% names(first_gen_interaction)) ){
                  
                               probs <- sapply(n_interactions_w_infected,
                                               FUN = function(y){
                                                 if(x %in% names(y)){
                                                   temp <- y[[x]]
                                                   ( (1 - (1-p_infected)^temp) * (1-p_infected)^(total_ints[[x]]-temp) )
                                                 } else {
                                                   0
                                                 }
                                               })
                               
                               probs <- probs/sum(probs)
                               actors[[x]]$first_gen_infection <- runif(1) < probs[[seed_index]]
                               
                             }
                             
                             if(is.null(mean_exposure_days)){
                               
                               actors[[x]]$state <- 2L
                               actors[[x]]$duration_infected<- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
                               
                             } else {
                               
                               actors[[x]]$state <- 1L
                               actors[[x]]$duration_exposed <- extraDistr::rtpois(n = 1, lambda = mean_exposure_days, a = 0.0)
                               
                             }
                             
                           }))
          
        }
        
      }
      
    }
    
    # update days exposed
    if(length(exp_t) > 0){
      
      invisible(lapply(names(exp_t),
                       function(x){actors[[x]]$duration_exposed <- actors[[x]]$duration_exposed - 1L}))
      
      dur_exposed_t <- sapply(names(exp_t),
                              function(x){actors[[x]]$duration_exposed})
      
      if(any(dur_exposed_t == 0)){
        
        invisible(lapply(names(exp_t)[which(dur_exposed_t == 0)],
                         function(x){
                           actors[[x]]$state <- 2L
                           actors[[x]]$duration_infected <- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
                         }))
        
      }  
      
    }
    
    # update days infected
    if(length(inf_t) > 0){
      
      invisible(lapply(names(inf_t),
                       function(x){actors[[x]]$duration_infected <- actors[[x]]$duration_infected - 1L}))
      
      dur_infect_t <- sapply(names(inf_t),
                             function(x){actors[[x]]$duration_infected})
      
      if(any(dur_infect_t == 0)){
        
        invisible(lapply(names(inf_t)[which(dur_infect_t == 0)],
                         function(x){
                           actors[[x]]$state <- 3L
                         }))
        
      }  
      
    }
    
    # Increment
    k <- k + 1
    
  }
  
}



plot(x = 1:ncol(res), y = res[1,], type = "l", col = "red", ylim = c(0,692))
lines(x = 1:ncol(res), y = res[3,], col = "black")
lines(x = 1:ncol(res), y = res[4,], col = "blue")

max(incidence)

mean(sapply(names(actors),
            function(x){actors[[x]]$ever_infected}))

sum(sapply(names(actors),
            function(x){actors[[x]]$first_gen_infection}), na.rm = T)


