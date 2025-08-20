library(tidyverse)
# set.seed(1234)
# prob of transmisson among infected
p_infected <- 0.07

# days exposed
exposure_days <- 0L

# days infected
mean_infected_days <- 12L

# Initilize actors
actors <- initialize_actors(user_values$name)

# permuted timestamps for intereactions
interactions_temp <- cns %>% 
  distinct(number_timestamp) %>%
  mutate(new_timestamp = (row_number() - 1) %/% 288 + 1) # each timestamp is 5 mins so 288 timestamps day

# interactions <- interactions_temp %>% 
#   distinct(new_timestamp) %>% 
#   mutate(order = order(runif(n()))) %>% 
#   arrange(order) %>% 
#   mutate(timestamp = 1:n()) %>% 
#   inner_join(interactions_temp, ., by = "new_timestamp") %>% 
#   inner_join(cns, ., by = "number_timestamp") %>% 
#   arrange(timestamp) %>% 
#   select(timestamp, number_timestamp, user_a, user_b) %>% 
#   group_by(timestamp) %>% distinct(user_a, user_b) %>% 
#   ungroup()

interactions <- interactions_temp %>% 
  inner_join(cns, ., by = "number_timestamp") %>%
  rename(timestamp = new_timestamp) %>% 
  count(timestamp, user_a, user_b)

# N time-steps
n_timesteps <- nrow(distinct(interactions, timestamp))

# set.seed(123)
# randomly seed infection in actor
seed_index <- sample.int(n = length(actors), size = 1)
actors[[as.character(seed_index)]]$state <- 2L # infected
actors[[as.character(seed_index)]]$duration_infected <- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
actors[[as.character(seed_index)]]$infected_by_seed <- NA
actors[[as.character(seed_index)]]$ever_infected <- T

# Store results
res <- matrix(NA, nrow = 4, ncol = n_timesteps*2)
incidence <- numeric(n_timesteps*2)

k <- 1
for (t_start in 1:2){
  for (t in 1:n_timesteps) {
    
    # States at time t
    state_t <- sapply(as.character(1:length(actors)), function(x){actors[[x]]$state})
    
    # Who is suspectible at time t
    sup_t <- which(state_t == 0)
    
    # Who is exposed at time t
    exp_t <- which(state_t == 1)
    
    # Who is infected at time t
    inf_t <- which(state_t == 2)
    
    if(length(c(sup_t,exp_t,inf_t)) == 0){
      break
    }
    
    # Count number of interactions by actor
    interactions_t <- interactions %>% 
      filter(timestamp == t)
    
    n_interactions <- tapply(X = c(interactions_t$n, interactions_t$n), 
                             INDEX = c(interactions_t$user_a, interactions_t$user_b),
                             FUN = sum)
    
    invisible(lapply(names(n_interactions),
                     function(x){actors[[x]]$n_contacts <- actors[[x]]$n_contacts + n_interactions[[x]]}))
    
    # Transmission
    if( length(sup_t) > 0 & length(inf_t) > 0 ){
      
      # Determine interactions with infected actors
      int_t_infected <- tibble(infect_id = names(inf_t)) %>% 
        mutate(temp = map(infect_id, ~{
          out <- interactions_t %>%
            filter(user_a == .x | user_b == .x )  %>% 
            mutate(user_a = ifelse(user_a == .x, user_b, user_a)) %>% 
            group_by(user_a) %>%
            summarise(n = sum(n), .groups = "drop") 
        }))
      
      # determine prob of transmisson
      int_t_infected_trans <- int_t_infected %>% 
        select(temp) %>% 
        unnest(temp) %>% 
        group_by(user_a) %>% 
        summarise(n = sum(n), .groups = "drop") %>% 
        filter(user_a %in% sup_t)
      
      n_interactions_w_inf <- int_t_infected_trans$user_a
      transmisson <- n_interactions_w_inf[runif(length(n_interactions_w_inf)) < 1 - (1-p_infected)^(int_t_infected_trans$n)]
      incidence[k] <- length(transmisson)
      
      if(length(transmisson) > 0){
        # invisible(lapply(as.character(transmisson),
        #                  function(x){
        #                    actors[[x]]$state <- 1L
        #                    actors[[x]]$duration_exposed <- exposure_days
        #                  }))
        
        infected_by <- int_t_infected %>% 
          unnest(temp) %>% 
          filter(user_a %in% transmisson) %>% 
          group_by(user_a) %>% 
          mutate(p = 1-(1-p_infected)^(n),
                 p = p/sum(p),
                 chosen = sample(infect_id, size = 1, prob = p)) %>% 
          ungroup() %>% 
          distinct(user_a, infect_id = chosen)
        
        
        invisible(lapply(as.character(transmisson),
                         function(x){
                           actors[[x]]$state <- 2L
                           actors[[x]]$duration_infected<- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
                           actors[[x]]$ever_infected <- T
                           actors[[x]]$infected_by <- infected_by %>% filter(user_a == x) %>% .$infect_id
                         }))
      }
      
    }
    
    # update days exposed
    if(length(exp_t) > 0){
      
      invisible(lapply(names(exp_t),
                       function(x){actors[[x]]$duration_exposed <- actors[[x]]$duration_exposed - 1L}))
      
      dur_exposed_t <- sapply(names(exp_t),
                              function(x){actors[[x]]$duration_exposed})
      
      if(any(dur_exposed_t ==0)){
        
        invisible(lapply(names(exp_t)[which(dur_exposed_t ==0)],
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
    
    res[, k] <- c(length(sup_t), length(exp_t), length(inf_t),  length(actors)- length(c(sup_t,exp_t,inf_t)))
    k <- k + 1
    
  }
}

plot(x = 1:ncol(res), y = res[1,], type = "l", col = "red")
lines(x = 1:ncol(res), y = res[3,], col = "black")
lines(x = 1:ncol(res), y = res[4,], col = "blue")

max(incidence)

mean(sapply(names(actors),
            function(x){actors[[x]]$ever_infected}))

sum(sapply(setdiff(names(actors), seed_index),
           function(x){actors[[x]]$infected_by})== as.character(seed_index), na.rm = T)

ever_infected
