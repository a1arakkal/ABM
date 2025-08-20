
# Load in interactions by time list --------------------------------------------
setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")

# Helper function --------------------------------------------------------------

# This function is to initialze the environment to keep track of the actors's states
# through each step of the ABM

initialize_actors <- function(actor_names){
  
  actor_states <- list(state = 0L, # 0L susceptible, 1L exposed, 2L infected, 3L recoverd
                       first_gen_infection = F,
                       n_contacts = 0,
                       ever_infected = F,
                       duration_exposed = 0L,
                       duration_infected = 0L)
  
  actors <- lapply(1:length(actor_names), function(x){actor_states})
  names(actors) <- as.character(actor_names)
  actors <- list2env(actors)
  
  return(actors)
  
}

# Function for single run of the ABM
run_single_ABM <- function(p_infected, mean_exposure_days, 
                           mean_infected_days, actor_labels,
                           timesteps, n_repeat){
  
  # initialize actors
  actors <- initialize_actors(actor_labels)
  
  # Randomly select patient 0
  seed_index <- as.character(sample.int(n = length(actors), size = 1))
  actors[[seed_index]]$state <- 2L 
  actors[[seed_index]]$duration_infected <- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
  actors[[seed_index]]$ever_infected <- T
  
  # Pre-allocate results storage
  n_timesteps <- length(timesteps)
  res <- matrix(NA, nrow = 4, ncol = n_timesteps * n_repeat) # we will just repeat the observed interactions 2 week 4 times 
  rownames(res) <- c("Susceptible", "Exposed", "Infected", "Recovered")
  incidence <- numeric(n_timesteps * n_repeat) # keep track of new infections by day
  
  # initialize counter
  k <- 1
  
  # For loop number of time to repeat weeks
  for(R in 1:n_repeat){
    
    # For loop number for cycling through days
    for (i in 1:n_timesteps){
      
      # data at time
      t <- timesteps[i]
      int_and_neighbors_t <- int_and_neighbors_by_t[[t]]
      
      # States at time t
      state_t <- sapply(as.character(1:length(actors)), function(x){actors[[x]]$state})
      
      # Who is suspectible at time t
      sup_t <- names(which(state_t == 0))
      
      # Who is exposed at time t
      exp_t <- names(which(state_t == 1))
      
      # Who is infected at time t
      inf_t <- names(which(state_t == 2))
      
      # If everyone recovered break
      if(length(c(sup_t,exp_t,inf_t)) == 0){
        break
      }
      
      # Store counts of SEIR
      res[, k] <- c(length(sup_t), length(exp_t), length(inf_t),  length(actors)- length(c(sup_t,exp_t,inf_t)))
      
      # Count number of interactions by actor
      n_interactions <- int_and_neighbors_t$n_total_t

      # # Check
      # interactions_t <- interactions[interactions$timestamp == as.numeric(t), , drop = F]
      #
      # n_interactions2 <- tapply(X = c(interactions_t$n, interactions_t$n),
      #                          INDEX = c(interactions_t$user_a, interactions_t$user_b),
      #                          FUN = sum)
      # invisible(lapply(names(n_interactions),
      #                  function(x){actors[[x]]$n_contacts <- actors[[x]]$n_contacts + n_interactions[[x]]}))
      # all.equal(n_interactions, setNames(as.vector(n_interactions2), names(n_interactions2)))
      
      # Transmission
      if( length(sup_t) > 0 & length(inf_t) > 0 ){
        names(int_and_neighbors_t)
        
        interaction_w_infect <- int_and_neighbors_t$neighbors[inf_t]
          
        interaction_w_infect <- lapply(interaction_w_infect,
                                       function(x){x[names(x) %in% sup_t]})
        
        # # check
        # interaction_w_infect2 <- ( (interactions_t$user_a %in% inf_t) & (interactions_t$user_b %in% sup_t) ) |  ( (interactions_t$user_b %in% inf_t) & (interactions_t$user_a %in% sup_t) )

        if(length(unlist(interaction_w_infect))>1){
          
          total_ints <- tapply(unlist(unname(interaction_w_infect)), names(unlist(unname(interaction_w_infect))), sum)
          prob_trans <- 1-(1-p_infected)^total_ints
          
          # # Check
          # temp <- interactions_t[interaction_w_infect2, , drop = F]
          # n_interactions_w_infected2 <- lapply(inf_t,
          #                                     FUN = function(x){
          #                                       temp2 <- temp[temp$user_a %in% x |  temp$user_b %in% x, , drop = F]
          #                                       temp2 <- tapply(X = c(temp2$n, temp2$n),
          #                                                       INDEX = c(temp2$user_a, temp2$user_b),
          #                                                       FUN = sum)
          #                                       return(temp2[setdiff(names(temp2), as.character(x))])})
          # 
          # total_ints2 <- tapply(unlist(unname(n_interactions_w_infected2)), names(unlist(unname(n_interactions_w_infected2))), sum)
          # prob_trans2 <- 1-(1-p_infected)^total_ints2
          # all.equal(prob_trans, prob_trans2)
          
          transmission <- runif(length(prob_trans)) < prob_trans
          
          # Compute number of new infections THINK!!!!!!!!!!
          incidence[k] <- sum(transmission)
          
          if(any(transmission)){
            
            first_gen_interaction <- interaction_w_infect[[seed_index]]
            
            invisible(lapply(names(which(transmission)),
                             
                             function(x){
                               
                               actors[[x]]$ever_infected <- T
                               
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
        
        invisible(lapply(exp_t,
                         function(x){actors[[x]]$duration_exposed <- actors[[x]]$duration_exposed - 1L}))
        
        dur_exposed_t <- sapply(exp_t,
                                function(x){actors[[x]]$duration_exposed})
        
        if(any(dur_exposed_t == 0)){
          
          invisible(lapply(exp_t[which(dur_exposed_t == 0)],
                           function(x){
                             actors[[x]]$state <- 2L
                             actors[[x]]$duration_infected <- extraDistr::rtpois(n = 1, lambda = mean_infected_days, a = 0.0)
                           }))
          
        }  
        
      }
      
      # update days infected
      if(length(inf_t) > 0){
        
        invisible(lapply(inf_t,
                         function(x){actors[[x]]$duration_infected <- actors[[x]]$duration_infected - 1L}))
        
        dur_infect_t <- sapply(inf_t,
                               function(x){actors[[x]]$duration_infected})
        
        if(any(dur_infect_t == 0)){
          
          invisible(lapply(inf_t[which(dur_infect_t == 0)],
                           function(x){
                             actors[[x]]$state <- 3L
                           }))
          
        }  
        
      }
      
      # Increment
      k <- k + 1
      
    }
    
  }
  
  return(list(res = res,
         incidence = incidence,
         attack_rate =  mean(sapply(names(actors),
                                     function(x){actors[[x]]$ever_infected})),
         R0 = sum(sapply(names(actors),
                         function(x){actors[[x]]$first_gen_infection}), na.rm = T)))

}

# Parameters ------------------------------------------------------------------

# prob of transmisson among infected
p_infected <- 0.005

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- 5L

# Initilize actors
actor_labels <- 1:total_actors

# N time-steps
timesteps <- names(int_and_neighbors_by_t)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks


test <- run_single_ABM(p_infected = p_infected,
                       mean_exposure_days = mean_exposure_days,
                       mean_infected_days = mean_infected_days,
                       actor_labels = actor_labels,
                       timesteps = timesteps,
                       n_repeat = n_repeat)

# microbenchmark::microbenchmark(run_single_ABM(p_infected = p_infected,
#                                               mean_exposure_days = mean_exposure_days,
#                                               mean_infected_days = mean_infected_days,
#                                               actor_labels = actor_labels,
#                                               timesteps = timesteps,
#                                               n_repeat = n_repeat),
# times = 100)

# Run ABM for multiple trials
n_trial <- 1e4
run_ABM <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat)},
                              mc.preschedule = TRUE,
                              mc.cores = 5)

# Average R0 across trials
R0 <- sapply(run_ABM, function(x){x$R0})
mean(R0)
quantile(R0, probs = c(0.025,0.975))
boxplot(R0)

# Average acttack across trials
attack_rate <- sapply(run_ABM, function(x){x$attack_rate})
mean(attack_rate)
quantile(attack_rate, probs = c(0.025,0.975))
boxplot(attack_rate)

# Average incidence across trials
incidence <- Reduce("+", lapply(run_ABM, function(x){x$incidence}))/n_trial
max(incidence)
data.frame(incidence = incidence,
           time = 1:length(incidence)) %>% 
ggplot(aes(x = time, y = incidence)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    title = "Incidence averaged over trials",
    x = "Time",
    y = "Mean Incidence"
  ) 

# Average trajectories across trials
res_array <- simplify2array(lapply(run_ABM, function(x) x$res))
res_mean <- apply(res_array, c(1,2), mean)
res_quantile_lower <- apply(res_array, c(1,2), quantile, probs = 0.025)
res_quantile_upper<- apply(res_array, c(1,2), quantile, probs = 0.975)

time <- 1:ncol(res_mean)
compartments <- rownames(res_mean)

# Convert to long format for ggplot
df <- data.frame(
  time = rep(time, each = nrow(res_mean)),
  Compartment = rep(compartments, times = ncol(res_mean)),
  mean = as.vector(res_mean),
  lower = as.vector(res_quantile_lower),
  upper = as.vector(res_quantile_upper)
)

library(ggplot2)
ggplot(df, aes(x = time, y = mean, color = Compartment, fill = Compartment)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    title = "ABM trajectories averaged over trials",
    x = "Time",
    y = "Mean number of actors"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")
