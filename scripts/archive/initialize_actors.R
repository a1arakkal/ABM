
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


       