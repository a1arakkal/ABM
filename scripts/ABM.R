# Load libraries ---------------------------------------------------------------

library(ggplot2)

# Load in interactions by time list --------------------------------------------

setwd("/Users/atlan/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")

# Source helper functions ------------------------------------------------------

source("https://raw.githubusercontent.com/a1arakkal/ABM/refs/heads/master/scripts/helper_functions.R")

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

# test <- run_single_ABM(p_infected = p_infected,
#                        mean_exposure_days = mean_exposure_days,
#                        mean_infected_days = mean_infected_days,
#                        actor_labels = actor_labels,
#                        timesteps = timesteps,
#                        n_repeat = n_repeat)

# microbenchmark::microbenchmark(run_single_ABM(p_infected = p_infected,
#                                               mean_exposure_days = mean_exposure_days,
#                                               mean_infected_days = mean_infected_days,
#                                               actor_labels = actor_labels,
#                                               timesteps = timesteps,
#                                               n_repeat = n_repeat),
# times = 100)

# Run ABM for multiple trials --------------------------------------------------

n_trial <- 1e4
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat)},
                              mc.preschedule = TRUE,
                              mc.cores = 40)

save(run_ABM, file = "test.run.RData")

# Load ABM res -----------------------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")
load("test.run.RData")
n_trial <- length(run_ABM)

# Average R0 across trials -----------------------------------------------------

R0 <- sapply(run_ABM, function(x){x$R0})
cat("Estimated R0:", mean(R0))
quantile(R0, probs = c(0.025,0.975))
boxplot(R0)

# Introducing exposure days can increase R0 as the seed will have
# to compete with fewer infectious agents at a given t thus more likely to have 
# 1st gen infections. But in the limit of exposure days it should give most 
# accurate estimate of R0 as only infectious agent will be seed and all
# infections will be 1st gen.

# Average acttack across trials ------------------------------------------------

attack_rate <- sapply(run_ABM, function(x){x$attack_rate})
mean(attack_rate)
quantile(attack_rate, probs = c(0.025,0.975))
boxplot(attack_rate)

# Average incidence across trials ----------------------------------------------

incidence <- Reduce("+", lapply(run_ABM, function(x){x$incidence}))/n_trial
cat("Estimated max number of incident cases:", max(incidence))
data <- data.frame(incidence = incidence,
           time = 1:length(incidence))

ggplot(data, aes(x = time, y = incidence)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    title = "Number of incident cases by time averaged over trials",
    x = "Time (days)",
    y = "Mean Incidence"
  ) 

# Average interactions per time across trials ----------------------------------

average_interactions_by_time <- Reduce("+", lapply(run_ABM, function(x){x$average_interactions_by_time}))/n_trial
data <- data.frame(incidence = average_interactions_by_time,
                   time = 1:length(average_interactions_by_time))

ggplot(data, aes(x = time, y = average_interactions_by_time)) +
  geom_line(size = 1) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(7, 56, 7),
                     labels = c(1:8)) +
  labs(
    title = "Average number of interactions per actor by time averaged over trials",
    x = "Time (weeks)",
    y = "Mean number of interactions"
  ) 

# Check
# interactions %>%
#   filter(timestamp>14) %>%
#   group_by(timestamp) %>%
#   summarise(ave = 2*sum(n)/692,
#             n_int = sum(n),
#             den = n_int/((total_actors*(total_actors-1))*0.5),
#             ave2 = den*691) %>%
#   ggplot(aes(x = timestamp, y = ave)) +
#   geom_line(size = 1) +
#   theme_minimal() +
#   labs(
#     title = "Average number of interactions by time averaged over trials",
#     x = "Time (weeks)",
#     y = "Mean number of interactions"
#   )
# 
# interactions %>%
#   filter(timestamp>14) %>% 
#   nest(data = c(user_a, user_b, n), .by = timestamp) %>% 
#   mutate(ave = map_dbl(data,
#                         ~{select(.x, actor = user_a, n) %>% 
#                             bind_rows(select(.x, actor = user_b, n)) %>% 
#                             group_by(actor) %>% 
#                             summarise(total = sum(n)) %>% 
#                             ungroup() %>% 
#                             summarise(ave = sum(total)/692) %>% .$ave})) %>% 
#   ggplot(aes(x = timestamp, y = ave)) +
#   geom_line(size = 1) +
#   theme_minimal() +
#   labs(
#     title = "Average number of interactions by time averaged over trials",
#     x = "Time (weeks)",
#     y = "Mean number of interactions"
#   )
# 
# plot(sapply(int_and_neighbors_by_t, function(x) sum(x$n_total_t)/total_actors), type = "l")


# Average cumulative interaction across trials ---------------------------------

average_cumulative_interactions_per_actor <- Reduce("+", lapply(run_ABM, function(x){x$average_cumulative_interactions_per_actor}))/n_trial
cat("Estimated average cumulative number of interactions per actor:", average_cumulative_interactions_per_actor)

# # Check
# interactions %>%
#   filter(timestamp>14) %>% 
#   select(actor = user_a, n) %>% 
#   bind_rows(interactions %>%
#               filter(timestamp>14) %>% 
#               select(actor = user_b, n)) %>% 
#   group_by(actor) %>% 
#   summarise(total = sum(n)) %>% 
#   ungroup() %>% 
#   summarise(sum(total)*4/692) # times 4 as we duplicated the 2 weeks 4 times

# Average trajectories across trials -------------------------------------------

res_array <- simplify2array(lapply(run_ABM, function(x) x$res))
res_mean <- apply(res_array, c(1,2), mean)
res_quantile_lower <- apply(res_array, c(1,2), quantile, probs = 0.025)
res_quantile_upper<- apply(res_array, c(1,2), quantile, probs = 0.975)

time <- 1:ncol(res_mean)
compartments <- rownames(res_mean)

df <- data.frame(
  time = rep(time, each = nrow(res_mean)),
  Compartment = rep(compartments, times = ncol(res_mean)),
  mean = as.vector(res_mean),
  lower = as.vector(res_quantile_lower),
  upper = as.vector(res_quantile_upper)
)

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
