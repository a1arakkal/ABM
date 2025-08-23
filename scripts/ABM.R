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
                          "clusters_ignore_noise",
                          "clusters_dichotomize_g_1")))

# Source helper functions ------------------------------------------------------

source("https://raw.githubusercontent.com/a1arakkal/ABM/refs/heads/master/scripts/helper_functions.R")

# Parameters ------------------------------------------------------------------

# prob of transmisson among infected
p_infected <- 0.005

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- c("asymptomatic" = 5L,
                        "pre_symptomatic" = 2L,
                        "symptomatic" = 3L)

# prob of being asymptomatic given infected
p_asym <- 0.2

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector eg 1 will isolate individuals not clusters
quarantine_days <- 7 # NULL no intervention

# N time-steps
timesteps <- names(int_and_neighbors_by_t)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks

# test <- run_single_ABM(p_infected = p_infected,
#                        mean_exposure_days = mean_exposure_days,
#                        mean_infected_days = mean_infected_days,
#                        actor_labels = actor_labels,
#                        timesteps = timesteps,
#                        n_repeat = n_repeat,
#                        p_asym = p_asym,
#                        clusters = clusters,
#                        quarantine_days = quarantine_days)

# microbenchmark::microbenchmark(run_single_ABM(p_infected = p_infected,
#                                               mean_exposure_days = mean_exposure_days,
#                                               mean_infected_days = mean_infected_days,
#                                               actor_labels = actor_labels,
#                                               timesteps = timesteps,
#                                               n_repeat = n_repeat,
#                                               p_asym = p_asym,
#                                               clusters = clusters,
#                                               quarantine_days = quarantine_days),
# times = 100)

# Run ABM for multiple trials --------------------------------------------------

n_trial <- 1e4
cores <- 100L

## Run ABM with no intervention 
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_no_intervention <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = NULL,
                                                               quarantine_days = quarantine_days)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with isolating individual actors not clusters
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_isolate_individuals <- parallel::mclapply(1:n_trial,
                                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                                               mean_exposure_days = mean_exposure_days, 
                                                                               mean_infected_days = mean_infected_days,
                                                                               actor_labels = actor_labels,
                                                                               timesteps = timesteps,
                                                                               n_repeat = n_repeat,
                                                                               p_asym = p_asym,
                                                                               clusters = 1,
                                                                               quarantine_days = quarantine_days)},
                                              mc.preschedule = TRUE,
                                              mc.cores = cores)

## Run ABM with clustering accounting for noise
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_accounting_for_noise <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_accounting_for_noise,
                                                               quarantine_days = quarantine_days)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with clustering ignoring for noise
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_ignore_noise <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_ignore_noise,
                                                               quarantine_days = quarantine_days)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with clustering dichotomize network with cutoff of 1
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_dichotomize_g_1 <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_dichotomize_g_1,
                                                               quarantine_days = quarantine_days)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

save(run_ABM_no_intervention, run_ABM_isolate_individuals,
     run_ABM_accounting_for_noise,
     run_ABM_ignore_noise, run_ABM_dichotomize_g_1,
     file = "test.run.RData")

# Load libraries ---------------------------------------------------------------

library(ggplot2)

# Load ABM res -----------------------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")
load("test.run.RData")

# run_ABM <- run_ABM_no_intervention
# run_ABM <- run_ABM_isolate_individuals
# run_ABM <- run_ABM_accounting_for_noise
# run_ABM <- run_ABM_ignore_noise
# run_ABM <- run_ABM_dichotomize_g_1

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

# Average quarantined across trials ----------------------------------------------

quarantined <- Reduce("+", lapply(run_ABM, function(x){x$quarantined}))/n_trial
cat("Estimated max number of quarantined actors:", max(quarantined))
data <- data.frame(quarantined = quarantined,
                   time = 1:length(quarantined))

ggplot(data, aes(x = time, y = quarantined)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    title = "Number of quarantined actors by time averaged over trials",
    x = "Time (days)",
    y = "Mean number of quarantined actors"
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
