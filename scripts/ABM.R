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
mean_infected_days <- c("asymptomatic" = 10L,
                        "pre_symptomatic" = 5L, # if 0 not pre_symptomatic period
                        "symptomatic" = 5L)

# prob of being asymptomatic given infected
p_asym <- 0.6 # if 0 all infections are symptomatic

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector eg 1 will isolate individuals not clusters
quarantine_days <- 7 # ignored if cluster is NULL (i.e., no intervention)

# N time-steps
timesteps <- names(int_and_neighbors_by_t)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks

# Minimum degree for seed options if 0 will randomly select from 1-692
min_degree_t1 <- 0
  
# How far to look back for digital contact tracing if used as intervention inclusive of current day
digital_contact_tracing_look_back <- 4

# test <- run_single_ABM(p_infected = p_infected,
#                        mean_exposure_days = mean_exposure_days,
#                        mean_infected_days = mean_infected_days,
#                        actor_labels = actor_labels,
#                        min_degree_t1 = min_degree_t1,
#                        timesteps = timesteps,
#                        n_repeat = n_repeat,
#                        p_asym = p_asym,
#                        clusters = clusters,
#                        quarantine_days = quarantine_days)

# microbenchmark::microbenchmark(run_single_ABM(p_infected = p_infected,
#                                               mean_exposure_days = mean_exposure_days,
#                                               mean_infected_days = mean_infected_days,
#                                               actor_labels = actor_labels,
#                                               min_degree_t1 = min_degree_t1,
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
                                                               min_degree_t1 = min_degree_t1,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = NULL,
                                                               quarantine_days = quarantine_days,
                                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with isolating individual actors not clusters
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_isolate_individuals <- parallel::mclapply(1:n_trial,
                                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                                               mean_exposure_days = mean_exposure_days, 
                                                                               mean_infected_days = mean_infected_days,
                                                                               actor_labels = actor_labels,
                                                                               min_degree_t1 = min_degree_t1,
                                                                               timesteps = timesteps,
                                                                               n_repeat = n_repeat,
                                                                               p_asym = p_asym,
                                                                               clusters = 1, # if vector of length 1 will used digital contact tracing approach 
                                                                               quarantine_days = quarantine_days,
                                                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back)},
                                              mc.preschedule = TRUE,
                                              mc.cores = cores)

## Run ABM with clustering accounting for noise
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_accounting_for_noise <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               min_degree_t1 = min_degree_t1,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_accounting_for_noise,
                                                               quarantine_days = quarantine_days,
                                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with clustering ignoring for noise
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_ignore_noise <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               min_degree_t1 = min_degree_t1,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_ignore_noise,
                                                               quarantine_days = quarantine_days,
                                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

## Run ABM with clustering dichotomize network with cutoff of 1
set.seed(1234, kind = "L'Ecuyer-CMRG")
run_ABM_dichotomize_g_1 <- parallel::mclapply(1:n_trial,
                              FUN = function(x){run_single_ABM(p_infected = p_infected,
                                                               mean_exposure_days = mean_exposure_days, 
                                                               mean_infected_days = mean_infected_days,
                                                               actor_labels = actor_labels,
                                                               min_degree_t1 = min_degree_t1,
                                                               timesteps = timesteps,
                                                               n_repeat = n_repeat,
                                                               p_asym = p_asym,
                                                               clusters = clusters_dichotomize_g_1,
                                                               quarantine_days = quarantine_days,
                                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back)},
                              mc.preschedule = TRUE,
                              mc.cores = cores)

save(run_ABM_no_intervention, run_ABM_isolate_individuals,
     run_ABM_accounting_for_noise,
     run_ABM_ignore_noise, run_ABM_dichotomize_g_1,
     file = "test.run.RData")

# Load libraries ---------------------------------------------------------------

library(ggplot2)
library(dplyr)

# Load ABM res -----------------------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")
load("test.run.RData")

extract_function <- function(run_ABM){
  
  n_trial <- length(run_ABM)
  
  # Average R0 across trials -----------------------------------------------------
  
  R0 <- sapply(run_ABM, function(x){x$R0})
  # cat("Estimated R0:", mean(R0))
  # quantile(R0, probs = c(0.025,0.975))
  # boxplot(R0)
  
  # Introducing exposure days can increase R0 as the seed will have
  # to compete with fewer infectious agents at a given t thus more likely to have 
  # 1st gen infections. But in the limit of exposure days it should give most 
  # accurate estimate of R0 as only infectious agent will be seed and all
  # infections will be 1st gen.
  
  # Average attack across trials ------------------------------------------------
  
  attack_rate <- sapply(run_ABM, function(x){x$attack_rate})
  # mean(attack_rate)
  # quantile(attack_rate, probs = c(0.025,0.975))
  # boxplot(attack_rate)
  
  # Average quarantine across trials ------------------------------------------------
  
  quarantined_rate <- sapply(run_ABM, function(x){x$quarantined_rate})
  # mean(quarantined_rate)
  # quantile(quarantined_rate, probs = c(0.025,0.975))
  # boxplot(quarantined_rate)
  
  # Average incidence across trials ----------------------------------------------
  
  incidence <- lapply(seq_along(run_ABM), function(i) {
    data.frame(
      run = i,
      time = 1:length( run_ABM[[i]]$incidence),          
      cases = run_ABM[[i]]$incidence        
    )
  }) %>% bind_rows()
  
  ave_incidence <- Reduce("+", lapply(run_ABM, function(x){x$incidence}))/n_trial
  incidence_max <- do.call("c", lapply(run_ABM, function(x){max(x$incidence)}))
  # cat("Estimated max number of incident cases:", max(incidence))
  data <- data.frame(incidence = ave_incidence,
                     time = 1:length(ave_incidence))
  
  incident_cases_plot <- ggplot(data, aes(x = time, y = ave_incidene)) +
    geom_line(size = 1) +
    theme_minimal() +
    labs(
      title = "Number of incident cases by time averaged over trials",
      x = "Time (days)",
      y = "Mean Incidence"
    ) 
  
  # Average quarantined across trials ----------------------------------------------
  
  quarantined <- lapply(seq_along(run_ABM), function(i) {
    data.frame(
      run = i,
      time = 1:length( run_ABM[[i]]$quarantined),          
      cases = run_ABM[[i]]$quarantined        
    )
  }) %>% bind_rows()

  quarantined_ave <- Reduce("+", lapply(run_ABM, function(x){x$quarantined}))/n_trial
  quarantined_max <- do.call("c", lapply(run_ABM, function(x){max(x$quarantined)}))
  data <- data.frame(quarantined = quarantined_ave,
                     time = 1:length(quarantined_ave))
  
  quarantined_plot <- ggplot(data, aes(x = time, y = quarantined_ave)) +
    geom_line(size = 1) +
    theme_minimal() +
    labs(
      title = "Number of quarantined actors by time averaged over trials",
      x = "Time (days)",
      y = "Mean number of quarantined actors"
    ) 
  
  # Average interactions per time across trials ----------------------------------
  interactions <- lapply(seq_along(run_ABM), function(i) {
    data.frame(
      run = i,
      time = 1:length( run_ABM[[i]]$average_interactions_by_time),           
      cases = run_ABM[[i]]$average_interactions_by_time        
    )
  }) %>% bind_rows()
  
  average_interactions_by_time <- Reduce("+", lapply(run_ABM, function(x){x$average_interactions_by_time}))/n_trial
  data <- data.frame(average_interactions_by_time = average_interactions_by_time,
                     time = 1:length(average_interactions_by_time))
  
  average_interactions_plot <- ggplot(data, aes(x = time, y = average_interactions_by_time)) +
    geom_line(size = 1) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(3.5, 56-3.5, 7),
                       labels = paste0("Week ", 1:8)) +
    geom_vline(xintercept = seq(7, 56, 7), color = "black", lty = "dashed") +
    labs(
      title = "Average number of interactions per actor by time averaged over trials",
      x = "Time (Weeks)",
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
  
  average_cumulative_interactions_per_actor <- do.call("c", lapply(run_ABM, function(x){x$average_cumulative_interactions_per_actor}))
  # cat("Estimated average cumulative number of interactions per actor:", average_cumulative_interactions_per_actor)
  
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
  trajectories <- lapply(seq_along(run_ABM), function(i) {
    data.frame(
      run = i,
      time = rep(1:length( run_ABM[[i]]$res), each = 4)) %>% 
      cbind(as.data.frame(t( run_ABM[[i]]$res)))
  }) %>% bind_rows()
  
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
  
  ABM_plot <- ggplot(df, aes(x = time, y = mean, color = Compartment, fill = Compartment)) +
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
  
  return(structure(list(R0 = R0,
                        attack_rate = attack_rate,
                        incidence_max = incidence_max,
                        quarantined_rate = quarantined_rate,
                        quarantined_max = quarantined_max,
                        average_interactions_by_time = average_interactions_by_time,
                        average_cumulative_interactions_per_actor = average_cumulative_interactions_per_actor,
                        incidence = incidence,
                        quarantined = quarantined,
                        interactions_by_time = interactions,
                        trajectories = trajectories,
                        plots = list(ABM_plot = ABM_plot,
                                     average_interactions_plot = average_interactions_plot,
                                     quarantined_plot = quarantined_plot,
                                     incident_cases_plot = incident_cases_plot)), class = "ABM"))
  
}

print.ABM <- function(x,...){
  
  cat("Estimated R0:\n")
  cat("   Mean:", mean(x$R0), "\n")
  cat("   Median:", median(x$R0), "\n")
  cat("   (Q1, Q3):", quantile(x$R0, probs = c(0.25,0.75)), "\n\n")
  
  cat("Estimated attack rate:\n")
  cat("   Mean:", mean(x$attack_rate), "\n")
  cat("   Median:", median(x$attack_rate), "\n")
  cat("   (Q1, Q3):", quantile(x$attack_rate, probs = c(0.25,0.75)), "\n\n")
  
  cat("Estimated max number of daily incident cases:\n")
  cat("   Mean:", mean(x$incidence_max), "\n")
  cat("   Median:", median(x$incidence_max), "\n")
  cat("   (Q1, Q3):", quantile(x$incidence_max, probs = c(0.25,0.75)), "\n\n")
  
  cat("Estimated quarantined rate:\n")
  cat("   Mean:", mean(x$quarantined_rate), "\n")
  cat("   Median:", median(x$quarantined_rate), "\n")
  cat("   (Q1, Q3):", quantile(x$quarantined_rate, probs = c(0.25,0.75)), "\n\n")
  
  cat("Estimated max number of quarantined actors:\n")
  cat("   Mean:", mean(x$quarantined_max), "\n")
  cat("   Median:", median(x$quarantined_max), "\n")
  cat("   (Q1, Q3):", quantile(x$quarantined_max, probs = c(0.25,0.75)), "\n\n")
 
  cat("Estimated average cumulative number of interactions per actor over length of ABM:\n")
  cat("   Mean:", mean(x$average_cumulative_interactions_per_actor), "\n")
  cat("   Median:", median(x$average_cumulative_interactions_per_actor), "\n")
  cat("   (Q1, Q3):", quantile(x$average_cumulative_interactions_per_actor, probs = c(0.25,0.75)), "\n\n")
  
}

plot.ABM <- function(x){
  
  available_plots <- names(x$plots)
  
  # Prompt user to select a plot
  choice <- utils::menu(available_plots, title = "Select a plot to view:")
  
  if (choice == 0) {
    message("No plot selected.")
    return(invisible(NULL))
  }
  
  # Return the chosen plot
  return(x$plots[[available_plots[choice]]])
  
}

no_intervention <- extract_function(run_ABM_no_intervention)
isolate_individuals <- extract_function(run_ABM_isolate_individuals)
ignore_noise <- extract_function(run_ABM_ignore_noise)
accounting_for_noise <- extract_function(run_ABM_accounting_for_noise)
ABM_dichotomize_g_1 <- extract_function(run_ABM_dichotomize_g_1)

combined_list <- list("No Intervention" = no_intervention,
                      "Isolated Individuals" = isolate_individuals,
                      "CD - Ignore Noise" = ignore_noise,
                      "CD - Dichotomize Network" = ABM_dichotomize_g_1,
                      "CD - Accouting for Noise" = accounting_for_noise)

combined_fun <- function(metric){
  do.call("rbind",
          lapply(names(combined_list),
                 FUN = function(x){
                   if(!is.data.frame(combined_list[[x]][[metric]])){
                     df <- data.frame(group = x,
                                      metric = metric,
                                      value = combined_list[[x]][[metric]])
                   } else{
                     df <- combined_list[[x]][[metric]]
                     df$group <- x
                   }
                   return(df)
                 }))
}

# Plot for R0 by group ---------------------------------------------------------
combined_fun("R0") %>% 
  mutate(group = factor(group, levels = c("No Intervention",
                                          "Isolated Individuals",
                                          "CD - Ignore Noise",
                                          "CD - Dichotomize Network",
                                          "CD - Accouting for Noise"))) %>% 
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black") +
  theme_minimal() +
  labs(
    title = "Boxplot of R0 by group",
    x = "",
    y = "R0"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Plot for attack_rate by group ------------------------------------------------
combined_fun("attack_rate") %>% 
  mutate(group = factor(group, levels = c("No Intervention",
                                          "Isolated Individuals",
                                          "CD - Ignore Noise",
                                          "CD - Dichotomize Network",
                                          "CD - Accouting for Noise"))) %>% 
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black") +
  theme_minimal() +
  labs(
    title = "Boxplot of attack rate by group",
    x = "",
    y = "Attack Rate"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Plot for incidence by time and by group --------------------------------------
{
  df <- combined_fun("incidence") %>% 
    mutate(group = factor(group, levels = c("No Intervention",
                                            "Isolated Individuals",
                                            "CD - Ignore Noise",
                                            "CD - Dichotomize Network",
                                            "CD - Accouting for Noise"))) 
  
  # Compute median and IQR for each time and group
  summary_df <- df %>%
    group_by(time, group) %>%
    summarize(
      median = median(cases),
      q25 = quantile(cases, 0.25),
      q75 = quantile(cases, 0.75),
      iqr = IQR(cases),
      lower_whisker = max(min(cases), q25 - 1.5*IQR(cases)),
      upper_whisker = min(max(cases), q75 + 1.5*IQR(cases)),
      .groups = "drop"
    )
  
  # Define offsets for staggering groups
  group_levels <- levels(summary_df$group)
  n_groups <- length(group_levels)
  offset <- 0.17  # adjust spacing
  
  summary_df <- summary_df %>%
    mutate(time_jitter = time + (as.numeric(group) - (n_groups+1)/2)*offset)
  
  palette <- c(
    "No Intervention" = "#E41A1C",       # red
    "Isolated Individuals" = "#377EB8",  # blue
    "CD - Ignore Noise" = "#4DAF4A",     # green
    "CD - Dichotomize Network" = "#984EA3", # purple
    "CD - Accouting for Noise" = "#FF7F00"  # orange
  )
  
  ggplot(summary_df, aes(x = time_jitter)) +
    # Whiskers (thin black lines outside IQR)
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = lower_whisker, yend = q25), color = "black", size = 0.5) +
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = q75, yend = upper_whisker), color = "black", size = 0.5) +
    # IQR region (thicker colored line)
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = q25, yend = q75, color = group), size = 1) +
    # Median dot (same color as group)
    geom_point(aes(y = median, color = group), size = 1) +
    scale_color_manual(values = palette) +
    theme_minimal() +
    labs(
      title = "Number infected over time by group",
      x = "Time (Days)",
      y = "Number infected",
      color = NULL,
      fill = NULL
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    scale_x_continuous(breaks = seq(7, 56, 7)) 
}

# Plot for quarantined_rate by group -------------------------------------------
combined_fun("quarantined_rate") %>% 
  mutate(group = factor(group, levels = c("No Intervention",
                                          "Isolated Individuals",
                                          "CD - Ignore Noise",
                                          "CD - Dichotomize Network",
                                          "CD - Accouting for Noise"))) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Boxplot of quarantined rate by group",
    x = "",
    y = "% Ever Quarantined"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot for quarantined by time and by group -------------------------------------------
{
df <- combined_fun("quarantined") %>% 
  mutate(group = factor(group, levels = c("No Intervention",
                                          "Isolated Individuals",
                                          "CD - Ignore Noise",
                                          "CD - Dichotomize Network",
                                          "CD - Accouting for Noise"))) 

# Compute median and IQR for each time and group
summary_df <- df %>%
  group_by(time, group) %>%
  summarize(
    median = median(cases),
    q25 = quantile(cases, 0.25),
    q75 = quantile(cases, 0.75),
    iqr = IQR(cases),
    lower_whisker = max(min(cases), q25 - 1.5*IQR(cases)),
    upper_whisker = min(max(cases), q75 + 1.5*IQR(cases)),
    .groups = "drop"
  )

# Define offsets for staggering groups
group_levels <- levels(summary_df$group)
n_groups <- length(group_levels)
offset <- 0.17  # adjust spacing

summary_df <- summary_df %>%
  mutate(time_jitter = time + (as.numeric(group) - (n_groups+1)/2)*offset)

palette <- c(
  "No Intervention" = "#E41A1C",       # red
  "Isolated Individuals" = "#377EB8",  # blue
  "CD - Ignore Noise" = "#4DAF4A",     # green
  "CD - Dichotomize Network" = "#984EA3", # purple
  "CD - Accouting for Noise" = "#FF7F00"  # orange
)

ggplot(summary_df, aes(x = time_jitter)) +
  # Whiskers (thin black lines outside IQR)
  geom_segment(aes(x = time_jitter, xend = time_jitter, y = lower_whisker, yend = q25), color = "black", size = 0.5) +
  geom_segment(aes(x = time_jitter, xend = time_jitter, y = q75, yend = upper_whisker), color = "black", size = 0.5) +
  # IQR region (thicker colored line)
  geom_segment(aes(x = time_jitter, xend = time_jitter, y = q25, yend = q75, color = group), size = 1) +
  # Median dot (same color as group)
  geom_point(aes(y = median, color = group), size = 1) +
  scale_color_manual(values = palette) +
  theme_minimal() +
  labs(
    title = "Number quarantined over time by group",
    x = "Time (Days)",
    y = "Number quarantined",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  scale_x_continuous(breaks = seq(7, 56, 7)) 
}

# Plot for average number of interactions by time and by group -----------------
{
  df <- combined_fun("interactions_by_time") %>% 
    mutate(group = factor(group, levels = c("No Intervention",
                                            "Isolated Individuals",
                                            "CD - Ignore Noise",
                                            "CD - Dichotomize Network",
                                            "CD - Accouting for Noise"))) 
  
  # Compute median and IQR for each time and group
  summary_df <- df %>%
    group_by(time, group) %>%
    summarize(
      median = median(cases),
      q25 = quantile(cases, 0.25),
      q75 = quantile(cases, 0.75),
      iqr = IQR(cases),
      lower_whisker = max(min(cases), q25 - 1.5*IQR(cases)),
      upper_whisker = min(max(cases), q75 + 1.5*IQR(cases)),
      .groups = "drop"
    )
  
  # Define offsets for staggering groups
  group_levels <- levels(summary_df$group)
  n_groups <- length(group_levels)
  offset <- 0.17  # adjust spacing
  
  summary_df <- summary_df %>%
    mutate(time_jitter = time + (as.numeric(group) - (n_groups+1)/2)*offset)
  
  palette <- c(
    "No Intervention" = "#E41A1C",       # red
    "Isolated Individuals" = "#377EB8",  # blue
    "CD - Ignore Noise" = "#4DAF4A",     # green
    "CD - Dichotomize Network" = "#984EA3", # purple
    "CD - Accouting for Noise" = "#FF7F00"  # orange
  )
  
  ggplot(summary_df, aes(x = time_jitter)) +
    # Whiskers (thin black lines outside IQR)
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = lower_whisker, yend = q25), color = "black", size = 0.5) +
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = q75, yend = upper_whisker), color = "black", size = 0.5) +
    # IQR region (thicker colored line)
    geom_segment(aes(x = time_jitter, xend = time_jitter, y = q25, yend = q75, color = group), size = 1) +
    # Median dot (same color as group)
    geom_point(aes(y = median, color = group), size = 1) +
    scale_color_manual(values = palette) +
    theme_minimal() +
    labs(
      title = "Average number of interactions over time by group",
      x = "Time (Days)",
      y = "Average number of interactions",
      color = NULL,
      fill = NULL
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    scale_x_continuous(breaks = seq(7, 56, 7)) 
}

# Plot for average cumulative number of interactions across actors (i.e., average strength) by group -----------------
combined_fun("average_cumulative_interactions_per_actor") %>% 
  mutate(group = factor(group, levels = c("No Intervention",
                                          "Isolated Individuals",
                                          "CD - Ignore Noise",
                                          "CD - Dichotomize Network",
                                          "CD - Accouting for Noise"))) %>% 
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Boxplot of average strength by group",
    x = "",
    y = "Average strength"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
