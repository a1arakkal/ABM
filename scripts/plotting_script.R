
# Load libraries ---------------------------------------------------------------

library(tidyverse)

# Load data --------------------------------------------------------------------

# setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out")
setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out_false_pos_degree")
load("main_res.RData")

# Create nice labels -----------------------------------------------------------

nice_labels <- tibble(
  type = c(
    "run_ABM_no_intervention",
    "run_ABM_isolate_individuals",
    "run_ABM_accounting_for_noise",
    "run_ABM_random_lshm",
    "run_ABM_fb_clusters",
    "run_ABM_random_fb",
    "run_ABM_dichotomize_g_1",
    "run_ABM_ignore_noise"
  ),
  label = c(
    "No Intervention",
    "Digital Contact Tracing",
    "Clustering - LSHM",
    "Clustering - Random LSHM",
    "Clustering - FB",
    "Clustering - Random FB",
    "Clustering - Ignore Noise",
    "Clustering - Thresholding"
  )
) %>% 
  mutate(label = factor(label, 
                        levels = c(
                          "No Intervention",
                          "Digital Contact Tracing",
                          "Clustering - LSHM",
                          "Clustering - Random LSHM",
                          "Clustering - FB",
                          "Clustering - Random FB",
                          "Clustering - Ignore Noise",
                          "Clustering - Thresholding"
                        )))


nice_metric <- tibble(
  metric = c("R0",
             "attack_rate",
             "average_cumulative_interactions_per_actor",
             "average_time_in_quarantine",
             "incidence_max",
             "quarantined_max",
             "quarantined_rate",
             "incidence_time",
             "quarantined_time",
             "tot_inf_by_by_contacts",
             "tot_inf_by_person_days_outside_qua",
             "total_inf",
             "quarantined_total_person_days",
             "total_contacts",
             "efficiency_per_contact"),
  label = c("R0",
            "Attack Rate (%)",
            "Average cumulative number of contacts per actor",
            "Average cumulative time in quarantine per actor",
            "Max number of incident cases",
            "Max number of actors quarantined at any one time",
            "Percent of actors ever quarantined",
            "Median number of incident cases by days",
            "Median number of actors in quarantined by days",
            "Ratio of total number of infections by total number of contacts",
            "Ratio of total number of infections by total person-days of freedom",
            "Total number of infections",
            "Total number of quarantined person-days",
            "Total number of contacts",
            "Total infections prevented by total contacts prevented (compared to baseline)")
)
  
# Plot function  ---------------------------------------------------------------

plot_fun <- function(metric_val,
                     quarantine_days_val = 5,
                     p_asym_val = c(0, 0.4, 0.8, 1),
                     DCT_sensitivity_val = c(0, 0.6, 0.8, 1),
                     DCT_specificity_val =  c(0, 0.6, 0.8, 1)){
  
  if(metric_val %in% c("R0",
                       "attack_rate",
                       "average_cumulative_interactions_per_actor",
                       "average_time_in_quarantine",
                       "incidence_max",
                       "quarantined_max",
                       "quarantined_rate",
                       "tot_inf_by_by_contacts",
                       "tot_inf_by_person_days_outside_qua",
                       "total_inf",
                       "quarantined_total_person_days",
                       "total_contacts",
                       "efficiency_per_contact")){
    
    if(length(quarantine_days_val) > 1){
      stop("Please supply single value for quarantine_days_val")
    }
    
    main_res %>% 
      filter(round(p_asym, 1) %in% p_asym_val) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity %in% DCT_specificity_val) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise")) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = factor(p_asym), color = label, fill = label)) +
      geom_boxplot(
        aes(
          ymin = ymin,
          lower = Q1,
          middle = median,
          upper = Q3,
          ymax = ymax
        ),
        stat = "identity",
        alpha = 0.4
      ) +
      ggh4x::facet_grid2(DCT_sensitivity ~ DCT_specificity,
                 labeller = labeller(
                   DCT_sensitivity = function(x) paste0("Sensitivity == ", x),
                   DCT_specificity = function(x) paste0("Specificity == ", x),
                   .default = label_parsed  # Needed to parse expressions
                 ),
                 scales = "free_y",
                 independent = "y") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
      ) +
      labs(
        x = "Asymptomatic Probability",
        y = filter(nice_metric, metric_val == metric) %>% .$label,
        title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
        subtitle = paste0("Quarantined days = ", quarantine_days_val)
                       
      )
  } else if (metric_val %in% c("quarantined_time", "incidence_time")){
    
    if(length(p_asym_val) > 1){
      stop("Please supply single value for p_asym_val")
    }
    
    if(length(quarantine_days_val) > 1){
      stop("Please supply single value for quarantine_days_val")
    }
    
    main_res %>% 
      filter(p_asym %in% p_asym_val) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity %in% DCT_specificity_val) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise",
                          "run_ABM_no_intervention")) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = t, y = median,
                 color = label, group = label)) +
      geom_line()+
      ggh4x::facet_grid2(DCT_sensitivity ~ DCT_specificity,
                 labeller = labeller(
                   DCT_sensitivity = function(x) paste0("Sensitivity == ", x),
                   DCT_specificity = function(x) paste0("Specificity == ", x),
                   .default = label_parsed  # Needed to parse expressions
                 ),
                 scales = "free_y",
                 independent = "y") +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
      ) +
      labs(
        x = "Time (days)",
        y = filter(nice_metric, metric_val == metric) %>% .$label,
        title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
        subtitle = paste0("Quarantined days = ",
                       quarantine_days_val, " & ",
                       "Probability asymptomatic = ",
                       p_asym_val)
      )
      
  } else {
    stop("Metric not avaliable")
  }
}

# Available metrics
# "R0"
# "attack_rate"
# "average_cumulative_interactions_per_actor"
# "average_time_in_quarantine"
# "incidence_max"
# "quarantined_max"
# "quarantined_rate"
# "incidence_time"
# "quarantined_time"
# "tot_inf_by_by_contacts"
# "tot_inf_by_person_days_outside_qua"
# "total_inf"
# "quarantined_total_person_days"
# "total_contacts"
# "efficiency_per_contact"

plot_fun("tot_inf_by_by_contacts",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("total_contacts",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("total_inf",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("quarantined_total_person_days",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("tot_inf_by_by_contacts",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("quarantined_total_person_days",
         # p_asym_val = 0,
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("quarantined_time",
         p_asym_val = 0,
         quarantine_days_val = 5)
