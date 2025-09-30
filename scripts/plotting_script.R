
# Load libraries ---------------------------------------------------------------

library(tidyverse)

# Load data --------------------------------------------------------------------

# setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out")
# setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out_false_pos_degree")
setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out_false_pos_degree_diff_random_clust")
load("main_res.RData")

# Create nice labels -----------------------------------------------------------

nice_labels <- tibble(
  type = c(
    "run_ABM_no_intervention",
    "run_ABM_isolate_individuals",
    "run_ABM_accounting_for_noise",
    "run_ABM_random_lshm",
    "run_ABM_random_lshm_outside",
    "run_ABM_random_lshm_size_and_number",
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
    "Clustering - Random (outside) LSHM",
    "Clustering - Random (number and size) LSHM",
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
                          "Clustering - Random (outside) LSHM",
                          "Clustering - Random (number and size) LSHM",
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
             "efficiency_per_quarantined_person_day",
             "mean_efficiency_per_quarantined_person_day",
             "average_times_quarantined",
             "attack_rate_vs_mean_quarantined_days",
             "total_infections_vs_person_days_in_quarantine",
             "infections_averted_vs_mean_quarantined_days"),
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
            "Total infections prevented per quarantined person-day (compared to baseline)",
            "Average total infections prevented per quarantined person-day (compared to baseline)",
            "Average number of times quarantied (among those ever quarantined)",
            "Attack rate by average time in quarantine",
            "Total infections by person time in quarantine",
            "Total infections prevented by average time in quarantine ")
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
                       "efficiency_per_quarantined_person_day",
                       "average_times_quarantined")){
    
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
      
  } else if (metric_val == "mean_efficiency_per_quarantined_person_day"){
    
    main_res %>% 
      filter(p_asym %in% p_asym_val) %>% 
      filter(p_asym<1) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity %in% DCT_specificity_val) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise")) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = factor(p_asym), y = mean,
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
        x = "Asymptomatic Probability",
        y = filter(nice_metric, metric_val == metric) %>% .$label,
        title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
        subtitle = paste0("Quarantined days = ", quarantine_days_val)
      )
    
  } else if (metric_val == "attack_rate_vs_mean_quarantined_days"){
    
    if(length(quarantine_days_val) > 1){
      stop("Please supply single value for quarantine_days_val")
    }
    
    main_res %>% 
      filter(round(p_asym, 1) %in% p_asym_val) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity == 1 ) %>% 
      mutate(p_asym = factor(p_asym)) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise")) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = mean_quarantined_days, y = attack_rates, color = label)) +
      geom_point(alpha = 0.6) +
      # geom_smooth(se = FALSE) +
      ggh4x::facet_grid2(DCT_sensitivity ~p_asym ,
                         labeller = labeller(
                           DCT_sensitivity = function(x) paste0("Sensitivity == ", x),
                           p_asym = function(x) paste0("Prob~asym. == ", x),
                           .default = label_parsed  # Needed to parse expressions
                         ),
                         scales = "free_y") +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
      ) +
      labs(x = "Mean Quarantine Days per Actor",
           y = "Attack Rate",
           title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
           subtitle = paste0("Quarantined days = ", quarantine_days_val))
    
  } else if (metric_val == "total_infections_vs_person_days_in_quarantine"){
    
    if(length(quarantine_days_val) > 1){
      stop("Please supply single value for quarantine_days_val")
    }
    
    main_res %>% 
      filter(round(p_asym, 1) %in% p_asym_val) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity == 1 ) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise",
                          "run_ABM_no_intervention")) %>% 
      mutate(p_asym = factor(p_asym)) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = quarantined_total, y = incidence_total, color = label)) +
      geom_point(alpha = 0.6) +
      # geom_smooth(se = FALSE) +
      ggh4x::facet_grid2(DCT_sensitivity ~p_asym ,
                         labeller = labeller(
                           DCT_sensitivity = function(x) paste0("Sensitivity == ", x),
                           p_asym = function(x) paste0("Prob~asym. == ", x),
                           .default = label_parsed  # Needed to parse expressions
                         ),
                         scales = "free_y") +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
      ) +
      labs(x = "Total Quarantined Person-days",
           y = "Total Number of Infections",
           title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
           subtitle = paste0("Quarantined days = ", quarantine_days_val))
    
  } else if (metric_val == "infections_averted_vs_mean_quarantined_days"){
    
    if(length(quarantine_days_val) > 1){
      stop("Please supply single value for quarantine_days_val")
    }
    
    main_res %>% 
      filter(round(p_asym, 1) %in% p_asym_val) %>% 
      filter(DCT_sensitivity %in% DCT_sensitivity_val) %>% 
      filter(DCT_specificity == 1 ) %>% 
      filter(!type %in% c("run_ABM_dichotomize_g_1",
                          "run_ABM_ignore_noise",
                          "run_ABM_no_intervention")) %>% 
      mutate(p_asym = factor(p_asym)) %>% 
      inner_join(nice_labels, by = "type") %>% 
      filter(quarantine_days == quarantine_days_val) %>% 
      filter(metric == metric_val) %>% 
      ggplot(aes(x = mean_quarantined_days, y = incidence_diff, color = label)) +
      geom_point(alpha = 0.6) +
      # geom_smooth(se = FALSE) +
      ggh4x::facet_grid2(DCT_sensitivity ~p_asym ,
                         labeller = labeller(
                           DCT_sensitivity = function(x) paste0("Sensitivity == ", x),
                           p_asym = function(x) paste0("Prob~asym. == ", x),
                           .default = label_parsed  # Needed to parse expressions
                         ),
                         scales = "free_y") +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
      ) +
      labs(x = "Mean Quarantine Days per Actor",
           y = "Total Infections Prevented (compared to no-intervention)",
           title = paste0(filter(nice_metric, metric_val == metric) %>% .$label),
           subtitle = paste0("Quarantined days = ", quarantine_days_val))
    
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
# "efficiency_per_quarantined_person_day"
# "average_times_quarantined"
# "attack_rate_vs_mean_quarantined_days"
# "infections_averted_vs_mean_quarantined_days"

plot_fun("attack_rate_vs_mean_quarantined_days",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 10)

plot_fun("infections_averted_vs_mean_quarantined_days",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("total_inf",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 10)

plot_fun("quarantined_total_person_days",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 10)

plot_fun("mean_efficiency_per_quarantined_person_day",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)


plot_fun("tot_inf_by_person_days_outside_qua",
         # p_asym_val = 0,
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("attack_rate",
         p_asym_val = c(0, .2, .4, .6, .8),
         DCT_specificity = 1,
         quarantine_days_val = 5)

plot_fun("quarantined_time",
         p_asym_val = 0,
         DCT_specificity = 1,
         quarantine_days_val = 15)
