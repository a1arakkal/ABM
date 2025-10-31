
# Load libraries ---------------------------------------------------------------

library(janitor)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(gridExtra)

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
    "DCT",
    "Clustering - LSHM",
    "Clustering - Random LSHM",
    "Clustering - Random (outside) LSHM",
    "Clustering - Random (number and size) LSHM",
    "Clustering - FB",
    "Clustering - Random FB",
    "Clustering - Ignore Noise",
    "Clustering - Thresholding"
  )
)
 

nice_labels <- nice_labels %>% 
  filter(type %in% c("run_ABM_no_intervention",
         "run_ABM_isolate_individuals",
         "run_ABM_accounting_for_noise",
         # "run_ABM_random_lshm",
         # "run_ABM_random_lshm_outside",
         # "run_ABM_random_lshm_size_and_number",
         "run_ABM_fb_clusters"
         # "run_ABM_random_fb",
         #"run_ABM_dichotomize_g_1",
         #"run_ABM_ignore_noise"
         ))

cbPalette <- c(   "#FFFF66", # bright yellow
                  "#00C78C", # bright bluish green
                  "#3399FF", # bright blue
                  "#FF6A33", # bright vermillion
                  "#E68CC7",  # bright reddish purple
                  "#FFB347" # brighter orange
                  
)

# Load data --------------------------------------------------------------------

# Figure 1 ---------------------------------------------------------------------

# Contact pattern over last 2 weeks repeated 4 time

# Pull in interaction data 
setwd("/Volumes/argon_home/dissertation/real_data_application/paper2/")

# https://figshare.com/articles/dataset/The_Copenhagen_Networks_Study_interaction_data/7267433/1
cns = 
  readr::read_csv("bt_symmetric.csv")

cns %<>%
  clean_names()

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

cns_fb = 
  readr::read_csv("fb_friends.csv")

cns_fb %<>%
  clean_names()

cns %<>%
  filter(user_b >= 0) # See readme.  Negative values are not meaningful.

# Unique user values
user_values <- tibble(user = sort(unique(c(cns$user_a,cns$user_b)))) %>% 
  mutate(name = 1:n())

total_actors <- nrow(user_values)

# Replace user values with new user values
cns %<>%
  left_join(user_values,
            by = join_by(user_a == user)) %>% 
  select(-user_a) %>% 
  rename(user_a = name) %>% 
  left_join(user_values,
            by = join_by(user_b == user)) %>% 
  select(-user_b) %>% 
  rename(user_b = name)

# collapse 5 min timestamps for intereactions into daily level
interactions_temp <- cns %>% 
  distinct(number_timestamp) %>%
  mutate(new_timestamp = (row_number() - 1) %/% 288 + 1) # each timestamp is 5 mins so 288 timestamps day, subtract 1 as we want groups to be closed on right interval
# table(cut(1:8064, breaks = c(0, seq(288, 288*4*7, 288)) )

# tibble with total actor pair interactions by day
interactions <- interactions_temp %>% 
  inner_join(cns, ., by = "number_timestamp") %>% 
  rename(timestamp = new_timestamp) %>% 
  count(timestamp, user_a, user_b)

data <- interactions %>% 
  group_by(timestamp) %>% 
  summarise(total = sum(n), .groups = "drop") %>% 
  mutate(time = row_number()) 

fig1 <- data %>% 
  ggplot(aes(x = time, y = total)) +
  geom_line() +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank(),
        strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
        strip.text.x = element_text(size = 14),
        text = element_text(size = 15)) +
  labs(
    x = "Time (weeks)",
    y = "Total number of interactions"
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_x_continuous(breaks = seq(3.5, max(data$time)-3.5, 7),
                     labels = paste0("Week ", 1:(max(data$time)/7))) +
  geom_vline(xintercept = seq(7, max(data$time), 7), color = "black", lty = "dashed") 


jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig1.jpeg",
     width = 4.8, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig1)
dev.off()


# Figure 2 ---------------------------------------------------------------------

# Contact pattern over last 2 weeks repeated 4 time

interactions_last_2_weeks <- interactions %>% 
  filter(timestamp > 28/2) %>% 
  group_by(timestamp) %>% 
  summarise(total = sum(n), .groups = "drop")

data <- rbind(interactions_last_2_weeks,
              interactions_last_2_weeks,
              interactions_last_2_weeks,
              interactions_last_2_weeks) %>% 
  mutate(time = row_number()) 

fig2 <- data %>% 
  ggplot(aes(x = time, y = total)) +
  geom_line() +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank(),
        strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
        strip.text.x = element_text(size = 14),
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Time (weeks)",
    y = "Total number of interactions"
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_x_continuous(breaks = seq(3.5, max(data$time)-3.5, 7),
                     labels = paste0("Week ", 1:(max(data$time)/7))) +
  geom_vline(xintercept = seq(7, max(data$time), 7), color = "black", lty = "dashed") 

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig2.jpeg",
     width = 4.8, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig2)
dev.off()

# Figure 3 ---------------------------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

load(file = "cns_week_1_2_res.RData")
cns_fits_lshm <- cns_fits
rm(cns_fits_ignore_noise, cns_fits_di, cns_fits)

load(file = "fb_net_LSM_res.RData")
cns_fits_fb <- cns_fits

# Plot network LSHM 

cns_el <- interactions %>% 
  filter(timestamp <= 28/2) %>% # subset to interaction in first 2 weeks
  group_by(user_a, user_b) %>% 
  summarise(n = sum(n), .groups = "drop")

cns_net =
  tbl_graph(nodes = 
              user_values %>% select(name),
            edges = cns_el,
            directed = FALSE)

cbPalette <- c(   "#FFFF66", # bright yellow
                  "#00C78C", # bright bluish green
                  "#3399FF", # bright blue
                  "#FF6A33", # bright vermillion
                  "#E68CC7",  # bright reddish purple
                  "#FFB347" # brighter orange
                  
)

# LSHM
p1 = 
  cns_net %>%
  activate(nodes) %>% 
  filter(!name %in% setdiff(user_values$name, colnames(cns_fits_lshm$A))) %>% 
  mutate(cluster = 
           factor(cns_fits_lshm$optimal_res$cluster_labels, 
                  levels = unique(cns_fits_lshm$optimal_res$cluster_labels), 
                  labels = paste0("Cluster ", 1:length(unique(cns_fits_lshm$optimal_res$cluster_labels))))) %>% 
  ggraph::ggraph(layout = cns_fits_lshm$optimal_res$U) + 
  ggraph::geom_edge_bundle_path(tension = 0.6,
                                alpha = 0.1) +
  # ggraph::geom_edge_link(alpha = 0.1) +
  ggraph::geom_node_point(aes(fill = cluster),
                          shape = 21,
                          size = 2,
                          color = "black",
                          stroke = 0.5) +
  scale_fill_manual(values = cbPalette) +   
  ggtitle("(A)")+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 12),
        legend.position = "none")


# Plot network FB 

cns_fb %<>% 
  rename(user_a = number_user_a) %>% 
  filter(user_a != user_b) # 11 cases where they are friends with themselves???

# Replace user values with new user values
cns_fb %<>%
  filter( (user_a %in% user_values$user) &  (user_b %in% user_values$user)) %>% 
  left_join(user_values,
            by = join_by(user_a == user)) %>% 
  select(-user_a) %>% 
  rename(user_a = name) %>% 
  left_join(user_values,
            by = join_by(user_b == user)) %>% 
  select(-user_b) %>% 
  rename(user_b = name) %>% 
  mutate(value = 1L)

cns_net =
  tbl_graph(nodes = 
              user_values %>% select(name),
            edges = cns_fb,
            directed = FALSE)

cbPalette <- c(   "#FFFF66", # bright yellow
                  "#00C78C", # bright bluish green
                  "#3399FF", # bright blue
                  "#FF6A33", # bright vermillion
                  "#E68CC7",  # bright reddish purple
                  "#FFB347" # brighter orange
                  
)

# fb
p2 = 
  cns_net %>%
  activate(nodes) %>% 
  filter(!name %in% setdiff(user_values$name, colnames(cns_fits_fb$A))) %>% 
  mutate(cluster = 
           factor(cns_fits_fb$optimal_res$cluster_labels, 
                  levels = unique(cns_fits_fb$optimal_res$cluster_labels), 
                  labels = paste0("Cluster ", 1:length(unique(cns_fits_fb$optimal_res$cluster_labels))))) %>% 
  ggraph::ggraph(layout = cns_fits_fb$optimal_res$U) + 
  ggraph::geom_edge_bundle_path(tension = 0.6,
                                alpha = 0.1) +
  # ggraph::geom_edge_link(alpha = 0.1) +
  ggraph::geom_node_point(aes(fill = cluster),
                          shape = 21,
                          size = 2,
                          color = "black",
                          stroke = 0.5) +
  scale_fill_manual(values = cbPalette) +   
  ggtitle("(B)")+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 12),
        legend.position = "none")


# Arrange
fig3 <- grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  ncol = 1,
  heights = 2
)

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig3.jpeg",
     width = 4.8, height = 3, units = "in", res = 1000)
grid::grid.draw(fig3)
dev.off()



# Figure 4 ---------------------------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_out_false_pos_degree_diff_random_clust")

colors_group <- c("No Intervention" = "#FF6A33", 
                  "DCT (Sen. = 0.6)" = "#00C78C",
                  "DCT (Sen. = 0.8)" = "#E68CC7",
                  "DCT (Sen. = 1.0)" = "#FFFF66",
                  "Clustering - LSHM" = "#3399FF",
                  "Clustering - FB" = "#FFB347")
load("main_res.RData")

# Attack rate (%) by intervention

fig4 <- main_res %>% 
  filter(quarantine_days == 5) %>% 
  filter(metric == "attack_rate") %>% 
  filter(round(p_asym, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>% 
  filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
  filter(DCT_specificity == 1) %>%
  inner_join(nice_labels, by = "type") %>% 
  mutate( label = ifelse(type == "run_ABM_isolate_individuals",
                         paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
                         label) ) %>% 
  select(2:7, p_asym, label) %>% 
  distinct() %>% 
  mutate(label = factor(label, 
                      levels = c(
                        "No Intervention",
                        "DCT (Sen. = 0.6)",
                        "DCT (Sen. = 0.8)",
                        "DCT (Sen. = 1.0)",
                        "Clustering - LSHM",
                        "Clustering - Random LSHM",
                        "Clustering - Random (outside) LSHM",
                        "Clustering - Random (number and size) LSHM",
                        "Clustering - FB",
                        "Clustering - Random FB",
                        "Clustering - Ignore Noise",
                        "Clustering - Thresholding"
                      ))) %>% 
  # count(label) %>% 
  ggplot(aes(x = factor(p_asym), fill = label)) +
  geom_boxplot(
    aes(
      ymin = ymin*100,
      lower = Q1*100,
      middle = median*100,
      upper = Q3*100,
      ymax = ymax*100
    ),
    stat = "identity",
    alpha = 1
  ) +
  scale_fill_manual(values = colors_group) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
    strip.text.x = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Asymptomatic probability",
    y = "Attack rate (%)"
  )+
  geom_vline(xintercept = seq(1.5, length(c(0, 0.2, 0.4, 0.6, 0.8)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig4.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig4)
dev.off()


# Figure 5 ---------------------------------------------------------------------

# Total Number of Contacts Across 56 Days

fig5 <- main_res %>% 
  filter(quarantine_days == 5) %>% 
  filter(metric == "average_cumulative_interactions_per_actor") %>% 
  filter(round(p_asym, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>% 
  filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
  filter(DCT_specificity == 1) %>%
  inner_join(nice_labels, by = "type") %>% 
  mutate( label = ifelse(type == "run_ABM_isolate_individuals",
                         paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
                         label) ) %>% 
  select(2:7, p_asym, label) %>% 
  distinct() %>% 
  mutate(label = factor(label, 
                        levels = c(
                          "No Intervention",
                          "DCT (Sen. = 0.6)",
                          "DCT (Sen. = 0.8)",
                          "DCT (Sen. = 1.0)",
                          "Clustering - LSHM",
                          "Clustering - Random LSHM",
                          "Clustering - Random (outside) LSHM",
                          "Clustering - Random (number and size) LSHM",
                          "Clustering - FB",
                          "Clustering - Random FB",
                          "Clustering - Ignore Noise",
                          "Clustering - Thresholding"
                        ))) %>% 
  # count(label) %>% 
  ggplot(aes(x = factor(p_asym), fill = label)) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = ymax
    ),
    stat = "identity",
    alpha = 1
  ) +
  scale_fill_manual(values = colors_group) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
    strip.text.x = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Asymptomatic probability",
    y = "Average number of contacts per actor"
  )+
  geom_vline(xintercept = seq(1.5, length(c(0, 0.2, 0.4, 0.6, 0.8)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig5.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig5)
dev.off()


# Figure 6 ---------------------------------------------------------------------

fig6 <- main_res %>% 
  filter(quarantine_days == 5) %>% 
  filter(metric == "tot_inf_by_by_contacts") %>% 
  filter(round(p_asym, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>% 
  filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
  filter(DCT_specificity == 1) %>%
  inner_join(nice_labels, by = "type") %>% 
  mutate( label = ifelse(type == "run_ABM_isolate_individuals",
                         paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
                         label) ) %>% 
  select(2:7, p_asym, label) %>% 
  distinct() %>% 
  mutate(label = factor(label, 
                        levels = c(
                          "No Intervention",
                          "DCT (Sen. = 0.6)",
                          "DCT (Sen. = 0.8)",
                          "DCT (Sen. = 1.0)",
                          "Clustering - LSHM",
                          "Clustering - Random LSHM",
                          "Clustering - Random (outside) LSHM",
                          "Clustering - Random (number and size) LSHM",
                          "Clustering - FB",
                          "Clustering - Random FB",
                          "Clustering - Ignore Noise",
                          "Clustering - Thresholding"
                        ))) %>% 
  # count(label) %>% 
  filter(label != "No Intervention") %>% 
  ggplot(aes(x = factor(p_asym), fill = label)) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = ymax
    ),
    stat = "identity",
    alpha = 1
  ) +
  scale_fill_manual(values = colors_group) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
    strip.text.x = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Asymptomatic probability",
    y = "Ratio of total infections to\ntotal number of contacts"
  )+
  geom_vline(xintercept = seq(1.5, length(c(0, 0.2, 0.4, 0.6, 0.8)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig6.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig6)
dev.off()

# Figure 7 ---------------------------------------------------------------------
setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/res_modularity")
rm(main_res)
load("main_res.RData")


colors_group <- c("No Intervention" = "#FF6A33", 
                  "DCT (Sen. = 0.6)" = "#00C78C",
                  "DCT (Sen. = 0.8)" = "#E68CC7",
                  "DCT (Sen. = 1.0)" = "#FFFF66",
                  "Clustering" = "#3399FF",
                  "Clustering - Random" = "#FFB347")

# main_res %>% 
#   filter(cluster_label == 'clusters_accounting_for_noise') %>% 
#   filter(quarantine_days == 5) %>% 
#   filter(metric == "tot_inf_by_by_contacts") %>% 
#   filter(round(p_asym, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>% 
#   filter(round(modularity, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8, 1)) %>% 
#   filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
#   filter(DCT_specificity == 1) %>%
#   mutate(label = case_when(
#     type == "run_ABM_no_intervention" ~ "No Intervention",
#     type == "run_ABM_isolate_individuals" ~ "DCT",
#     type == "run_ABM_clustering" ~ "Clustering",
#     type == "run_ABM_random_clustering" ~ "Clustering - Random",
#   )) %>% 
#   mutate( label = ifelse(type == "run_ABM_isolate_individuals",
#                          paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
#                          label) ) %>% 
#   filter(label != "No Intervention") %>% 
#   select(2:7, p_asym, modularity, label) %>% 
#   distinct() %>% 
#  # group_by(p_asym, modularity, label) %>%
#  # filter(median == min(median)) %>%
#  # ungroup() %>% 
#   mutate(label = factor(label, 
#                         levels = c(
#                           "No Intervention",
#                           "DCT (Sen. = 0.6)",
#                           "DCT (Sen. = 0.8)",
#                           "DCT (Sen. = 1.0)",
#                           "Clustering",
#                           "Clustering - Random"
#                         ))) %>% 
#   # count(label) %>%
#   ggplot(aes(x = factor(modularity), fill = label)) +
#   geom_boxplot(
#     aes(
#       ymin = ymin,
#       lower = Q1,
#       middle = median,
#       upper = Q3,
#       ymax = ymax
#     ),
#     stat = "identity",
#     alpha = 1
#   ) +
#   scale_fill_manual(values = colors_group) +
#   ggh4x::facet_grid2(~p_asym ,
#                      labeller = labeller(
#                        # modularity = function(x) paste0("Modularity == ", x),
#                        p_asym = function(x) paste0("Prob~asym. == ", x),
#                        .default = label_parsed  # Needed to parse expressions
#                      ),
#                      scales = "free_y") +
#   theme_bw() +
#   theme(
#     # axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid.major.x = element_blank(),
#     strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
#     strip.text.x = element_text(size = 14),
#     text = element_text(size = 14)
#   ) +
#   labs(
#     x = "Modularity",
#     y = "Ratio of total infections versus total number of contacts"
#   )+
#   geom_vline(xintercept = seq(1.5, length(c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) - 0.5, by = 1),
#              linetype = "dashed", color = "grey50")+
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))



fig7 <- main_res %>% 
  filter(cluster_label == 'clusters_accounting_for_noise') %>% 
  filter(quarantine_days == 5) %>% 
  filter(metric == "tot_inf_by_by_contacts") %>% 
  filter(round(p_asym, 1) %in% c(0.4)) %>% 
  filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
  filter(DCT_specificity == 1) %>%
  mutate(label = case_when(
    type == "run_ABM_no_intervention" ~ "No Intervention",
    type == "run_ABM_isolate_individuals" ~ "DCT",
    type == "run_ABM_clustering" ~ "Clustering",
    type == "run_ABM_random_clustering" ~ "Clustering - Random",
  )) %>% 
  mutate( label = ifelse(type == "run_ABM_isolate_individuals",
                         paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
                         label) ) %>% 
  filter(label != "No Intervention") %>% 
  select(2:7, p_asym, modularity, label) %>% 
  # mutate(modularity = round(modularity*100)) %>% 
  distinct() %>% 
  mutate(label = factor(label, 
                        levels = c(
                          "No Intervention",
                          "DCT (Sen. = 0.6)",
                          "DCT (Sen. = 0.8)",
                          "DCT (Sen. = 1.0)",
                          "Clustering",
                          "Clustering - Random"
                        ))) %>% 
  # count(label) %>%
  ggplot(aes(x = factor(modularity), fill = label)) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = ymax
    ),
    stat = "identity",
    alpha = 1
  ) +
  scale_fill_manual(values = colors_group) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
    strip.text.x = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Proportion of edges within cluster",
    y = "Ratio of total infections to\ntotal number of contacts"
  )+
  geom_vline(xintercept = seq(1.5, length(seq(0, 1, 0.1)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_y_continuous(labels = scales::scientific)


jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig7.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig7)
dev.off()

# Figure 8 ---------------------------------------------------------------------
# 
# main_res %>%
#   filter(cluster_label == 'clusters_fb_network') %>%
#   filter(quarantine_days == 5) %>%
#   filter(metric == "tot_inf_by_by_contacts") %>%
#   filter(round(p_asym, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8)) %>%
#   filter(round(modularity, 1) %in% c(0, 0.2, 0.4, 0.6, 0.8, 1)) %>%
#   filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
#   filter(DCT_specificity == 1) %>%
#   mutate(label = case_when(
#     type == "run_ABM_no_intervention" ~ "No Intervention",
#     type == "run_ABM_isolate_individuals" ~ "DCT",
#     type == "run_ABM_clustering" ~ "Clustering",
#     type == "run_ABM_random_clustering" ~ "Clustering - Random",
#   )) %>%
#   mutate( label = ifelse(type == "run_ABM_isolate_individuals",
#                          paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
#                          label) ) %>%
#   filter(label != "No Intervention") %>%
#   select(2:7, p_asym, modularity, label) %>%
#   distinct() %>%
#  # group_by(p_asym, modularity, label) %>%
#  # filter(median == min(median)) %>%
#  # ungroup() %>% 
#   mutate(label = factor(label,
#                         levels = c(
#                           "No Intervention",
#                           "DCT (Sen. = 0.6)",
#                           "DCT (Sen. = 0.8)",
#                           "DCT (Sen. = 1.0)",
#                           "Clustering",
#                           "Clustering - Random"
#                         ))) %>%
#   # count(label) %>%
#   ggplot(aes(x = factor(modularity), fill = label)) +
#   geom_boxplot(
#     aes(
#       ymin = ymin,
#       lower = Q1,
#       middle = median,
#       upper = Q3,
#       ymax = ymax
#     ),
#     stat = "identity",
#     alpha = 1
#   ) +
#   scale_fill_manual(values = colors_group) +
#   ggh4x::facet_grid2(~p_asym ,
#                      labeller = labeller(
#                        # modularity = function(x) paste0("Modularity == ", x),
#                        p_asym = function(x) paste0("Prob~asym. == ", x),
#                        .default = label_parsed  # Needed to parse expressions
#                      ),
#                      scales = "free_y") +
#   theme_bw() +
#   theme(
#     # axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid.major.x = element_blank(),
#     strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
#     strip.text.x = element_text(size = 14),
#     text = element_text(size = 14)
#   ) +
#   labs(
#     x = "Modularity",
#     y = "Ratio of total infections versus total number of contacts"
#   )+
#   geom_vline(xintercept = seq(1.5, length(c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) - 0.5, by = 1),
#              linetype = "dashed", color = "grey50")+
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE))



fig8 <- main_res %>% 
  filter(cluster_label == 'clusters_fb_network') %>% 
  filter(quarantine_days == 5) %>% 
  filter(metric == "tot_inf_by_by_contacts") %>% 
  filter(round(p_asym, 1) %in% c(0.4)) %>% 
  filter(DCT_sensitivity %in% c(0.6, 0.8, 1)) %>%
  filter(DCT_specificity == 1) %>%
  mutate(label = case_when(
    type == "run_ABM_no_intervention" ~ "No Intervention",
    type == "run_ABM_isolate_individuals" ~ "DCT",
    type == "run_ABM_clustering" ~ "Clustering",
    type == "run_ABM_random_clustering" ~ "Clustering - Random",
  )) %>% 
  mutate( label = ifelse(type == "run_ABM_isolate_individuals",
                         paste0(label, " (Sen. = ", ifelse(DCT_sensitivity == 1, "1.0", DCT_sensitivity), ")"),
                         label) ) %>% 
  filter(label != "No Intervention") %>% 
  select(2:7, p_asym, modularity, label) %>% 
  distinct() %>% 
  # mutate(modularity = round(modularity*100)) %>%  
  mutate(label = factor(label, 
                        levels = c(
                          "No Intervention",
                          "DCT (Sen. = 0.6)",
                          "DCT (Sen. = 0.8)",
                          "DCT (Sen. = 1.0)",
                          "Clustering",
                          "Clustering - Random"
                        ))) %>%
  # count(label) %>%
  ggplot(aes(x = factor(modularity), fill = label)) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = ymax
    ),
    stat = "identity",
    alpha = 1
  ) +
  scale_fill_manual(values = colors_group) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    strip.text.y = element_text(face = "bold", size = 14),  # bold + bigger for row labels (family)
    strip.text.x = element_text(size = 14),
    text = element_text(size = 14)
  ) +
  labs(
    x = "Proportion of edges within cluster",
    y = "Ratio of total infections to\ntotal number of contacts"
  )+
  geom_vline(xintercept = seq(1.5, length(seq(0, 1, 0.1)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_y_continuous(labels = scales::scientific)


jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig8.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig8)
dev.off()

# Example of change modularity -------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")

int_and_neighbors_by_t_true <- int_and_neighbors_by_t["20"]

# Load in LSHM results

load(file = "cns_week_1_2_res.RData")

clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise")))

load(file = "fb_net_LSM_res.RData")
clusters_fb_network <- JANE:::summary.JANE(cns_fits)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise",
                          "clusters_fb_network")))

# Source helper functions 

modularity <- function(network,
                       total_actors,
                       assumed_true_communities){
  
  out <-  sapply(names(network), 
                 function(x){
                   
                   A <- matrix(0, nrow = total_actors, ncol= total_actors)
                   colnames(A) <- 1:total_actors
                   rownames(A) <- 1:total_actors
                   dat <- network[[x]]$neighbors
                   
                   for(i in names(dat)){
                     A[i, names(dat[[i]])] <- 1
                   }
                   
                   sum_within <- 0
                   for (i in unique(assumed_true_communities)){
                     sum_within <- sum_within + sum(A[names(which(assumed_true_communities == i)),
                                                      names(which(assumed_true_communities == i))])
                   }
                   
                   sum_between <- sum(A) - sum_within
                   
                   return(c(sum_within, sum_between, sum_within/(sum(A))))
                 })
  
}

modularity_clusters_accounting_for_noise <- modularity(int_and_neighbors_by_t_true,
                                                       total_actors,
                                                       clusters_accounting_for_noise)

modularity_clusters_fb_network <- modularity(int_and_neighbors_by_t_true,
                                             total_actors,
                                             clusters_fb_network)

# Function to adjust modularity ------------------------------------------------

adjust_modularity <- function(network,
                              total_actors,
                              assumed_true_communities,
                              prob_within,
                              tol = 1e-2){
  
  new_network <- vector("list", length = length(network))
  names(new_network) <- names(network)
  modularity <- numeric(length = length(network))
  names(modularity) <- names(network)
  
  for (t in names(network)){
    
    # extract within and between edges -----------------------------------------
    A <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(A) <- 1:total_actors
    rownames(A) <- 1:total_actors
    dat <- network[[t]]$neighbors
    
    for(i in names(dat)){
      A[i, names(dat[[i]])] <- unname(dat[[i]])
    }
    
    # isSymmetric(A)
    indices_w_ints <- which(A>0, arr.ind = TRUE)
    indices_w_ints <- cbind(indices_w_ints, 0)
    colnames(indices_w_ints) <- c("i", "j", "within")
    indices_w_ints <- indices_w_ints[indices_w_ints[, "j"] > indices_w_ints[, "i"], ]
    indices_w_ints <- cbind(indices_w_ints, A[indices_w_ints[, 1:2]])
    colnames(indices_w_ints) <- c("i", "j", "within", "value")
    
    within_loc <- (as.character(indices_w_ints[,"i"]) %in% names(assumed_true_communities)) & 
      (as.character(indices_w_ints[,"j"]) %in% names(assumed_true_communities)) &
      (assumed_true_communities[as.character(indices_w_ints[,"i"])] ==  assumed_true_communities[as.character(indices_w_ints[,"j"])])
    within_loc[is.na(within_loc)] <- FALSE
    indices_w_ints[within_loc, "within"] <- 1
    
    obs_p_within <- prop.table(table(indices_w_ints[, "within"]))["1"]
    within_edges  <- indices_w_ints[indices_w_ints[, "within"] == 1, ]
    between_edges <-  indices_w_ints[indices_w_ints[, "within"] == 0, ]
    
    if(!is.null(prob_within)){
      
      if(prob_within == 1){
        
        out <- within_edges
        
      } else if (prob_within == 0){
        
        out <- between_edges
        
      } else if ( abs(prob_within - obs_p_within) < tol){
        
        out <- indices_w_ints
        
      } else {
        
        if (obs_p_within < prob_within) {
          
          # need to remove some between edges
          within <- within_edges
          N_between <- round((nrow(within) * (1 - prob_within)) / prob_within)
          between <- between_edges[order(runif(nrow(between_edges))), ]
          out <- rbind(within, between[1:N_between, ])
          
        } else {
          
          # need to remove some within edges
          between <- between_edges
          N_within <- round((prob_within * nrow(between)) / (1 - prob_within))
          within <- within_edges[order(runif(nrow(within_edges))), ]
          out <- rbind(within[1:N_within, ], between)
          
        }
        
      } 
      
    } else {
      
      out <- indices_w_ints
      
    }
    
    B <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(B) <- 1:total_actors
    rownames(B) <- 1:total_actors
    
    B[out[, 1:2]] <- out[, "value"]
    B[out[, 2:1]] <- out[, "value"]
    # isSymmetric(B)
    n_total_t <- rowSums(B)
    check_isolates <- n_total_t == 0
    
    sum_within <- 0
    for (i in unique(assumed_true_communities)){
      sum_within <- sum_within + sum((B>0)[names(which(assumed_true_communities == i)),
                                           names(which(assumed_true_communities == i))])
    }
    modularity[t] <- sum_within/sum(B>0)
    # modularity[t] <- prop.table(table(out[, "within"]))[["1"]]
    
    temp2 <- lapply(rownames(B)[!check_isolates], function(y){
      z <- B[y, ]
      z <- z[z!=0]
      return(z)
    })
    
    names(temp2) <- rownames(B)[!check_isolates]
    
    new_network[[t]] <- list(neighbors = temp2,
                             n_total_t = n_total_t[!check_isolates])
    
  }
  
  return(list(network = new_network,
              modularity = modularity))
  
}

add_cluster_boundaries <- function(p, cluster_vec) {
  boundaries <- cumsum(table(cluster_vec))
  n <- length(cluster_vec)
  
  for (b in boundaries[-length(boundaries)]) {
    p <- p +
      geom_vline(xintercept = b, linetype="solid", 
                 color = "black", size=0.5)+
      geom_hline(yintercept = n - b, linetype="solid", 
                 color = "black", size=0.5)
  }
  return(p)
}

plot_fun <- function(network,
                     total_actors,
                     assumed_true_communities,
                     prob_within,
                     tol = 1e-2,
                     t){

  if (is.null(prob_within)){

    A <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(A) <- 1:total_actors
    rownames(A) <- 1:total_actors
    dat <- network[[t]]$neighbors

    for(i in names(dat)){
      A[i, names(dat[[i]])] <- unname(dat[[i]])
    }
    
    } else {

    temp <- adjust_modularity(network = network,
                              total_actors = total_actors,
                              assumed_true_communities = assumed_true_communities,
                              prob_within = prob_within,
                              tol = tol-2)

    A <- matrix(0, nrow = total_actors, ncol= total_actors)
    colnames(A) <- 1:total_actors
    rownames(A) <- 1:total_actors
    dat <- temp$network[[t]]$neighbors

    for(i in names(dat)){
      A[i, names(dat[[i]])] <- unname(dat[[i]])
    }

  }

  #cns_fits_di
  plot_mat <- as.matrix(A)
  order_clust <- sort(assumed_true_communities)
  
  adj_plot <- as.data.frame(as.table(plot_mat)) %>%
    rename(row = Var1, col = Var2) %>% 
    mutate(value = cut(Freq, breaks = c(-1, 0, Inf),
                       labels = c("0", ">0"))) %>% 
    mutate(
      row = factor(row, levels = rev(names(order_clust))),
      col = factor(col, levels = names(order_clust))
    ) %>% 
    ggplot(aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", ">0" = "#FF6A33")) + 
    geom_tile() +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5,
                                size = 12)
    ) +
    ylab("")+
    xlab("")
  
  adj_plot <- add_cluster_boundaries(adj_plot, order_clust) 
  
  return(adj_plot)
}

p1 <- plot_fun(network = int_and_neighbors_by_t_true,
         total_actors = total_actors,
         assumed_true_communities = clusters_accounting_for_noise,
         prob_within = NULL,
         t = "20") +
  ggtitle("(A)")


p2 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_accounting_for_noise,
               prob_within = 0,
               t = "20") +
  ggtitle("(B)")


p3 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_accounting_for_noise,
               prob_within = .5,
               t = "20") +
  ggtitle("(C)")


p4 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_accounting_for_noise,
               prob_within = 1,
               t = "20") +
  ggtitle("(D)")


# Arrange 
fig9 <- grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  arrangeGrob(p3, p4, ncol = 2),
  ncol = 1,
  heights = c(2, 2)
)

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig9.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig9)
dev.off()



p1 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_fb_network,
               prob_within = NULL,
               t = "20") +
  ggtitle("(A)")


p2 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_fb_network,
               prob_within = 0,
               t = "20") +
  ggtitle("(B)")


p3 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_fb_network,
               prob_within = .5,
               t = "20") +
  ggtitle("(C)")


p4 <- plot_fun(network = int_and_neighbors_by_t_true,
               total_actors = total_actors,
               assumed_true_communities = clusters_fb_network,
               prob_within = 1,
               t = "20") +
  ggtitle("(D)")


# Arrange 
fig10 <- grid.arrange(
  arrangeGrob(p1, p2, ncol = 2),
  arrangeGrob(p3, p4, ncol = 2),
  ncol = 1,
  heights = c(2, 2)
)

jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig10.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig10)
dev.off()


### Data for text --------------------------------------------------------------

# Load in interactions by time list --------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

# Load in LSHM results --------------------------------------------

load(file = "cns_week_1_2_res.RData")

clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

isolates <- 692-nrow(cns_fits$A)
sum(cns_fits$A>0)/2
sum(cns_fits$A>0)/(nrow(cns_fits$A)*(nrow(cns_fits$A)-1))
length(table((clusters_accounting_for_noise)))
sort(table(clusters_accounting_for_noise))
rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise")))

load(file = "fb_net_LSM_res.RData")
clusters_fb_network <- JANE:::summary.JANE(cns_fits)$cluster_labels

isolates <- 692-nrow(cns_fits$A)
sum(cns_fits$A>0)/2
sum(cns_fits$A>0)/(nrow(cns_fits$A)*(nrow(cns_fits$A)-1))
length(table((clusters_fb_network)))
sort(table(clusters_fb_network))
