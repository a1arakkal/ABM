
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
    y = "Average number of contact per actor"
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
    y = "Ratio of total infections versus\ntotal number of contacts"
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
  filter(round(p_asym, 1) %in% c(0.6)) %>% 
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
    x = "Modularity",
    y = "Ratio of total infections versus\ntotal number of contacts"
  )+
  geom_vline(xintercept = seq(1.5, length(seq(0, 1, 0.1)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


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
  filter(round(p_asym, 1) %in% c(0.6)) %>% 
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
  group_by(p_asym, modularity, label) %>%
  filter(median == min(median)) %>%
  ungroup() %>%
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
    x = "Modularity",
    y = "Ratio of total infections versus\ntotal number of contacts"
  )+
  geom_vline(xintercept = seq(1.5, length(seq(0, 1, 0.1)) - 0.5, by = 1),
             linetype = "dashed", color = "grey50")+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


jpeg("/Volumes/argon_home/dissertation/real_data_application/paper3/figures/fig8.jpeg",
     width = 6.1, height = 4.8, units = "in", res = 1000)
grid::grid.draw(fig8)
dev.off()
