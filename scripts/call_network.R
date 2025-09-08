
library(janitor)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(Matrix)

# Pull in interaction data ------------------------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper2/")

# https://figshare.com/articles/dataset/The_Copenhagen_Networks_Study_interaction_data/7267433/1
cns = 
  readr::read_csv("bt_symmetric.csv")

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

cns %<>%
  clean_names()

cns %<>%
  filter(user_b >= 0) # See readme.  Negative values are not meaningful.

# Unique user values
user_values <- tibble(user = sort(unique(c(cns$user_a,cns$user_b)))) %>% 
  mutate(name = 1:n())

rm(cns)

cns = 
  readr::read_csv("calls.csv")

cns %<>%
  clean_names()

cns %<>%
  filter(duration >= 0) # See readme.  -1 means missed call

total_actors <- nrow(user_values)

# Replace user values with new user values
cns %<>%
  filter( (caller %in% user_values$user) &  (callee %in% user_values$user)) %>% 
  left_join(user_values,
            by = join_by(caller == user)) %>% 
  select(-caller) %>% 
  rename(caller = name) %>% 
  left_join(user_values,
            by = join_by(callee == user)) %>% 
  select(-callee) %>% 
  rename(callee = name)

cns_aggregated <- cns %>% 
  group_by(caller, callee) %>% 
  summarise(duration = sum(duration), .groups = "drop") %>% 
  filter(duration>0)

# cns_net <-  tbl_graph(nodes =
#                         user_values %>% select(name),
#                       edges = cns_aggregated,
#                       directed = T)
# 
# B <-  cns_net %>%
#   igraph::as_adjacency_matrix(attr = "duration")

A <- Matrix::sparseMatrix(i = cns_aggregated$caller, 
                          j = cns_aggregated$callee,
                          x = cns_aggregated$duration,
                          dims = c(total_actors, total_actors),
                          symmetric = FALSE)
# all.equal(A, B, check.attributes = F)
rownames(A) <- 1:total_actors
colnames(A) <- 1:total_actors

# Check for isolates. These will be removed from LSHM approach so they will no
# be assigned into a cluster. Thus, in intervention if these actors are the seed
# can only isolated the individual actors.
which(rowSums(A) == 0 & colSums(A) == 0) #236

saveRDS(A, file = "call_network.RDS")


setwd("/Users/atlan/dissertation/real_data_application/paper3/")
A <- readRDS("call_network.RDS")

library(JANE)
future::plan(future::multisession, workers = 16)
time_cns_fits <- system.time({
  cns_fits = 
    JANE(A,
         D = 2,
         K = 2:10,
         family = "poisson",
         noise_weights = TRUE,
         guess_noise_weights = 1.0,
         seed = 03091962,
         model = "RSR")
})
future::plan(future::sequential)
gc()

save(cns_fits, time_cns_fits, A,
     file = "call_net_LSHM_res.RData")
