
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
  readr::read_csv("fb_friends.csv")

cns %<>%
  clean_names()

total_actors <- nrow(user_values)

cns %<>% 
  rename(user_a = number_user_a) %>% 
  filter(user_a != user_b) # 11 cases where they are friends with themselves???

# Replace user values with new user values
cns %<>%
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

# cns_net <-  tbl_graph(nodes =
#                         user_values %>% select(name),
#                       edges = cns,
#                       directed = F)
# 
# B <-  cns_net %>%
#   igraph::as_adjacency_matrix(attr = "value")

A <- Matrix::sparseMatrix(i = cns$user_a, 
                          j = cns$user_b,
                          x = cns$value,
                          dims = c(total_actors, total_actors),
                          symmetric = T)
# all.equal(as.matrix(A), as.matrix(B), check.attributes = F)
rownames(A) <- 1:total_actors
colnames(A) <- 1:total_actors

# Check for isolates. These will be removed from LSHM approach so they will no
# be assigned into a cluster. Thus, in intervention if these actors are the seed
# can only isolated the individual actors.
which(rowSums(A) == 0) #36

saveRDS(A, file = "fb_network.RDS")


setwd("/Users/atlan/dissertation/real_data_application/paper3/")
A <- readRDS("fb_network.RDS")

library(JANE)
future::plan(future::multisession, workers = 16)
time_cns_fits <- system.time({
  cns_fits = 
    JANE(A,
         D = 2,
         K = 2:10,
         noise_weights = FALSE,
         seed = 03091962,
         model = "RS")
})
future::plan(future::sequential)
gc()

save(cns_fits, time_cns_fits, A,
     file = "fb_net_LSM_res.RData")
