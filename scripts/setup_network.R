
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

# # Count interactions by between actors by timestamp
# net_time <- cns %>% 
#   group_by(number_timestamp) %>% 
#   count(user_a,user_b)
# 
# # Compute density by timestamp
# net_time %>% 
#   group_by(number_timestamp) %>% 
#   summarise(n_int = sum(n),
#             den = sum(n)/((total_actors*(total_actors-1))*0.5)) %>%  # we only have lower diag so sum is size of edge set
#   # we only have lower diag so sum is size of edge set
#   .$den %>% range

# Collapse network by day ------------------------------------------------------

# collapse 5 min timestamps for intereactions into daily level
interactions_temp <- cns %>% 
  distinct(number_timestamp) %>%
  mutate(new_timestamp = (row_number() - 1) %/% 288 + 1) # each timestamp is 5 mins so 288 timestamps day, subtract 1 as we want groups to be closed on right interval
# table(cut(1:8064, breaks = c(0, seq(288, 288*4*7, 288)) ))

# tibble with total actor pair interactions by day
interactions <- interactions_temp %>% 
  inner_join(cns, ., by = "number_timestamp") %>% 
  rename(timestamp = new_timestamp) %>% 
  count(timestamp, user_a, user_b)

# Count interactions by between actors by timestamp
interactions %>% 
  group_by(timestamp) %>%  
  summarise(n_int = n(),
            den = n()/((total_actors*(total_actors-1))*0.5)) %>%  # we only have lower diag so sum is size of edge set
  # we only have lower diag so sum is size of edge set
  .$den %>% range

# Create nested list with first layer being timestamp and within each t
# the neighbors of actors with interactions at time t.
# Do this for only the last two weeks as the interactions during the last 
# two weeks is what we will use in ABM

int_and_neighbors_by_t <- parallel::mclapply((1:28)[(1:28) > 28/2], mc.preschedule = TRUE, mc.cores = 5,
                                     
                                     FUN = function(x){
                                       
                                       cns_el <- interactions %>% 
                                         filter(timestamp == x) %>% # subset to interaction in first w weeks
                                         group_by(user_a, user_b) %>% 
                                         summarise(n = sum(n), .groups = "drop")
                                       
                                       A <- Matrix::sparseMatrix(i = cns_el$user_a, 
                                                                 j = cns_el$user_b,
                                                                 x = cns_el$n,
                                                                 dims = c(total_actors, total_actors),
                                                                 symmetric = TRUE)
                                       
                                       rownames(A) <- 1:total_actors
                                       colnames(A) <- rownames(A) 
                                       n_total_t <- rowSums(A)
                                       check_isolates <- n_total_t == 0
                                       
                                       temp2 <- lapply(rownames(A)[!check_isolates], function(y){
                                         z <- A[y, ]
                                         z <- z[z!=0]
                                         return(z)
                                       })
                                       
                                       names(temp2) <- rownames(A)[!check_isolates]
                                       
                                       return(list(neighbors = temp2,
                                                   n_total_t = n_total_t[!check_isolates]))
                                       
                                     })

names(int_and_neighbors_by_t) <- (1:28)[(1:28) > 28/2]

#############
# Version 2 #
#############

interactions_last_2_weeks <- interactions %>% 
  filter(timestamp > 28/2) %>%
  nest(data = c(user_a, user_b, n)) %>% 
  mutate(n_total_t = map(data, ~{
    pivot_longer(.x, cols = c("user_a", "user_b"),
                 names_to = "label",
                 values_to = "actor") %>% 
      group_by(actor) %>% 
      summarise(n_total_t = sum(n), .groups = "drop")
  })) %>% 
  mutate(neighbors = map(timestamp,
                         \(y){
                           distinct(bind_rows(interactions %>% 
                                                filter(timestamp == y) %>% select(actor = user_a),
                                              interactions %>% 
                                                filter(timestamp == y) %>% select(actor = user_b))) %>% 
                             mutate(neighbors = map(actor, 
                                                    ~{
                                                      interactions %>% 
                                                        filter(timestamp == y) %>% 
                                                        filter(user_a == .x| user_b == .x )  %>% 
                                                        rowwise() %>% 
                                                        mutate(neighbors = setdiff(c(user_a, user_b), .x)) %>% 
                                                        ungroup() %>% 
                                                        distinct(neighbors, n) 
                                                    }))
                           
                         })) %>% 
  select(-data)

neighbors_by_t2 = interactions_last_2_weeks %>%
  select(timestamp, neighbors) %>%     
  split(.$timestamp) %>%           
  map(~{
    .x[["neighbors"]] %>% .[[1]] %>% 
      split(.$actor) %>%  
      map(~ {
        temp <- .x[["neighbors"]] %>% .[[1]]
        setNames(temp$n, nm = temp$neighbors)})
  })
  
n_total_by_t2 <- interactions_last_2_weeks %>% 
  select(timestamp, n_total_t) %>% 
  split(.$timestamp) %>%           
  map(~{
     temp <- .x[["n_total_t"]] %>% .[[1]] 
     setNames(temp$n_total_t, nm = temp$actor)
  })

all.equal(neighbors_by_t2, lapply(int_and_neighbors_by_t, function(x) x$neighbors))
all.equal(n_total_by_t2, lapply(int_and_neighbors_by_t, function(x) x$n_total_t))

save(int_and_neighbors_by_t, total_actors, file = "int_and_neighbors_by_t.RData")

# Create aggregated network for first 2 weeks to run LSHM ----------------------

# cns_el1 <- cns %>% # this aggregates interactions over 4 weeks
#   count(user_a,user_b)

cns_el <- interactions %>% 
  filter(timestamp <= 28/2) %>% # subset to interaction in first w weeks
  group_by(user_a, user_b) %>% 
  summarise(n = sum(n), .groups = "drop")

# cns_net <-  tbl_graph(nodes = 
#                         user_values %>% select(name),
#                       edges = cns_el,
#                       directed = FALSE)
# 
# B <-  cns_net %>% 
#   igraph::as_adjacency_matrix(attr = "n")

A <- Matrix::sparseMatrix(i = cns_el$user_a, 
                          j = cns_el$user_b,
                          x = cns_el$n,
                          dims = c(total_actors, total_actors),
                          symmetric = TRUE)
# all.equal(A, B, check.attributes = F)

# Check for isolates. These will be removed from LSHM approach so they will no
# be assigned into a cluster. Thus, in intervention if these actors are the seed
# can only isolated the individual actors.
which(rowSums(A) == 0)

saveRDS(A, file = "adjacency_mat_week_1_2.RDS")
