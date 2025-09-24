library(janitor)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(Matrix)

setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")
load(file = "cns_week_1_2_res.RData")

library(JANE)
clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

# rm(list = setdiff(ls(), c("int_and_neighbors_by_t", "A",
#                           "total_actors",
#                           "clusters_accounting_for_noise",
#                           "clusters_ignore_noise",
#                           "clusters_dichotomize_g_1")))
colnames(A)


#random clusters same size as lshm clusters
non_isolates <- as.character(1:total_actors)[ as.character(1:total_actors) %in% names(clusters_accounting_for_noise) ]
clusters_random_lshm <- rmultinom(length(non_isolates), 1, prob = prop.table(table(clusters_accounting_for_noise)))
clusters_random_lshm <- apply(clusters_random_lshm , 2, which.max)
names(clusters_random_lshm) <- non_isolates

table(clusters_random_lshm)
table(clusters_accounting_for_noise)
length(unique(names(clusters_random_lshm)))
length(unique(names(clusters_accounting_for_noise)))
setdiff(names(clusters_random_lshm), names(clusters_accounting_for_noise))
setdiff(names(clusters_accounting_for_noise), names(clusters_random_lshm))

order_LSHM <- c(names(sort((clusters_accounting_for_noise))), setdiff(colnames(A), names(clusters_accounting_for_noise)))
order_random <- c(names(sort((clusters_random_lshm))), setdiff(colnames(A), names(clusters_random_lshm)))

B <- as.matrix((A>0)*1)
image(B[order_LSHM,order_LSHM])
image(B[order_random,order_random])

# Within cluster total weight
test <- 0
for (i in unique(clusters_accounting_for_noise)){
test = test + sum(A[names(which(clusters_accounting_for_noise == i)), names(which(clusters_accounting_for_noise == i))])
}

test2 <- 0
for (i in unique(clusters_random_lshm)){
  test2 = test2 + sum(A[names(which(clusters_random_lshm == i)), names(which(clusters_random_lshm == i))])
}

# Between cluster total weight
test <- 0
for (i in unique(clusters_accounting_for_noise)){
  test = test + sum(A[names(which(clusters_accounting_for_noise == i)), setdiff(colnames(A), names(which(clusters_accounting_for_noise == i)))])
}

test2 <- 0
for (i in unique(clusters_random_lshm)){
  test2 = test2 + sum(A[names(which(clusters_random_lshm == i)), setdiff(colnames(A), names(which(clusters_random_lshm == i)))])
}

library(Matrix)
test_A <- cns_fits$A
test_A[summary(cns_fits)$Z_W[summary(cns_fits)$Z_W[,6]==2, 1:2], drop = F] <- 0
test_A[summary(cns_fits)$Z_W[summary(cns_fits)$Z_W[,6]==2, 2:1], drop = F] <- 0
test_A <- cbind(test_A, B[rownames(test_A), setdiff(colnames(B), colnames(test_A))])
test_A <- rbind(test_A, B[setdiff(rownames(B), rownames(test_A)),colnames(test_A)])
test_A <- as.matrix(test_A)
#cns_fits_ignore_noise
assortnet::assortment.discrete(test_A,
                               c(clusters_accounting_for_noise, 
                                 seq(max(clusters_accounting_for_noise), max(clusters_accounting_for_noise)+length(setdiff(colnames(test_A), names(clusters_accounting_for_noise)))-1)),
                               weighted=F)$r

assortnet::assortment.discrete(cns_fits$A,
                               clusters_accounting_for_noise,
                               weighted=F)$r
image(test_A[order,order])



#cns_fits_di
assortnet::assortment.discrete(cns_fits$A,
                               clusters_random_lshm,
                               weighted=F)$r


# How well does it do with last 2 weeks of data --------------------------------

setwd("/Volumes/argon_home/dissertation/real_data_application/paper2/")

# https://figshare.com/articles/dataset/The_Copenhagen_Networks_Study_interaction_data/7267433/1
cns = 
  readr::read_csv("bt_symmetric.csv")

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

t= 16
cns_el <- interactions %>% 
  filter(timestamp == t) %>% # subset to interaction in last 2 weeks
  group_by(user_a, user_b) %>% 
  summarise(n = sum(n), .groups = "drop")


A <- Matrix::sparseMatrix(i = cns_el$user_a, 
                          j = cns_el$user_b,
                          x = cns_el$n,
                          dims = c(total_actors, total_actors),
                          symmetric = TRUE)
# all.equal(A, B, check.attributes = F)
rownames(A) <- 1:total_actors
colnames(A) <- 1:total_actors


order_LSHM <- c(names(sort((clusters_accounting_for_noise))), setdiff(colnames(A), names(clusters_accounting_for_noise)))
order_random <- c(names(sort((clusters_random_lshm))), setdiff(colnames(A), names(clusters_random_lshm)))
B <- as.matrix((A>0)*1)
image(B[order_LSHM,order_LSHM])
image(B[order_random,order_random])

# Within cluster total weight
test <- 0
for (i in unique(clusters_accounting_for_noise)){
  test = test + sum(A[names(which(clusters_accounting_for_noise == i)), names(which(clusters_accounting_for_noise == i))])
}
test

test2 <- 0
for (i in unique(clusters_random_lshm)){
  test2 = test2 + sum(A[names(which(clusters_random_lshm == i)), names(which(clusters_random_lshm == i))])
}
test2

# Between cluster total weight
test <- 0
for (i in unique(clusters_accounting_for_noise)){
  test = test + sum(A[names(which(clusters_accounting_for_noise == i)), setdiff(colnames(A), names(which(clusters_accounting_for_noise == i)))])
}
test

test2 <- 0
for (i in unique(clusters_random_lshm)){
  test2 = test2 + sum(A[names(which(clusters_random_lshm == i)), setdiff(colnames(A), names(which(clusters_random_lshm == i)))])
}
test2


int_and_neighbors_t <- int_and_neighbors_by_t[[as.character(t)]]
symptom_infect_t = "16"
clusters <- clusters_accounting_for_noise
# clusters <- clusters_random_lshm

## Find those that are in a cluster
symptom_infect_t_in_cluster <- symptom_infect_t[symptom_infect_t %in% names(clusters)]

## Find those that not in a cluster as they are isolates removed before running LSHM
symptom_infect_t_not_in_cluster <- setdiff(symptom_infect_t, symptom_infect_t_in_cluster)

## Get the cluster lables for those in clusters
cluster_infected <- unique(unname(clusters[symptom_infect_t_in_cluster]))

## Get all the actors in those clusters
quarantine_in_cluster <- names(clusters[clusters %in% cluster_infected])

## Quarantine actors include both those in an infected cluster and isolates
qua_t <- unique(c(quarantine_in_cluster,
                       symptom_infect_t_not_in_cluster)) # deals with isolates that are removed from LSHM implementation


int_and_neighbors_t$neighbors <- int_and_neighbors_t$neighbors[setdiff(names(int_and_neighbors_t$neighbors), qua_t)] 

if(length(int_and_neighbors_t$neighbors) > 0){
  
  temp <- lapply(int_and_neighbors_t$neighbors,
                 FUN = function(x){
                   p1 <- x[setdiff(names(x), qua_t)] # remove the interactions with quarantined actors
                   p2 <- sum(p1) # total interactions excluding quarantined actors
                   return(list(p1 = p1,
                               p2 = p2))
                 })
  
  int_and_neighbors_t$neighbors <- lapply(temp,
                                          FUN = function(x){x$p1})
  
  int_and_neighbors_t$n_total_t <- sapply(temp,
                                          FUN = function(x){x$p2})
  
}

C <- A
C[qua_t, ] <- 0
C[ , qua_t ] <- 0

all.equal(colSums(C)[colSums(C)>0],   int_and_neighbors_t$n_total_t[ int_and_neighbors_t$n_total_t>0])
all.equal(rowSums(C)[colSums(C)>0],   int_and_neighbors_t$n_total_t[ int_and_neighbors_t$n_total_t>0])

all.equal(C["5", ][C["5", ]>0],   int_and_neighbors_t$neighbors[["5"]][int_and_neighbors_t$neighbors[["5"]]>0])
sum(C)
sum(C[C>1])
order<- c(names(sort((clusters))), setdiff(colnames(A), names(clusters)))
image(as.matrix((C>0)*1)[order,order])