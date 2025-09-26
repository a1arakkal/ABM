
library(JANE)
setwd("/Users/atlan/dissertation/real_data_application/paper3/")

# Pull in aggregated contacts over the first 2 weeks
A <- readRDS("adjacency_mat_week_1_2.RDS")

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
         model = "RS")
})
future::plan(future::sequential)
gc()


future::plan(future::multisession, workers = 16)
time_cns_fits_ignore_noise <- system.time({
  cns_fits_ignore_noise = 
    JANE(A,
         D = 2,
         K = 2:10,
         family = "poisson",
         noise_weights = FALSE,
         seed = 03091962,
         model = "RS")
})
future::plan(future::sequential)
gc()


future::plan(future::multisession, workers = 16)
time_cns_fits_di <- system.time({
  cns_fits_di= 
    JANE(A = as.matrix((A>1.0)*1.0),
         D = 2,
         K = 2:10,
         family = "bernoulli",
         noise_weights = FALSE,
         seed = 03091962,
         model = "RS")
})
future::plan(future::sequential)
gc()

save(cns_fits, time_cns_fits, cns_fits_ignore_noise, time_cns_fits_ignore_noise, 
     time_cns_fits_di, cns_fits_di, A,
     file = "cns_week_1_2_res.RData")


## Assortativity coeff ---------------------------------------------------------

#cns_fits_ignore_noise
assortnet::assortment.discrete(cns_fits_ignore_noise$A,
                               cns_fits_ignore_noise$optimal_res$cluster_labels,
                               weighted=F)$r

#cns_fits_di
assortnet::assortment.discrete(cns_fits_di$A,
                               cns_fits_di$optimal_res$cluster_labels,
                               weighted=F)$r

#cns_fits
# Remove noise edges from cns_fits based on hard clustering rule
test_A <- cns_fits$A
test_A[summary(cns_fits)$Z_W[summary(cns_fits)$Z_W[,6]==2, 1:2], drop = F] <- 0
test_A[summary(cns_fits)$Z_W[summary(cns_fits)$Z_W[,6]==2, 2:1], drop = F] <- 0

assortnet::assortment.discrete(test_A,
                               cns_fits$optimal_res$cluster_labels,
                               weighted=F)$r
