
library(JANE)
setwd("/Volumes/argon_home/dissertation/real_data_application/paper3/")

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
     file = "cns_res.RData")

