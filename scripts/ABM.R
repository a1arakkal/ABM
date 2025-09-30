# Load in interactions by time list --------------------------------------------

setwd("/Users/atlan/dissertation/real_data_application/paper3/")

load(file = "int_and_neighbors_by_t.RData")
int_and_neighbors_by_t_true <- int_and_neighbors_by_t

# Load in LSHM results --------------------------------------------

load(file = "cns_week_1_2_res.RData")

clusters_accounting_for_noise <- JANE:::summary.JANE(cns_fits)$cluster_labels
clusters_ignore_noise <- JANE:::summary.JANE(cns_fits_ignore_noise)$cluster_labels
clusters_dichotomize_g_1 <- JANE:::summary.JANE(cns_fits_di)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise",
                          "clusters_ignore_noise",
                          "clusters_dichotomize_g_1")))

load(file = "fb_net_LSM_res.RData")
clusters_fb_network <- JANE:::summary.JANE(cns_fits)$cluster_labels

rm(list = setdiff(ls(), c("int_and_neighbors_by_t_true",
                          "total_actors",
                          "clusters_accounting_for_noise",
                          "clusters_ignore_noise",
                          "clusters_dichotomize_g_1",
                          "clusters_fb_network")))

# Source helper functions ------------------------------------------------------

source("https://raw.githubusercontent.com/a1arakkal/ABM/refs/heads/master/scripts/helper_functions.R")

# Parameters ------------------------------------------------------------------

# prob of transmisson among infected
p_infected <- 0.0025

# mean days exposed
mean_exposure_days <- NULL # NULL is SIR model

# days infected
mean_infected_days <- c("asymptomatic" = 10L,
                        "pre_symptomatic" = 5L, # if 0 not pre_symptomatic period
                        "symptomatic" = 5L)

# prob of being asymptomatic given infected
p_asym_seq <- seq(0, 1, by = 0.2) # if 0 all infections are symptomatic

# Initilize actors
actor_labels <- 1:total_actors

# Clusters for quarantine
# clusters <- clusters_accounting_for_noise # NULL no intervention, any unnamed vector of length 1 will isolate infected and their 1-hop neighbors
quarantine_days_seq  <- c(5, 10, 15) # ignored if cluster is NULL (i.e., no intervention)

# N time-steps
timesteps <- names(int_and_neighbors_by_t_true)

# Number of times to repeat weeks of interaction
n_repeat <- 4 # 2*n_repeat weeks

# Minimum degree for seed options if 0 will randomly select from 1-692
min_degree_t1 <- 0

# How far to look back for digital contact tracing if used as intervention, inclusive of current day
digital_contact_tracing_look_back <- 4

# Sensitivity and specificity of digital contact tracing if used as intervention (allows for false positives and false negatives)
DCT_sensitivity_seq  <- c(0, seq(.6, 1, by = .1))
DCT_specificity_seq  <- c(0, seq(.6, 1, by = .1))

# test <- run_single_ABM(p_infected = p_infected,
#                        mean_exposure_days = mean_exposure_days,
#                        mean_infected_days = mean_infected_days,
#                        actor_labels = actor_labels,
#                        min_degree_t1 = min_degree_t1,
#                        timesteps = timesteps,
#                        n_repeat = n_repeat,
#                        p_asym = p_asym,
#                        clusters = clusters,
#                        quarantine_days = quarantine_days,
#                        digital_contact_tracing_look_back = digital_contact_tracing_look_back,
#                        DCT_sensitivity = DCT_sensitivity,
#                        DCT_specificity = DCT_specificity,
#                        int_and_neighbors_by_t = int_and_neighbors_by_t_true)
# 
# microbenchmark::microbenchmark(run_single_ABM(p_infected = p_infected,
#                                               mean_exposure_days = mean_exposure_days,
#                                               mean_infected_days = mean_infected_days,
#                                               actor_labels = actor_labels,
#                                               min_degree_t1 = min_degree_t1,
#                                               timesteps = timesteps,
#                                               n_repeat = n_repeat,
#                                               p_asym = p_asym,
#                                               clusters = clusters,
#                                               quarantine_days = quarantine_days,
#                                               digital_contact_tracing_look_back = digital_contact_tracing_look_back,
#                                               DCT_sensitivity = DCT_sensitivity,
#                                               DCT_specificity = DCT_specificity,
#                                               int_and_neighbors_by_t = int_and_neighbors_by_t_true),
#                                times = 100)

# Run ABM for multiple trials --------------------------------------------------

n_trial <- 1e3
cores <- 125L

set.seed(1234, kind = "L'Ecuyer-CMRG")
seeds <- sample.int(1e6, size = n_trial, replace = FALSE)

for (p_asym in p_asym_seq){
  for (quarantine_days in quarantine_days_seq){
    for (DCT_sensitivity in DCT_sensitivity_seq){
      for (DCT_specificity in DCT_specificity_seq){
        
        ## Run ABM with no intervention
        run_ABM_for_R0 <- parallel::mclapply(1:n_trial,
                                             FUN = function(x){
                                               
                                               set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                               
                                               run_single_ABM(p_infected = p_infected,
                                                              mean_exposure_days = 1000, # makes it so seed is the only infective for duration of ABM
                                                              mean_infected_days = mean_infected_days,
                                                              actor_labels = actor_labels,
                                                              min_degree_t1 = min_degree_t1,
                                                              timesteps = timesteps,
                                                              n_repeat = n_repeat,
                                                              p_asym = p_asym,
                                                              clusters = NULL,
                                                              quarantine_days = quarantine_days,
                                                              digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                              DCT_sensitivity = DCT_sensitivity,
                                                              DCT_specificity = DCT_specificity,
                                                              int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                             mc.preschedule = TRUE,
                                             mc.cores = cores)
        
        R0_est_only_seed_infected <- mean(sapply(run_ABM_for_R0, function(x){x$R0}))
        
        ## Run ABM with no intervention
        run_ABM_no_intervention <- parallel::mclapply(1:n_trial,
                                                      FUN = function(x){
                                                        
                                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                        
                                                        run_single_ABM(p_infected = p_infected,
                                                                       mean_exposure_days = mean_exposure_days,
                                                                       mean_infected_days = mean_infected_days,
                                                                       actor_labels = actor_labels,
                                                                       min_degree_t1 = min_degree_t1,
                                                                       timesteps = timesteps,
                                                                       n_repeat = n_repeat,
                                                                       p_asym = p_asym,
                                                                       clusters = NULL,
                                                                       quarantine_days = quarantine_days,
                                                                       digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                       DCT_sensitivity = DCT_sensitivity,
                                                                       DCT_specificity = DCT_specificity,
                                                                       int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                      mc.preschedule = TRUE,
                                                      mc.cores = cores)
        
        ## Run ABM with isolating individual actors not clusters
        run_ABM_isolate_individuals <- parallel::mclapply(1:n_trial,
                                                          FUN = function(x){
                                                            
                                                            set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                            
                                                            run_single_ABM(p_infected = p_infected,
                                                                           mean_exposure_days = mean_exposure_days,
                                                                           mean_infected_days = mean_infected_days,
                                                                           actor_labels = actor_labels,
                                                                           min_degree_t1 = min_degree_t1,
                                                                           timesteps = timesteps,
                                                                           n_repeat = n_repeat,
                                                                           p_asym = p_asym,
                                                                           clusters = 1, # if vector of length 1 will use digital contact tracing approach
                                                                           quarantine_days = quarantine_days,
                                                                           digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                           DCT_sensitivity = DCT_sensitivity,
                                                                           DCT_specificity = DCT_specificity,
                                                                           int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                          mc.preschedule = TRUE,
                                                          mc.cores = cores)
        
        ## Run ABM with clustering accounting for noise
        run_ABM_accounting_for_noise <- parallel::mclapply(1:n_trial,
                                                           FUN = function(x){
                                                             
                                                             set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                             
                                                             run_single_ABM(p_infected = p_infected,
                                                                            mean_exposure_days = mean_exposure_days,
                                                                            mean_infected_days = mean_infected_days,
                                                                            actor_labels = actor_labels,
                                                                            min_degree_t1 = min_degree_t1,
                                                                            timesteps = timesteps,
                                                                            n_repeat = n_repeat,
                                                                            p_asym = p_asym,
                                                                            clusters = clusters_accounting_for_noise,
                                                                            quarantine_days = quarantine_days,
                                                                            digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                            DCT_sensitivity = DCT_sensitivity,
                                                                            DCT_specificity = DCT_specificity,
                                                                            int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                           mc.preschedule = TRUE,
                                                           mc.cores = cores)
        
        ## Run ABM with random clusters similar size to clustering accounting for noise (allow different assignment for each run)
        run_ABM_random_lshm <- parallel::mclapply(1:n_trial,
                                                  FUN = function(x){
                                                    
                                                    set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                    
                                                    #random clusters same size as lshm clusters
                                                    clusters_random_lshm <- clusters_accounting_for_noise
                                                    permute <- sample(names(clusters_random_lshm), size = length(clusters_random_lshm), replace = FALSE)
                                                    names(clusters_random_lshm) <- permute
                                                    
                                                    run_single_ABM(p_infected = p_infected,
                                                                   mean_exposure_days = mean_exposure_days,
                                                                   mean_infected_days = mean_infected_days,
                                                                   actor_labels = actor_labels,
                                                                   min_degree_t1 = min_degree_t1,
                                                                   timesteps = timesteps,
                                                                   n_repeat = n_repeat,
                                                                   p_asym = p_asym,
                                                                   clusters = clusters_random_lshm,
                                                                   quarantine_days = quarantine_days,
                                                                   digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                   DCT_sensitivity = DCT_sensitivity,
                                                                   DCT_specificity = DCT_specificity,
                                                                   int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                  
                                                  mc.preschedule = TRUE,
                                                  mc.cores = cores)
        
        ## Run ABM with random clusters similar size to clustering accounting for noise (require same assignment for each run)
        #random clusters same size as lshm clusters
        set.seed(1234, kind = "L'Ecuyer-CMRG")
        clusters_random_lshm2 <- clusters_accounting_for_noise
        permute <- sample(names(clusters_random_lshm2), size = length(clusters_random_lshm2), replace = FALSE)
        names(clusters_random_lshm2) <- permute
        
        run_ABM_random_lshm_outside <- parallel::mclapply(1:n_trial,
                                                          FUN = function(x){
                                                            
                                                            set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                            
                                                            run_single_ABM(p_infected = p_infected,
                                                                           mean_exposure_days = mean_exposure_days, 
                                                                           mean_infected_days = mean_infected_days,
                                                                           actor_labels = actor_labels,
                                                                           min_degree_t1 = min_degree_t1,
                                                                           timesteps = timesteps,
                                                                           n_repeat = n_repeat,
                                                                           p_asym = p_asym,
                                                                           clusters = clusters_random_lshm2,
                                                                           quarantine_days = quarantine_days,
                                                                           digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                           DCT_sensitivity = DCT_sensitivity,
                                                                           DCT_specificity = DCT_specificity,
                                                                           int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                          
                                                          mc.preschedule = TRUE,
                                                          mc.cores = cores)
        
        ## Run ABM with random clusters size and number of clusters
        run_ABM_random_lshm_size_and_number <- parallel::mclapply(1:n_trial,
                                                                  FUN = function(x){
                                                                    
                                                                    set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                                    
                                                                    n_clusters <- sample(2:10, size = 1, replace = F)
                                                                    p_clust <- extraDistr::rdirichlet(1, rep(3, n_clusters))
                                                                    non_isolates <- names(clusters_accounting_for_noise)
                                                                    clusters_random_lshm <- rmultinom(length(non_isolates), 1, prob = p_clust)
                                                                    clusters_random_lshm <- apply(clusters_random_lshm , 2, which.max)
                                                                    names(clusters_random_lshm) <- non_isolates
                                                                    
                                                                    run_single_ABM(p_infected = p_infected,
                                                                                   mean_exposure_days = mean_exposure_days, 
                                                                                   mean_infected_days = mean_infected_days,
                                                                                   actor_labels = actor_labels,
                                                                                   min_degree_t1 = min_degree_t1,
                                                                                   timesteps = timesteps,
                                                                                   n_repeat = n_repeat,
                                                                                   p_asym = p_asym,
                                                                                   clusters = clusters_random_lshm,
                                                                                   quarantine_days = quarantine_days,
                                                                                   digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                                   DCT_sensitivity = DCT_sensitivity,
                                                                                   DCT_specificity = DCT_specificity,
                                                                                   int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                                  
                                                                  mc.preschedule = TRUE,
                                                                  mc.cores = cores)
        
        ## Run ABM with clustering ignoring for noise
        run_ABM_ignore_noise <- parallel::mclapply(1:n_trial,
                                                   FUN = function(x){
                                                     
                                                     set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                     
                                                     run_single_ABM(p_infected = p_infected,
                                                                    mean_exposure_days = mean_exposure_days,
                                                                    mean_infected_days = mean_infected_days,
                                                                    actor_labels = actor_labels,
                                                                    min_degree_t1 = min_degree_t1,
                                                                    timesteps = timesteps,
                                                                    n_repeat = n_repeat,
                                                                    p_asym = p_asym,
                                                                    clusters = clusters_ignore_noise,
                                                                    quarantine_days = quarantine_days,
                                                                    digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                    DCT_sensitivity = DCT_sensitivity,
                                                                    DCT_specificity = DCT_specificity,
                                                                    int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                   mc.preschedule = TRUE,
                                                   mc.cores = cores)
        
        ## Run ABM with clustering dichotomize network with cutoff of 1
        run_ABM_dichotomize_g_1 <- parallel::mclapply(1:n_trial,
                                                      FUN = function(x){
                                                        
                                                        set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                        
                                                        run_single_ABM(p_infected = p_infected,
                                                                       mean_exposure_days = mean_exposure_days,
                                                                       mean_infected_days = mean_infected_days,
                                                                       actor_labels = actor_labels,
                                                                       min_degree_t1 = min_degree_t1,
                                                                       timesteps = timesteps,
                                                                       n_repeat = n_repeat,
                                                                       p_asym = p_asym,
                                                                       clusters = clusters_dichotomize_g_1,
                                                                       quarantine_days = quarantine_days,
                                                                       digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                       DCT_sensitivity = DCT_sensitivity,
                                                                       DCT_specificity = DCT_specificity,
                                                                       int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                      mc.preschedule = TRUE,
                                                      mc.cores = cores)
        
        ## Run ABM with clustering based on FB data
        run_ABM_fb_clusters <- parallel::mclapply(1:n_trial,
                                                  FUN = function(x){
                                                    
                                                    set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                    
                                                    run_single_ABM(p_infected = p_infected,
                                                                   mean_exposure_days = mean_exposure_days,
                                                                   mean_infected_days = mean_infected_days,
                                                                   actor_labels = actor_labels,
                                                                   min_degree_t1 = min_degree_t1,
                                                                   timesteps = timesteps,
                                                                   n_repeat = n_repeat,
                                                                   p_asym = p_asym,
                                                                   clusters = clusters_fb_network,
                                                                   quarantine_days = quarantine_days,
                                                                   digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                   DCT_sensitivity = DCT_sensitivity,
                                                                   DCT_specificity = DCT_specificity,
                                                                   int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                  mc.preschedule = TRUE,
                                                  mc.cores = cores)
        
        ## Run ABM with clustering on random partition similar size as FB clusters
        run_ABM_random_fb <- parallel::mclapply(1:n_trial,
                                                FUN = function(x){
                                                  
                                                  set.seed(seeds[x], kind = "L'Ecuyer-CMRG")
                                                  
                                                  #random clusters same size as fb clusters
                                                  clusters_random_fb <- clusters_fb_network
                                                  permute <- sample(names(clusters_random_fb), size = length(clusters_random_fb), replace = FALSE)
                                                  names(clusters_random_fb) <- permute
                                                  
                                                  run_single_ABM(p_infected = p_infected,
                                                                 mean_exposure_days = mean_exposure_days,
                                                                 mean_infected_days = mean_infected_days,
                                                                 actor_labels = actor_labels,
                                                                 min_degree_t1 = min_degree_t1,
                                                                 timesteps = timesteps,
                                                                 n_repeat = n_repeat,
                                                                 p_asym = p_asym,
                                                                 clusters = clusters_random_fb,
                                                                 quarantine_days = quarantine_days,
                                                                 digital_contact_tracing_look_back = digital_contact_tracing_look_back,
                                                                 DCT_sensitivity = DCT_sensitivity,
                                                                 DCT_specificity = DCT_specificity,
                                                                 int_and_neighbors_by_t = int_and_neighbors_by_t_true)},
                                                
                                                mc.preschedule = TRUE,
                                                mc.cores = cores)
        
        res <- list(run_ABM_no_intervention = run_ABM_no_intervention,
                    run_ABM_isolate_individuals = run_ABM_isolate_individuals,
                    run_ABM_accounting_for_noise = run_ABM_accounting_for_noise,
                    run_ABM_random_lshm = run_ABM_random_lshm,
                    run_ABM_random_lshm_outside = run_ABM_random_lshm_outside,
                    run_ABM_random_lshm_size_and_number = run_ABM_random_lshm_size_and_number,
                    run_ABM_ignore_noise = run_ABM_ignore_noise,
                    run_ABM_dichotomize_g_1 = run_ABM_dichotomize_g_1,
                    run_ABM_fb_clusters = run_ABM_fb_clusters,
                    run_ABM_random_fb = run_ABM_random_fb,
                    R0_est_only_seed_infected = R0_est_only_seed_infected,
                    quarantine_days = quarantine_days,
                    p_asym = p_asym,
                    DCT_sensitivity = DCT_sensitivity,
                    DCT_specificity = DCT_specificity)
        
        #note res_out folder assumes prob false postive = 1-p_sen same for all non-ints.
        
        saveRDS(res,
                file = paste0("res_out_false_pos_degree_diff_random_clust/qdays_", quarantine_days, "_pasym_", p_asym, "_DTCspec_", DCT_specificity, "_DTCsen_", DCT_sensitivity, ".RDS"))
        
        print(paste0("qdays_", quarantine_days, "_pasym_", p_asym, "_DTCspec_", DCT_specificity, "_DTCsen_", DCT_sensitivity, " saved"))
        # save(run_ABM_no_intervention, run_ABM_isolate_individuals,
        #      run_ABM_accounting_for_noise, run_ABM_random_lshm,
        #      run_ABM_ignore_noise, run_ABM_dichotomize_g_1,
        #      run_ABM_fb_clusters, run_ABM_random_fb,
        #      R0_est_only_seed_infected,
        #      file = "test.run.RData")
        
        # res <- readRDS(file = paste0("res_out_false_pos_degree/qdays_", quarantine_days, "_pasym_", p_asym, "_DTCspec_", DCT_specificity, "_DTCsen_", DCT_sensitivity, ".RDS"))
        # 
        # res$run_ABM_random_lshm_outside <- run_ABM_random_lshm_outside
        # res$run_ABM_random_lshm_size_and_number <- run_ABM_random_lshm_size_and_number
        # res$run_ABM_random_lshm <- run_ABM_random_lshm
        # res$run_ABM_random_fb <- run_ABM_random_fb
        # 
        # saveRDS(res,
        #         file = paste0("res_out_false_pos_degree_diff_random_clust/qdays_", quarantine_days, "_pasym_", p_asym, "_DTCspec_", DCT_specificity, "_DTCsen_", DCT_sensitivity, ".RDS"))
        # 
        
      }
    }
  }
}