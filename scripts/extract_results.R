
# Load libraries ---------------------------------------------------------------

library(tidyverse)

# Helper function --------------------------------------------------------------

box_plot_fun <- function(y, name){
  bp_stats <- boxplot.stats(do.call("c", lapply(y, function(x){max(x[[name]])})))
  tibble(
    ymin = bp_stats$stats[1],  # lower whisker
    Q1   = bp_stats$stats[2],  # first quartile
    median = bp_stats$stats[3], # median
    Q3   = bp_stats$stats[4],  # third quartile
    ymax = bp_stats$stats[5]  # upper whisker
  ) %>% 
    mutate(metric = name)
}

extract_function <- function(run_ABM){
  
  n_trial <- length(run_ABM)
  
  out <- tibble()
  
  for (i in c("R0", "attack_rate", 
              "quarantined_rate",
              "average_time_in_quarantine", 
              "average_cumulative_interactions_per_actor")){
    
    out <- bind_rows(out, box_plot_fun(run_ABM, i))
    
  }
    
  # Max quarantined across trials ----------------------------------------------
  quarantined_max_stats <- boxplot.stats(do.call("c", lapply(run_ABM, function(x){max(x$quarantined)})))
  quarantined_max <- tibble(
    ymin = quarantined_max_stats$stats[1],  # lower whisker
    Q1   = quarantined_max_stats$stats[2],  # first quartile
    median = quarantined_max_stats$stats[3], # median
    Q3   = quarantined_max_stats$stats[4],  # third quartile
    ymax = quarantined_max_stats$stats[5]  # upper whisker
  ) %>% 
    mutate(metric = "quarantined_max")
  
  out <- bind_rows(out, quarantined_max)
  
  # Quarantined by time across trials ----------------------------------------------
  quarantined_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$quarantined}))
  quarantined_time_stats <- apply(quarantined_time_stats, 1, FUN = function(x){boxplot.stats(x)$stats})
  quarantined_time <- as_tibble(t(quarantined_time_stats))
  colnames(quarantined_time) <- c("ymin", "Q1", "median", "Q3", "ymax")
  quarantined_time <- quarantined_time%>% 
    mutate(metric = "quarantined_time",
           t = row_number())
  
  out <- bind_rows(out, quarantined_time)
  
  # Max incidence across trials ----------------------------------------------
  incidence_max_stats <- boxplot.stats(do.call("c", lapply(run_ABM, function(x){max(x$incidence)})))
  incidence_max <- tibble(
    ymin = incidence_max_stats$stats[1],  # lower whisker
    Q1   = incidence_max_stats$stats[2],  # first quartile
    median = incidence_max_stats$stats[3], # median
    Q3   = incidence_max_stats$stats[4],  # third quartile
    ymax = incidence_max_stats$stats[5]  # upper whisker
  ) %>% 
    mutate(metric = "incidence_max")
  
  out <- bind_rows(out, incidence_max)
  
  # Incidence by time across trials ----------------------------------------------
  incidence_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$incidence}))
  incidence_time_stats <- apply(incidence_time_stats, 1, FUN = function(x){boxplot.stats(x)$stats})
  incidence_time <- as_tibble(t(incidence_time_stats))
  colnames(incidence_time) <- c("ymin", "Q1", "median", "Q3", "ymax")
  incidence_time <- incidence_time%>% 
    mutate(metric = "incidence_time",
           t = row_number())
  
  out <- bind_rows(out, incidence_time)
  
  # Loss function 1: total infections/total person days outside quarnatine
  incidence_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$incidence}))
  incidence_total <- apply(incidence_time_stats, 2, FUN = sum)
  
  quarantined_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$quarantined}))
  quarantined_total <- apply(quarantined_time_stats, 2, FUN = sum)
  person_time_outside_quarantine <- (692*nrow(quarantined_time_stats))-quarantined_total
  
  tot_inf_by_person_days_outside_qua <- boxplot.stats(incidence_total/person_time_outside_quarantine)$stats
  tot_inf_by_person_days_outside_qua <- as_tibble(t(tot_inf_by_person_days_outside_qua))
  colnames(tot_inf_by_person_days_outside_qua) <- c("ymin", "Q1", "median", "Q3", "ymax")
  tot_inf_by_person_days_outside_qua <- tot_inf_by_person_days_outside_qua%>% 
    mutate(metric = "tot_inf_by_person_days_outside_qua")
  
  out <- bind_rows(out, tot_inf_by_person_days_outside_qua)
  
  # Loss function 2: total infections/total number of contacts 
  
  contacts_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$average_interaction*692})) # average_interaction[k] <- sum(n_interactions)/length(actors) 
  contacts_total <- apply(contacts_time_stats, 2, FUN = sum)
  
  tot_inf_by_by_contacts <- boxplot.stats(incidence_total/contacts_total)$stats
  tot_inf_by_by_contacts <- as_tibble(t(tot_inf_by_by_contacts))
  colnames(tot_inf_by_by_contacts) <- c("ymin", "Q1", "median", "Q3", "ymax")
  tot_inf_by_by_contacts <- tot_inf_by_by_contacts%>% 
    mutate(metric = "tot_inf_by_by_contacts")
  
  out <- bind_rows(out, tot_inf_by_by_contacts)
  
  return(out)
  
}

# Extract Results --------------------------------------------------------------

# setwd("/Users/atlan/dissertation/real_data_application/paper3/res_out")
setwd("/Users/atlan/dissertation/real_data_application/paper3/res_out_false_pos_degree")
files <- list.files()
files <- setdiff(files, "main_res.RData")

outer_fun <- function(file_name){
  
  tryCatch({
    data_list <- readRDS(file_name)
    
    out <- tibble(type = names(data_list)[stringr::str_detect(names(data_list), pattern = "run_ABM")]) %>% 
      mutate(res = map(type, ~extract_function(data_list[[.]]))) %>% 
      unnest(res) %>% 
      mutate(quarantine_days = data_list[["quarantine_days"]],
             p_asym = data_list[["p_asym"]],
             DCT_sensitivity = data_list[["DCT_sensitivity"]],
             DCT_specificity = data_list[["DCT_specificity"]],
             R0_est_only_seed_infected = data_list[["R0_est_only_seed_infected"]])
    
    rm(data_list)
    
    return(out)
  }, 
  error = function(e){
    cat("Error with", file_name)
    return(tibble(type = NA,
                  ymin = NA,
                  Q1 = NA,
                  median = NA,
                  Q3 = NA,
                  ymax = NA,
                  metric = NA,
                  t = NA, 
                  quarantine_days = NA,
                  p_asym = NA,
                  DCT_sensitivity = NA,
                  DCT_specificity = NA,
                  R0_est_only_seed_infected = NA))
  }
  )
}

main_res <- parallel::mclapply(files,
                               FUN = function(x){outer_fun(x)},
                               mc.cores = 20)

main_res <- do.call("bind_rows", main_res)

save(main_res, file = "main_res.RData")

