
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

extract_function <- function(run_ABM, baseline){
  
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
  
  tot_inf <- boxplot.stats(incidence_total)$stats
  tot_inf <- as_tibble(t(tot_inf))
  colnames(tot_inf) <- c("ymin", "Q1", "median", "Q3", "ymax")
  tot_inf <- tot_inf%>% 
    mutate(metric = "total_inf")
  
  quarantined_total <- boxplot.stats(quarantined_total)$stats
  quarantined_total <- as_tibble(t(quarantined_total))
  colnames(quarantined_total) <- c("ymin", "Q1", "median", "Q3", "ymax")
  quarantined_total <- quarantined_total%>% 
    mutate(metric = "quarantined_total_person_days")
  
  out <- bind_rows(out, 
                   tot_inf_by_person_days_outside_qua,
                   tot_inf,
                   quarantined_total)
  
  # Loss function 2: total infections/total number of contacts 
  
  contacts_time_stats <- do.call("cbind", lapply(run_ABM, function(x){x$average_interaction*692})) # average_interaction[k] <- sum(n_interactions)/length(actors) 
  contacts_total <- apply(contacts_time_stats, 2, FUN = sum)
  
  tot_inf_by_by_contacts <- boxplot.stats(incidence_total/contacts_total)$stats
  tot_inf_by_by_contacts <- as_tibble(t(tot_inf_by_by_contacts))
  colnames(tot_inf_by_by_contacts) <- c("ymin", "Q1", "median", "Q3", "ymax")
  tot_inf_by_by_contacts <- tot_inf_by_by_contacts%>% 
    mutate(metric = "tot_inf_by_by_contacts")
  
  tot_contacts <- boxplot.stats(contacts_total)$stats
  tot_contacts <- as_tibble(t(tot_contacts))
  colnames(tot_contacts) <- c("ymin", "Q1", "median", "Q3", "ymax")
  tot_contacts <- tot_contacts%>% 
    mutate(metric = "total_contacts")
  
  out <- bind_rows(out, 
                   tot_inf_by_by_contacts,
                   tot_contacts)
  
  # Loss function 3: total infections prevent (i.e., I no inverention - I intevention)/ total contact prevented  (i.e., C no inverention - C intevention)
  incidence_time_stats_baseline <- do.call("cbind", lapply(baseline, function(x){x$incidence}))
  incidence_total_baseline <- apply(incidence_time_stats_baseline, 2, FUN = sum)

  contacts_time_stats_baseline <- do.call("cbind", lapply(baseline, function(x){x$average_interaction*692})) # average_interaction[k] <- sum(n_interactions)/length(actors) 
  contacts_total_baseline <- apply(contacts_time_stats_baseline, 2, FUN = sum)
  
  incidence_diff <- (incidence_total_baseline)-(incidence_total)
  contact_diff <- (contacts_total_baseline)-(contacts_total) + 1e-6
  efficiency_per_contact <- incidence_diff[incidence_diff>=0]/contact_diff[incidence_diff>=0]
  
  efficiency_per_contact <- boxplot.stats(efficiency_per_contact)$stats
  efficiency_per_contact <- as_tibble(t(efficiency_per_contact))
  colnames(efficiency_per_contact) <- c("ymin", "Q1", "median", "Q3", "ymax")
  efficiency_per_contact <- efficiency_per_contact%>% 
    mutate(metric = "efficiency_per_contact")
  
  mean_inc_baseline <- mean(incidence_total_baseline)
  mean_inc_intervention <- mean(incidence_total)
  mean_contacts_baseline <- mean(contacts_total_baseline)
  mean_contacts_intervention <- mean(contacts_total)
  mean_efficiency <- (mean_inc_baseline - mean_inc_intervention) / 
    (mean_contacts_baseline - mean_contacts_intervention + 1e-6)
  
  mean_efficiency_per_contact <- tibble(mean = mean_efficiency) %>% 
    mutate(metric = "mean_efficiency_per_contact")
  
  out <- bind_rows(out, 
                   efficiency_per_contact,
                   mean_efficiency_per_contact)
  
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
      mutate(res = map(type, ~extract_function(run_ABM = data_list[[.]],
                                               baseline = data_list[["run_ABM_no_intervention"]]))) %>% 
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

