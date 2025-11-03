# ABM
 
This repository contains code for the ABM implemented in Chapter 4 of the thesis

The repository is structured as follows:

# `scripts/`
 - ## `prepare_data/`
    - `setup_network.R`: Script to construct daily contact networks for the ABM, where the last 2 weeks of data are repeated 8 times. Additionally, the script generates the aggregated network over the first two weeks for the LSHM approach
 - ## `run_clustering/`
    - `run_JANE.R`: Script to fit the LSHM on the contact network over the first two weeks of data
    - `fb_network.R`: Script to fit the LSHM, without noise, on the Facebook friendship network
    - `call_network.R`: Script exploring fitting the LSHM on the phone call network
 - ## `ABM.R`: Main script to run the ABM
 - ## `helper_functions.R`: Script with helper functions used in `ABM.R`
 - ## `modularity_sims.R`: Script to run the ABM for the varying percentage of within-cluster edges analyses
 - ## `extract_results.R`: Script to extract results from ABM runs
 - ## `plots_paper.R`: Script to generate plots presented in Chapter 4 of the thesis
 - ## `plotting_script.R`: Plotting script
