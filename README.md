# ssd_paper_code

This is the code that can reproduce the results in the paper:
"Spatially Selected and Dependent Random Effects for Small Area Estimation with Application to Rent Burden".


**Main Manuscript**

Code: 
* `paper_material.Rmd` all code used to create the plots, table, and results referenced in the paper
* `download_data.R` use this file to download the data used for the paper. it will be stored in `data/`
* `nc_simulation.R` run this with `n.sims=300` simulations to reproduce the empirical simulation study with 300 simulations
* `fitting_chains.Rmd` run this to fit the various models for the data analysis and verify the computation times 


Folders:
* `samplers/` include the MCMC samplers for all models used
* `data/` includes all data used.
* `sim_analysis_functions/` functions used to download data, set up & run simulations, and process results are included here
*  `sim_results_NC/` stores results for the North Carolina simulation study. NOTE: only partial results (4 simulations) are stored due to file size. 
* `data_analysis/` stores MCMC chains for the data analysis
* `figs/` png copies of the figures in the paper


***

**Supplementary Materials**
 
* `suppl_report.Rmd` a notebook to recreate the Supplemental Report
* `preamble.tex` used for the supplemental report 
* `sa_simulation.R` can be used to run a simulation study with the South Atlantic Census Division. This will be time & memory intensive
* `sim_results_SA/` stores results for the S.Atlantic Simulation study. NOTE: only partial results (1 simulation) are stored due to file size. 

