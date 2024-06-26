# ssd_paper_code

This is the code for the paper:
"Spatially Selected and Dependent Random Effects for Small Area Estimation with Application to Rent Burden"

Description of core code:

* `download_data.R` use this file to download the data used for the paper. it will be stored in `data/`
* `nc_simulation.R` run this with `n.sims=100` simulations to reproduce the empirical simulation study with 100 simulations
* `fit_models_data_analysis` run this to fit the various models for the data analysis (WARNING: will take a long time)
* `paper_material.Rmd` all code used to create the plots, table, and results referenced in the paper


Description of folders:

* `samplers/` include the MCMC samplers for all models used
* `data/` includes all data used.
* `sim_analysis_functions/` functions used to download data, set up & run simulations, and process results are included here
*  `sim_results_NC` stores results for the simulation study
* `data_analysis` stores MCMC chains for the data analysis
* `figs` png copies of the figures in the paper

Full results are not included:

* Due to file sizes, the repo does not have the full simulation study results saved (only 4 simulations are included). However, the summary file  `sim_results_NC/results_by_sim.RDA` contain results for full study with 100 simulations.  
* Due to file sizes, the repo does not have the full MCMC chains saved (only small chains with 25 iterations are included for testing)
