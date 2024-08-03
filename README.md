# ssd_paper_code

This is the code that can reproduce the results in the paper:
"Spatially Selected and Dependent Random Effects for Small Area Estimation with Application to Rent Burden".

Note: the MCMC sampler for the SSD model have since been improved, making it *significantly* faster. This was tested using samplers in `test_speed_new/` and the test_speed_new.Rmd.


*Partial results only:* Due to file sizes, full results are not stored here.

* The repo does not have the full simulation study results, only 4 simulations are included. However, the summary file  `sim_results_NC/results_by_sim.RDA` contain results for full study with 100 simulations.  
* The repo does not have the full MCMC chains saved. Only small chains with 25 iterations are included for testing.

***

**Code**

* `download_data.R` use this file to download the data used for the paper. it will be stored in `data/`
* `nc_simulation.R` run this with `n.sims=100` simulations to reproduce the empirical simulation study with 100 simulations
* `fit_models_data_analysis` run this to fit the various models for the data analysis (WARNING: will take a long time)
* `paper_material.Rmd` all code used to create the plots, table, and results referenced in the paper

**Folders**

* `samplers/` include the MCMC samplers for all models used
* `data/` includes all data used.
* `sim_analysis_functions/` functions used to download data, set up & run simulations, and process results are included here
*  `sim_results_NC/` stores results for the simulation study
* `data_analysis/` stores MCMC chains for the data analysis
* `figs/` png copies of the figures in the paper
