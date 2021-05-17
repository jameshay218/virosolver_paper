# virosolver paper additional code and scripts
------------
Scripts to accompany virosolver package and paper. This is a fairly long README structured as follows:
1. Setup: required libraries and external code.
2. General tips: general principles for the repo structure and initial scripts to run.
3. Analysis of Massachusetts nursing homes.
4. Simulated nursing home analyses.
5. Simulation study of changing testing capacities.
6. Simulated Massachusetts analyses.
7. Analysis of BWH hospital data.
 
## 1. Setup
------------
There are 3 central R-packages to these analyses that are not on CRAN:
1. The `lazymcmc` R package, which is used for the MCMC procedure. This is easy to do with `devtools::install_github("jameshay218/lazymcmc")`. *However*, for many of the analyses, a separate branch implementing parallel tempering is needed. I'd recommend you set this up as follows:
  - Install the `lazymcmc` base package using `devtools::install_github` as above. Any time this version is used, `library(lazymcmc)` is called.
  - Clone the `parallel_tempering` branch from [here](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering). Whenever this version is needed, then `devtools::load_all("PATH TO LAZYMCMC PARALLEL TEMPERING REPO")` is called instead.
2. The `virosolver` R-package, containing all of the Ct models and associated functions used here. This is the R-package you would use for your own datasets. As above, simply install with `devtools::install_github("jameshay218/virosolver")`.
3. The `EpiNow2` R-package. This is available with instructions at [here](https://github.com/epiforecasts/EpiNow2), maintained and developed by colleagues at LSHTM. This package is simply used to generate R(t) estimates for the various comparisons throughout the paper. Note that this uses [RStan](https://mc-stan.org/users/interfaces/rstan).

A number of generic R packages are also used throughout:
```r
c("tidyverse","ggthemes","ggpubr","ggsci","data.table","patchwork",
"fitdistrplus","deSolve","lazymcmc","odin","doParallel","coda")
```
 
## 2. General tips
Clone this git repo to your PC and change the `HOME_WD` variable wherever you see it at the top of the scripts.

### Viral kinetics parameters
Assumptions about viral/Ct kinetics are based on our data, so may not be the same for yours. For example, we assume a limit of detection at a Ct of 40 (the "intercept" parameter). You can change these parameters in the appropriate `partab` file in the `pars` directory. 

The `partab` structure is used throughout the analyses and by the MCMC framework. This specifies the parameter bounds, whether the parameter is estimated or not (if `fixed=0`, then we estimate it. If `fixed=1`, then it is, unsurprisingly, kept fixed throughout).

### Outputs
The code saves the MCMC chains to the `mcmc_chains` folder. Some plots are also automatically generated in the `plots` folder. Otherwise, some output plots are just generated as R objects in this R session.

*NOTE:* due to the file sizes of the generated MCMC chains and diagnostic plots, these are not included in this repository. Please either regenerate the relevant analyses or contact `jameshay218` here.

### MCMC
It is really important to check the chain convergence. Look for the lines that call the function `gelman_diagnostics`. If chains are not converged, then they can be run for longer by changing the `mcmcPars` variables in the fitting scripts. Increase the `interations` and `adaptive_period` values.

### Pre-scripts
These scripts are used to generate inputs and parameter values for the main analyses.
* a.fit_viral_load_model.R: run this script to find the viral kinetics parameters that fit the Borremans et al. data set. Part of the optimization is fitting to the proportion detectable curve by least squares, but optionally the optimization also accounts for the lower 1% quantile of Ct values at peak viral load and the upper 1% at day 30 post infection. This can be done separately for the nursing home data and BWH data, which show slightly different distributions. This saves the resulting parameter estimates to pars/XX/partab_fitted.csv, where XX is either `nursing_home` or `massachusetts`. Note that a sensitivity analysis where the Borremans et al. data are shifted down in the y-axis systematically can be run by adjusting the proportion detectable data.
* b.fit_rt_epinow2.R: estimates Rt using the EpiNow2 package based on case counts reported by NYT for Massachusetts and the county of Suffolk.

## 1. Analysis of Massachusetts nursing homes
All scripts for these analyses are in the `1.nursing_homes`. From here, run the `1.nursing_home_analyses_combined.R` script which compiles all of the distinct parts. Scripts are then sourced as follows:
  - 1.nursing_home_analyses_headers.R: Loads in the nursing home Ct data, does some cleaning and plotting, and sets up the parameters for the subsequent analyses.
  - 1.nursing_home_analyses_SEEIRR.R: This fits a compartmental SEEIRR model to the point prevalence estimates for each nursing home. The model trajectories are then used to generate daily growth rate estimates. Uses the normal version of `lazymcmc`.
  - 1.nursing_home_analyses_ct_model.R: Takes the cross sections for each nursing home and separately fits exponential and SEIR Ct models to the Ct distributions. Note that for each incidence model version, the model is fitted using all Ct values, but flags can be set to fit to only positives. Also creates the `figures/supplement/nh_posteriors` figure. Uses the normal version of `lazymcmc`.
  - Figure2.R: this is run at the end of `1.nursing_home_analyses_combined.R`. It uses the MCMC chains generated and read in from the above scripts to create Figure 2, saved in the `figures` folder.
  - SuppFigure_fits.R: creates the corresponding Figure 2 plots for nursing homes 2-4. Saved as two plots in the `figures/supplement` folder.

## 2. Simulated nursing home analyses
All scripts for these analyses are in the `2.sim_nursing_homes`. Most of these analyses were run on the Harvard FASRC using a shell script. You'll see that L35 of `2.simulated_nursing_homes.R` reads in a control table specifying model parameters and assumptions, and the SLURM task ID indicates which row of the control table to run. You can set `simno` on L37 to run a single job locally.
  - 2.simulated_nursing_homes.R: Simulates an SEIR model and Ct values under the settings and parameters defined in `pars/nursing_homes/sim_control_nh.csv`. Uses the normal version of `lazymcmc`.
  - 2.combine_sim_NH_results.R: Once all of the jobs in the previous script are run, this script is used to compile all of the different jobs and generates all of the supplementary figures for these analyses.
  - nursing_home_jobs.sh: The shell script used to run `2.simulated_nursing_homes.R` in parallel on the cluster. Generates all of the MCMC chains and plots in the `mcmc_chains` and `plots` directories, which are then read in in the `2.combine_sim_NH_results.R` script.
  
## 3. Simulation study of changing testing capacities.
Analysis of how `EpiNow2` (Rt estimation) and `virosolver` (Ct model) fare at re-estimating the true infection incidence when testing capacity is changing over time (e.g., number of tests carried out per week increases exponentially). First simulates **N** distinct populations from an SEIR model to give the true infections over time, and then simulates different testing strategies for these populations. We investigated both testing of only symptomatic individuals (people are tested with some delay after their incubation period) or random cross-sections of the population (people are randomly sampled from the population and tested). `EpiNow2` and `virosolver` are then run separately on the same simulations for comparison. Results are save to the `results`, `mcmc_chains` and `plots` directories under the `ReportSims` sub-directories.
  - simulate_linelists.R: simulates the base R objects for the subsequent SEIR simulations, saved in the `data` folder. Uses the settings in the script and `pars/massachusetts/partab_seir_model.csv`. Also demonstrates how the accompanying testing capacity simulation functions are run and fed into the viral load simulation function. The main result of this script is saving the `SEIR_dynamics.Rda` object, which forms the basis of the SEIR simulation for all analyses in this section. The `_switch` version of this script uses the bimodal peak SEIR model instead.
  - 3.symptomatic_testing_simulate_linelists.R: Simulates full linelists for each simulation run using the SEIR model with epi parameters specified in `pars/massachusetts/partab_seir_model.csv`. Run locally. These simulations are used for the symptomatic testing simulations.
  - 3.changing_test_rates.R: Fits the exponential growth model to a time-series dataset of randomly sampled Ct values. Testing strategies tested are flat, rising and falling testing schemes. The estimated 35-day growth rate is compared to the growth rate based on the raw number of detected infections.
  - changing_testing_simulations.sh: shell script to run `3.changing_test_rates.R` on the cluster.
  - 3.symptomatic_testing_epinow2.R: Uses the simulations from `3.symptomatic_testing_simulate_linelists.R` and generates R(t) estimates based on reported symptomatic cases.
  - 3.symptomatic_testing_virosolver_new.R: Uses the simulations from `3.symptomatic_testing_simulate_linelists.R` and generates growth rate estimates using the `virosolver` Ct models.  The `_samp_size` version of this script runs the same thing but varies the number of samples taken overall. Additionally, the `_samp_size` script also fits a "prevalence-only" model to each timepoint/simulation, which fits either an SEIR or SEEIRR model to cross-sectional prevalence assuming that the E+I compartment sizes are equivalent to RT-qPCR prevalence.
  - testing_simulations_epinow2.sh: shell script to run `3.symptomatic_testing_epinow2.R`.
  - testing_simulations_virosolver.sh: shell script to run `3.symptomatic_testing_virosolver.R`.

## 4. Simulated Massachusetts analyses.
Simulation of repeated cross-sections arising from an SEIR epidemic. One script is provided to run a single simulation-recovery experiment locally, however, for our analyses we used the included shell scripts to perform multiple simulation-recovery experiments. These simulations match the subsequent analysis of the BWH data, and use the parameter tables in `pars/massachusetts`.
  - 4.full_sim_recovery.R: Simulates an SEIR model in a large population, with parameters specified in the script. This runs a single simulation-recovery. Uses the normal version of `lazymcmc`.
  - simulate_linelists_cluster.R: Simulates data from the SEIR model, to be called by the shell script `simulate_ma_seir_data.sh`. Can also be run directly locally.
  - 4.full_sim_recovery_cluster.R: As above, but called via the shell script `fit_sim_ma_seir_data.sh`. Uses simulation settings defined in `pars/massachusetts/sim_ma_control_use.csv`.
  - 4.full_sim_recovery_freegp.R: Identical to `4.full_sim_recovery_cluster.R`, but relaxes the assumption that the Gaussian process parameters are fixed. This could be done by editing the code to `4.full_sim_recovery_cluster.R` instead.
  - 4.full_sim_recovery_cluster_pos.R: Identical to `4.full_sim_recovery_cluster.R`, but uses only positive Ct values, thereby not using any information on the proportion positive over time.
  - 4.full_sim_recovery_cluster_late.R: Identical to `4.full_sim_recovery_cluster.R`, but assumes that sampling only begins part way through the epidemic.
  - 4.full_sim_recovery_cluster_sampsize.R: Identical to `4.full_sim_recovery_cluster.R`, but also varies the total number of samples observed from the simulation.
  - 4.full_sim_recovery_cluster_sampsize_prior.R: Identical to `4.full_sim_recovery_cluster_sampsize.R`, but does NOT use the likelihood (ie. it samples directly from the prior).
  - simulate_ma_seir_data.sh: calls the `simulate_linelists_cluster.R` script when run on the cluster.
  - fit_sim_ma_seir_data.sh: calls the `4.full_sim_recovery_cluster.R` script on the cluster.
  - fit_sim_ma_seir_freegp.sh: calls the `4.full_sim_recovery_freegp.R` script on the cluster.
  - fit_sim_ma_seir_late.sh: calls the `4.full_sim_recovery_cluster_late.R` script on the cluster.
  - fit_sim_ma_seir_pos.sh: calls the `4.full_sim_recovery_cluster_pos.R` script on the cluster.
  - fit_sim_ma_seir_sampsize.sh: calls the `.full_sim_recovery_cluster_sampsize.R` script on the cluster.
  - fit_sim_ma_seir_sampsize_prior.sh: calls the `4.full_sim_recovery_cluster_sampsize_prior.R` script on the cluster.
  - plot_gp_sim.R: specify the output location of one simulation-recovery from `4.full_sim_recovery_cluster.R` and generates Movie S1, the gif showing the Gaussian process model fit over repeated cross sections.
   - plot_gp_sim_samplesize.R: specify the output location of one simulation-recovery from `4.full_sim_recovery_cluster_sampsize.R` and generates Movie S2, the gif showing the Gaussian process model fit to simulations of varying sample sizes.
   - plot_gp_sim_late.R: specify the output location of one simulation-recovery from `4.full_sim_recovery_cluster_late.R` and generates Movie S3, the gif showing the Gaussian process model fit over repeated cross sections with sampling beginning partway through.
  - plot_for_sim_linelist.R: plots what was previously Figure 2 in the first preprint version. This just shows lots of plots from a single simulation when run after `4.full_sim_recovery.R`.
  
## 5. Analysis of BWH hospital data.
Fit the various Ct models (single cross sectional SEIR and exponential growth models, and the full Gaussian process model) to the detectable Panther Ct values from BWH. The scripts can be run locally, but some shell scripts are also included to run longer jobs on the cluster.
  - 5.fit_ma_gp.R: The main script for the BWH Gaussian process analyses. Reads in and cleans the Panther Ct data from testing done at BWH. Fits the GP model to the entire timeseries and generates some corresponding plots in the figures folder. Uses the normal version of `lazymcmc`.
  - 5.fit_ma_free_gp.R: Identical to the above script, but allows the Gaussian Process parameters to be estimated within the MCMC fitting.
  - 5.fit_ma_free_reduce_samp.R: Identical to `5.fit_ma_gp.R`, but sub-samples the overall dataset to demonstrate robustness to smaller sample sizes. NOTE that you need to edit the arguments on L31, L32 and L171/172. This script is to be run locally, with arguments changed for each assumed sampling strategy.
  - fit_real_ma.sh: shell script to call `5.fit_ma_gp.R`.
  - fit_real_free_gp.sh: shell script to call `5.fit_ma_free_gp.R`.
  - 5.fit_ma_single_timepoints_exp.R: fits the exponential growth Ct model to each week of data from BWH. Note, uses parallel tempering branch of `lazymcmc`.
  - 5.fit_ma_single_timepoints_seir.R: fits the SEIR Ct model to each week of data from BWH. Note, uses parallel tempering branch of `lazymcmc`.
  - 5.add_prior.R: script called as part of `Figure4.R`.
  - Figure4.R: uses the outputs from the above scripts in this section to create Figure 4. Note that this script uses data from the NY times [covid-19 data github page](https://github.com/nytimes/covid-19-data). Also plots all associated supplementary figures (Fig S14, S15 and S17).