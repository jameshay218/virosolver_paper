########################################
## SCRIPT TO SIMULATE LINE LIST DATA ON CLUSTER
## 6th November 2020
## ------ Aims ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
#' 4. Subset the simulated line list based on a specified observation model, accounting for eg. changing testing
#' 5. Simulate Ct values for the sub-setted line list, giving a Ct value for each observed individual
#' NOTE: can run this using the cluster script bash_scripts/simulate_ma_seir_data.sh with the for loop removed, or can just run as one task
########################################
########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(extraDistr)
library(lazymcmc)
library(patchwork)
library(ggthemes)
library(odin)
HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Where to perform the simulations
setwd(paste0(HOME_WD,"/virosolver_paper"))

## Seed based on SLURM ID
#task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
for(task_id in 1:10){
  set.seed(task_id)
  
  ## Where to save the simulations
  run_name <- paste0("sim_MA_gp_",task_id)
  save_wd <- paste0(HOME_WD, "/virosolver_paper/data/MA_SEIR_sim/")
  save_file <- paste0(save_wd, run_name)
  if(!file.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
  
  ## Load functions for line list simulation
  source("code/linelist_sim_funcs.R")
  source("code/odin_funcs.R")
  source("code/plot_funcs.R")
  ## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
  source("code/priors.R")
  
  ########################################
  ## 2. Model parameters and simulation settings
  ########################################
  ## Model parameters
  model_pars <- read.csv("pars/massachusetts/partab_seir_switch_model.csv")
  pars <- model_pars$values
  names(pars) <- model_pars$names
  
  ## Simulation parameters
  population_n <- 1000000
  ## Over the course of 200 days
  times <- 0:200
  
  ## Sampling
  sampling_frequency <- 14
  sampling_number <- 2000
  
  ########################################
  ## 3. Full simulated line list
  ########################################
  ## Simulate SEIR dynamics, incidence, growth rates, Rt etc
  ## Use "ode" for deterministic
  ## Use "odin" for stochastic
  seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                         pars=pars, ver="odin",switch_model=TRUE)
  
  ## Simulate onset times, confirmation delays etc
  ## This returns a tibble with line list entries for **every** individual in the population
  complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                                    population_n=population_n)
  
  ########################################
  ## 4. Simulate observation process
  ########################################
  sample_probs <- c(rep(0, sampling_frequency-1),sampling_number/population_n)
  sample_probs <- rep(sample_probs, length(times)/sampling_frequency +1)
  sample_probs <- sample_probs[1:length(times)]
  frac_report <- tibble(t=times,prob=sample_probs)
  
  observed_linelist <- simulate_reporting(complete_linelist, 
                                          frac_report=NULL,
                                          timevarying_prob=frac_report,
                                          solve_times=times, 
                                          symptomatic=FALSE)
  
  simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,
                                                        kinetics_pars=pars)
  
  obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
    rename(t = sampled_time, ct=ct_obs) %>% arrange(t)
  
  ## Save simulated line list, Cts and SEIR dynamics
  write_csv(seir_dynamics$seir_outputs, path=paste0(save_file,"_seir_outputs.csv"))
  write_csv(complete_linelist, path=paste0(save_file,"_full_linelist.csv"))
  write_csv(obs_dat, path=paste0(save_file,"_cts.csv"))
  
  ## Save SEIR plots
  ## Ct distribution plot
  p_dat <- ggplot(obs_dat %>% filter(ct < pars["intercept"])) + 
    geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
    geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
    scale_y_continuous(trans="reverse") +
    export_theme +
    scale_x_continuous(limits=c(min(times),max(times)+50)) +
    ylab("Ct value") +
    xlab("Observation time")
  
  ggsave(paste0(save_file,"_ct_dist.png"),p_dat,width=7,height=4,dpi=150)
  ggsave(paste0(save_file,"_seir_model.png"),seir_dynamics$plot,width=7,height=7,dpi=150)
  ggsave(paste0(save_file,"_growth_rates.png"),seir_dynamics$growth_rate_p,width=7,height=6,dpi=150)
  
  
  print("Job completed successfully")
}
