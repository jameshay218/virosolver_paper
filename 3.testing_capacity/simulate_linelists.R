########################################
## SCRIPT TO SIMULATE LINE LIST DATA
## 2nd November 2020
## ------ Aims ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
#' 4. Subset the simulated line list based on a specified observation model, accounting for eg. changing testing
#' 5. Simulate Ct values for the sub-setted line list, giving a Ct value for each observed individual
########################################

########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lazymcmc)
library(extraDistr)
devtools::load_all("~/Documents/GitHub/virosolver")

## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/"
setwd(paste0(HOME_WD,"/virosolver_paper/"))

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

set.seed(1)

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
model_pars <- read.csv("pars/massachusetts/partab_seir_model.csv")
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulation parameters
population_n <- 1000000
times <- 0:200
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 

########################################
## 3. Full simulated line list
########################################
## Simulate SEIR dynamics, incidence, growth rates, Rt etc
## Use "ode" for deterministic
## Use "odin" for stochastic
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, ver="ode")
data_wd <- paste0(HOME_WD, "/virosolver_paper/data/ReportSims/Symptom")
save(seir_dynamics,file=paste0(data_wd,"/SEIR_dynamics.Rda"))

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                                  population_n=population_n)

########################################
## 4. Simulate observation process
########################################
## Choose which process to use
## You could combine these as well, just make sure that each individual is only sampled once.
## ie. use frac_report_increase on observed_individuals; subset observed_individuals by individuals 
## who weren't sampled with frac_report_increase; use prob_report_increase on the remaining individuals

## Simulate time-varying reporting **fraction**
## This represents a fraction of the population being tested each day
frac_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.00001,end_prob=0.001,
                                                                   growth_rate=0.1,switch_point=120), ver="increase")
frac_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.0001,end_prob=0.00001,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Simulate time-varying reporting **probability**
## This gives the probability of an individual getting reported if they have symptoms
prob_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.01,end_prob=0.1,
                                                                   growth_rate=0.05,switch_point=250), ver="increase")
prob_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.1,end_prob=0.01,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Choose which version you want
reporting_prob <- 1
#observed_linelist <- simulate_reporting(complete_linelist, 
#                                        frac_report=reporting_prob,
#                                        timevarying_prob=NULL,
#                                        solve_times=times, 
#                                        symptomatic=FALSE)

## How to use other versions
# ## Do simulation recovery
# observed_indivs_increasing <- simulate_reporting(complete_linelist, timevarying_prob=frac_report_increase, 
#                                                  solve_times=times, symptomatic=FALSE)
# observed_indivs_decreasing <- simulate_reporting(complete_linelist, timevarying_prob=frac_report_decrease, 
#                                                  solve_times=times, symptomatic=FALSE)
# 
# ## Symptomatic surveillance
 observed_indivs_symptom <- simulate_reporting(complete_linelist %>% filter(is_infected == 1), frac_report=reporting_prob, 
                                               solve_times=times, symptomatic=TRUE)
# observed_indivs_symptom_increasing <- simulate_reporting(complete_linelist, timevarying_prob=prob_report_increase, 
#                                                          solve_times=times, symptomatic=TRUE)
# observed_indivs_symptom_decreasing <- simulate_reporting(complete_linelist, timevarying_prob=prob_report_decrease, 
#                                                          solve_times=times, symptomatic=TRUE)


########################################
## 5. Simulate Ct values
########################################
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_symptom$sampled_individuals,
                                                      kinetics_pars=pars)

obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)

obs_dat %>% filter(ct < 40) %>% group_by(t) %>% summarize(ct=mean(ct)) %>% ggplot() + geom_point(aes(x=t,y=ct))

## Note you can now source plot_for_sim_linelist.R to view a figure like in the preprint
#write_csv(obs_dat, "data/sim_linelist.csv")
