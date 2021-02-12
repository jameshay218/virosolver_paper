########################################
## SCRIPT TO SIMULATE LINE LIST DATA FOR TESTING ANALYSES
## 7th December 2020
## ------ Aims ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
########################################

########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lazymcmc)
library(extraDistr)
library(ggthemes)

HOME_WD <- "~/Documents/GitHub/"
#HOME_WD <- "~"

set.seed(1000)

devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Load functions for line list simulation
source(paste0(HOME_WD,"/virosolver_paper/code/priors.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/plot_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/linelist_sim_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/odin_funcs.R"))

## Some plot settings
export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                     axis.text.y=element_text(size=7),
                                     axis.title.x=element_text(size=8),
                                     axis.title.y=element_text(size=8))

## Creating and Setting Directories:
main_wd <- paste0(HOME_WD,"/virosolver_paper/")
chainwd <- paste0(HOME_WD, "/virosolver_paper/mcmc_chains/ReportSims/Symptom")
plot_wd <- paste0(HOME_WD, "/virosolver_paper/plots/ReportSims/Symptom")
data_wd <- paste0(HOME_WD, "/virosolver_paper/data/ReportSims/Symptom")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(data_wd)) dir.create(data_wd,recursive = TRUE)

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
model_pars <- read.csv(paste0(main_wd,"pars/massachusetts/partab_seir_model.csv"))
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

save(seir_dynamics,file=paste0(data_wd,"/SEIR_dynamics.Rda"))

for(Sim in 1:100){
  print(Sim)
  complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,symp_frac=0.35,population_n=population_n)
  write_csv(x=complete_linelist, path=paste0(data_wd,"/CompleteLineList_",Sim,".csv"))
  
  ## Create summary file that has infection days, true growth rates, and true Rts for each simulation:
  infection_days <- complete_linelist %>% group_by(infection_time) %>% 
    summarize(N=n()) 
  infection_days$GR <- c(NA,log(infection_days$N[2:length(infection_days$N)]/infection_days$N[1:(length(infection_days$N)-1)]))
  infection_days$logN <- log(infection_days$N)
  infection_days$GR14 <- NA
  infection_days$GR28 <- NA
  infection_days$GR35 <- NA
  d14 <- 1:14
  d28 <- 1:28
  d35 <- 1:35
  for (i in 1:length(infection_days$GR)) {
    if (i-13 > 1) {
      logNs <- infection_days$logN[(i-13):i]
      res <- lm(logNs~d14)
      infection_days$GR14[i] <- unname(coef(res)["d14"])
      if (i-27 > 1) {
        logNs <- infection_days$logN[(i-27):i]
        res <- lm(logNs~d28)
        infection_days$GR28[i] <- unname(coef(res)["d28"])
        if (i-34 > 1) {
          logNs <- infection_days$logN[(i-34):i]
          res <- lm(logNs~d35)
          infection_days$GR35[i] <- unname(coef(res)["d35"])
        }
      }
    }
  }
  infection_days <- infection_days %>% left_join(seir_dynamics$seir_outputs %>% dplyr::select(step, Rt) %>% rename(infection_time=step))
  write_csv(infection_days, path=paste0(data_wd,"/TrueGRs_",Sim,".csv"))
}
