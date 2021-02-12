########################################
## SCRIPT TO SIMULATE AND ANALYZE LINE LIST DATA FOR SYMPTOM-BASED TESTING SCHEMES
## 12th November 2020
## Based on simulate_linelists.R and fit_simulated_linelists.R
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
library(extraDistr)
library(ggthemes)
library(ggpubr)
library(data.table)
library(fitdistrplus)
library(deSolve)
library(EpiNow2) # install.packages("EpiNow2")
# library(lazymcmc) devtools::install_github("jameshay218/lazymcmc) ## Note that in this script, the parallel tempering version is used. See README

HOME_WD <- "~"
#HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))
devtools::load_all(paste0(HOME_WD,"/lazymcmc"))

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
results_wd <- paste0(HOME_WD, "/virosolver_paper/results/ReportSims/Symptom/EpiNow2")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(data_wd)) dir.create(data_wd,recursive = TRUE)
if(!file.exists(results_wd)) dir.create(results_wd,recursive = TRUE)

## Where to perform the simulations
setwd(main_wd)

########################################
## 2. Simulation settings
########################################
## Arguments for this run. Can generat automatically here
control_table <- expand.grid(pop_no=1:100,strategy=LETTERS[1:5],testing_day=as.character(c(7,8,9))) %>% arrange(pop_no, testing_day, strategy)
control_table$run_name <- 1:nrow(control_table)

## Set Simulation Number and get sim settings
Sim <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(paste0("Starting Simulation Number: ",Sim))

## Name of this run
runname <- control_table$run_name[Sim] 
## Which simulated population to use?
pop_no <- control_table$pop_no[Sim] 
## Which testing day to use?
Days <- as.character(control_table$testing_day[Sim])
## Which testing strategy to use?
Probs <- control_table$strategy[Sim]

########################################
## 3. Model parameters and simulation settings
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
## 4. Read in underlying SEIR dynamics and simulated linelist data
########################################
load(file=paste0(data_wd,"/SEIR_dynamics.Rda"))

True_SEIR <- tibble(t=rep(times,3), variable=rep(c("infections","growth_rate","R"), each=length(times)),
                    median=c(seir_dynamics$seir_outputs$I, c(NA,seir_dynamics$growth_rates[seir_dynamics$growth_rates[,"ver"]=="daily","GR"]), 
                             seir_dynamics$seir_outputs$Rt),
                    lower_95=NA, lower_50=NA, upper_50=NA, upper_95=NA)

complete_linelist <- read.csv(paste0(data_wd,"/CompleteLineList_",pop_no,".csv"))
infection_days <- read.csv(paste0(data_wd,"/TrueGRs_",pop_no,".csv"))

########################################
## 5. Specify different observation processes
########################################
## Symptom-Based Testing Schemes:
days <- 36

t_7 <- 1:(which.max(seir_dynamics$incidence)-2*7) # Pre peak
t_8 <- 1:(which.max(seir_dynamics$incidence)+2*7) # Across peak
t_9 <- 1:(which.max(seir_dynamics$incidence)+6*7) # Decline phase

lastdays <- c(max(t_7),max(t_8),max(t_9))
names(lastdays) <- c("t_7","t_8","t_9")

p_A <- rep(.1, length.out=days)
p_B <- seq(from=.1,to=.2, length.out=days)
p_C <- exp(seq(from=log(.1),to=log(.25),length.out=days))
p_D <- seq(from=.1,to=.01, length.out=days)
p_E <- 0.2-exp(seq(from=log(.1),to=log(.195),length.out=days))

p_R <- c(0.003,0.003,0.003)
days_R1 <- c(0) # On last day
days_R2 <- c(7,0) # On week before
days_R3 <- c(14,7,0) # On last day, week before and 2 weeks before

########################################
## 6. Analyze symptom-based case report data
########################################
## Distributions for generation time, incubation period, and reporting delays, using defaults from EpiNow2:
# Incubation Period: uses mean 1.621 and SD 0.418 for lognormal distribution. Cf simulate_observations_wrapper in virosolver::simulation_functions.R
# This is what we simulated from
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

# DON'T USE: INCORRECT FOR SEIR SIM
# Generation time: uses mean 3.635 and SD 3.075 for gamma distribution based on Ganyani et al.
# generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")

## Use generation time from SEIR model. Convolution of two exponential distributions has mean = to sum of means from both and
## variance is sum of variances of both
SI_mean <- pars["infectious"] + pars["incubation"]
SI_sd <- sqrt(pars["infectious"]^2 + pars["incubation"]^2)
generation_time <- list(mean=SI_mean,mean_sd=3,sd=SI_sd,sd_sd=3,max=100)

# Reporting Delays
# Uses shape 5 and rate 2in discrete gamma distribution. Cf simulate_observations_wrapper in virosolver::simulation_functions.R
# This is what we simulated from
reporting_delay <- bootstrapped_dist_fit(extraDistr::rdgamma(1000, 5, 2), max_value = 15,bootstraps = 1)


## Within each simulation, simulate case reporting for each testing scenario:
st <- proc.time()

## Days to use in inference
testdays <- get(paste0("t_",Days))
lastday <- max(testdays)

## For each test day and testing scenario combination, simulate case reporting:
print(paste0("Generating Reporting for and Analyzing Scenario: ",Days,Probs))
probscheme <- get(paste0("p_",Probs))
tvp <- tibble(t=testdays,prob=c(rep(.1,length(testdays)-length(probscheme)),probscheme),ver=Probs)
LL <- simulate_reporting(complete_linelist %>% filter(onset_time <= max(testdays)),
                         timevarying_prob=tvp,
                         solve_times=testdays,
                         symptomatic=TRUE)

## Record number of confirmed cases per day:
Cases <- LL$sampled_individuals %>% dplyr::select(confirmed_time) %>%
  group_by(confirmed_time) %>% summarize(n=n()) %>% rename(confirm=n)

## Avoid right censoring problem for this, just a nuisance
Cases <- Cases %>% filter(confirmed_time <= lastday)

## Fill cases vector with 0 counts
Cases_zero <- tibble(confirmed_time=setdiff(seq(min(testdays),max(Cases$confirmed_time),by=1),unique(Cases$confirmed_time)),
                     confirm=0)
Cases <- bind_rows(Cases, Cases_zero) %>% arrange(confirmed_time)

StartDate <- as.Date("03/15/2020", "%m/%d/%Y")
Cases$date <- StartDate + Cases$confirmed_time

ggplot(Cases) + geom_line(aes(x=date,y=confirm))

## Get estimates from reported case counts:
report_ests <- estimate_infections(reported_cases=Cases,
                                   generation_time=generation_time,
                                   delays=delay_opts(incubation_period, reporting_delay),
                                   rt=rt_opts(prior=list(mean=2, sd=2)),
                                   CrIs = c(0.5, 0.95),
                                   stan=stan_opts(cores=4,control=list(adapt_delta=0.95)),
                                   gp = gp_opts(ls_min = 10, basis_prop = 0.2),
                                   verbose=TRUE)

## Have a look at fitted/predicted trajectories
plot_R <- report_ests$summarised %>% 
  filter(variable == "R") %>% 
  ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower_95,ymax=upper_95),alpha=0.1) + 
  geom_ribbon(aes(x=date,ymin=lower_50,ymax=upper_50),alpha=0.3) + 
  geom_line(data=infection_days %>% filter(infection_time < lastday),
            aes(x=as.Date(infection_time,origin=StartDate),y=Rt))

report_ests$summarised %>% 
  filter(variable == "infections") %>% 
  ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower_95,ymax=upper_95),alpha=0.25) + 
  geom_line(data=infection_days %>% filter(infection_time < lastday),
            aes(x=as.Date(infection_time,origin=StartDate),y=N*0.1*0.35))

report_ests$summarised %>% 
  filter(variable == "reported_cases") %>% 
  ggplot() + 
  geom_ribbon(aes(x=date,ymin=lower_95,ymax=upper_95),alpha=0.25) + 
  geom_line(data=Cases,aes(x=date,y=confirm))

## Save EpiNow2 estimates
Ests_output <- report_ests$summarised %>% dplyr::select(!strat) %>%
  dplyr::filter(variable %in% c("R","infections","prior_infections","growth_rate"))
Ests_output$t <- as.numeric(Ests_output$date-StartDate)

## Compile EpiNow2 estimates, simulation case counts and true SEIR dynamics and save to disk
Ests_Full <- bind_cols(Sim=Sim, TestDay=Days, TestProbs=Probs, Ests_output %>% dplyr::select(!date))
Cases_Full <- bind_cols(Sim=Sim, TestDay=Days, TestProbs=Probs, Cases)
True_SEIR_sim <- bind_cols(True_SEIR %>% filter(t==lastday),Sim=Sim, TestDay=Days, TestProbs="T")

save(list=c("Cases_Full","Ests_Full","True_SEIR_sim"), file=paste0(results_wd,"/Sim_EN2_pop",pop_no,"_",Days,"", Probs,".Rda"))
ggsave(paste0(plot_wd,"/Plot_Sim",pop_no,"_",Days,"", Probs,".png"), plot_R, width=4, height=3,dpi=90)
