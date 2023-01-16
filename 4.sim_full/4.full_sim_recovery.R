########################################
## SCRIPT TO PERFORM SIM-RECOVERY OF MASSACHUSSETTS EPIDEMIC
## James Hay 4 November 2020
##  - Simulation-recovery of nursing homes with different population sizes and prior strengths
##  - First part of the script simulates an outbreak in a nursing home with specified population size and R0
##  - Then, simulates observed Ct values at the specified observation times. Here, use 3 obs times for growth, peak and decline phases
##  - Finally, re-fit the SEIR model to these simulations

########################################
## 1. Headers
########################################
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
devtools::load_all("~/Documents/GitHub/virosolver")

## Arguments for this run
index <- 1
set.seed(index)
population_n <- 6900000
n_samp <- 100
runname <- "sim_ma_test"
run_version <- "gp" ##gp, seir or exp##
rerun_mcmc <- TRUE

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/3.sim_ma_ct/",
                  runname,"_",population_n,"/",index)
plot_wd <- paste0("~/Documents/GitHub/virosolver_paper/plots/3.sim_ma_ct/",
                  runname,"_",population_n,"/",index)
setwd(main_wd)

## Manage MCMC runs and parallel runs
nchains <- 1
n_clusters <- 3
#cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
#registerDoParallel(cl)

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=10,"adaptive_period"=20000,"save_block"=1000)

## Code for plotting
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc_ct <- TRUE

########################################
## 2. Model parameters and simulation settings
########################################

inc_func_use <- gaussian_process_model
prior_func_use <- prior_func_hinge_gp

## SEIR model parameters
parTab_seir <- read.csv("pars/massachusetts/partab_seir_switch_model.csv")
seir_pars <- parTab_seir$values
names(seir_pars) <- parTab_seir$names

## GP model parameters for fitting
parTab <- read.csv("pars/massachusetts/partab_gp_model.csv")
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

## Simulation parameters
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
                                       pars=seir_pars, switch_model=TRUE)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(floor(seir_dynamics$incidence),times=times,
                                                               population_n=population_n)
## Save simulated line list
write_csv(complete_linelist, path=paste0(plot_wd,"/",runname,"_",population_n,"_",index,"_full_linelist.csv"))

weekly_prob <- c(rep(0, 13),1000/population_n)
probs <- rep(weekly_prob, length(times_extended)/14 +1)
probs <- probs[1:length(times_extended)]
frac_report <- tibble(t=times_extended,prob=probs)

observed_linelist <- simulate_reporting(complete_linelist, 
                                               frac_report=0.1,
                                               timevarying_prob=NULL,
                                               solve_times=times, 
                                               symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,
                                                      kinetics_pars=pars)

obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)
#write_csv(obs_dat, path=paste0(plot_wd,"/",runname,"_",lab1,"_",obs_time,"_",index,"_cts.csv"))


obs_dat1 <- obs_dat %>% filter(t <= 150 & t >= 30)


p_dat <- ggplot(obs_dat1 %>% filter(ct < 40)) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_continuous(limits=c(0,200))

## GP model parameters for fitting
parTab <- read.csv("pars/massachusetts//partab_gp_model.csv")
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

ages <- 1:max(obs_dat1$t)
times <- 0:max(obs_dat1$t)

## This is for the GP version
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
if(run_version == "gp"){
  parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
}

## Epidemic cannot start after first observation time
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat$t)


f <- create_posterior_func(parTab, obs_dat1, prior_func_use, inc_func_use,solve_ver="likelihood",t_dist=t_dist)
f(pars)

## Run for each chain
chains <- NULL
if(rerun_mcmc){
#for(j in 1:nchains){
  j <- 1
  ## Get random starting values
  startTab <- generate_viable_start_pars(parTab,obs_dat1,
                                         create_posterior_func,
                                         inc_func_use,
                                         prior_func_use)
  covMat <- diag(nrow(startTab))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
  
  output <- run_MCMC(parTab=startTab,
                     data=obs_dat1,
                     INCIDENCE_FUNC=inc_func_use,
                     PRIOR_FUNC = prior_func_use,
                     solve_likelihood=TRUE,
                     mcmcPars=mcmcPars_ct,
                     filename=paste0(chainwd,"/",runname,"_",1,"_",1,"_chainno_",1),
                     CREATE_POSTERIOR_FUNC=create_posterior_func,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=FALSE,
                     t_dist=t_dist)
  
  ## Read in chain and remove burn in period
  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
  chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
  chains[[j]] <- chain
  chain <- do.call("bind_rows",chains)
} else {
  chain <- load_mcmc_chains
}
chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
chain_comb <- chain
chain_comb$sampno <- 1:nrow(chain_comb)

## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
trajs <- matrix(0, nrow=n_samp,ncol=length(times))
for(ii in seq_along(samps)){
  trajs[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
  #trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
}

trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs1_quants <- as.data.frame(trajs1_quants)
trajs1_quants$t <- 1:nrow(trajs1_quants)
colnames(trajs1_quants) <- c("lower","median","upper","t")

## Growth rate plot
p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  coord_cartesian(ylim=c(-0.5,0.5))


trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs_quants <- as.data.frame(trajs_quants)
trajs_quants$t <- 1:nrow(trajs_quants)
colnames(trajs_quants) <- c("lower","median","upper","t")

## Growth rate plot
p_inc <- ggplot(trajs_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
  geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
  export_theme +
  ylab("Per capita incidence") +
  xlab("Days since start") +
  scale_x_continuous(limits=c(0,220)) +
  coord_cartesian(ylim=c(0,0.003))
p_dat/p_inc

