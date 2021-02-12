########################################
## SCRIPT TO PERFORM SIM-RECOVERY ON NURSING HOME CTS
## James Hay 3 November 2020
##  - Simulation-recovery of nursing homes with different population sizes and prior strengths
##  - First part of the script simulates an outbreak in a nursing home with specified population size and R0
##  - Then, simulates observed Ct values at the specified observation times. Here, use 3 obs times for growth, peak and decline phases
##  - Finally, re-fit the SEIR model to these simulations

########################################
## 1. Headers
########################################
library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
library(coda)

HOME_WD <- "~/Documents/GitHub/"
#HOME_WD <- "~"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

source(paste0(HOME_WD,"/virosolver_paper/code/priors.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/plot_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/linelist_sim_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/odin_funcs.R"))


## Arguments for this run. Note that this will automatically find when Rt is 1 (for the peak) and fit to 20 days either side
control_table <- read_csv(paste0(HOME_WD,"/virosolver_paper/pars/nursing_homes/sim_control_nh.csv"))
## Get task ID, used to read options from control table
simno <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
simno <- 1
## Set random seed
set.seed(simno)

########################################
## 2. Read in run control parameters
########################################
## Name of this run
runname <- control_table$run_name[simno] 
## Which subfolder to save in? Usually different, but should remain unchanged if running multiple chains
run_index <- control_table$run_index[simno] 
## Using all Cts or just positives?
use_pos <- control_table$use_pos[simno] 
## Which incidence model are we assuming?
model_version <- control_table$model_version[simno] 
## Where to store MCMC chains?
top_chainwd <- control_table$chainwd[simno] 
## Where to save output plots?
top_plotwd <- control_table$plotwd[simno] 
## Pointer to prior function to use
prior_func_use <- get(control_table$prior_func[simno])
## Pointer to incidence function to use
inc_func_use <- get(control_table$inc_func[simno])
## Maximum allowable age of infection?
age_max <- control_table$age_max[simno]
## Population size
population_n <- control_table$population_n[simno]
## Prior sd mod
prior_sd_mod <- control_table$prior_mod[simno]


## Most importantly, where is the parameter control table stored?
parTab <- read.csv(control_table$parTab_file[simno], stringsAsFactors=FALSE)
parTab_seir <- read.csv(control_table$seir_parTab_file[simno], stringsAsFactors=FALSE)
n_samp <- 1000

## Manage MCMC runs and parallel runs
nchains <- 3
#n_clusters <- 3
#cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
#registerDoParallel(cl)

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=30000,"popt"=0.234,"opt_freq"=1000,
                 "thin"=50,"adaptive_period"=20000,"save_block"=1000)

## Arrange run name and working directories
full_runname <- paste0(runname,"_",model_version,"_pos",use_pos,"_n", population_n,"_prior",prior_sd_mod,"_run",run_index)
chainwd <- paste0(top_chainwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod,"/",run_index)
plot_wd <- paste0(top_plotwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod,"/",run_index)
if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc_ct <- TRUE

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
pars <- parTab$values
names(pars) <- parTab$names
## Means for priors
means <- parTab$values
names(means) <- parTab$names

## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
if(!is.na(prior_sd_mod)){
  sds_seir <- sds_seir*prior_sd_mod
  sds_exp <- sds_exp*prior_sd_mod
} else {
  prior_func_use <- NULL
}

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
seir_pars <- parTab_seir$values
names(seir_pars) <- parTab_seir$names
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=seir_pars, version="odin",switch_model=FALSE)

## Estimate growth rates on day of peak incidence and 3 weeks either side
peak_time <- times[which.max(seir_dynamics$incidence)] #seir_dynamics$seir_outputs$step[which.min(abs(1-seir_dynamics$seir_outputs$Rt))]
growth_time <- peak_time - 21
decline_time <- peak_time + 21
obs_times <- c(growth_time, peak_time, decline_time)
run_table <- tibble(label=c("growth","peak","decline"), obs_time=obs_times)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                               population_n=population_n)
## Save simulated line list
write_csv(complete_linelist, path=paste0(plot_wd,"/full_linelist.csv"))


observed_linelist_growth <- simulate_reporting(complete_linelist, 
                                               frac_report=1,
                                               timevarying_prob=NULL,
                                               solve_times=times, 
                                               symptomatic=FALSE)

## For growth, peak and decline phases
if(rerun_mcmc_ct){
  for(i in 1:nrow(run_table)){
  #res <- foreach(i=1:nrow(run_table),.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
    #devtools::load_all(paste0(HOME_WD,"/virosolver"))
    
    obs_time <- run_table$obs_time[i]
    lab1 <- run_table$label[i]
    
    if(!file.exists(paste0(chainwd,"/",lab1))) dir.create(paste0(chainwd,"/",lab1),recursive = TRUE)
    
    ## Growth phase observation
    observed_linelist_growth$sampled_individuals$sampled_time <- obs_time
    simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist_growth$sampled_individuals,
                                                          kinetics_pars=seir_pars)
    
    obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
      rename(t = sampled_time, ct=ct_obs) %>% arrange(t)
    
    write_csv(obs_dat, path=paste0(plot_wd,"/",lab1,"_",obs_time,"_cts.csv"))
    
    ## Observation times
    if(!is.na(age_max)){
      obs_dat <- obs_dat %>% mutate(t = t - max(t),t = t + age_max)
    }
    ## Only from start to observation time
    ages <- 1:max(obs_dat$t)
    times <- 0:max(obs_dat$t)
    
    ## Epidemic cannot start after first observation time
    parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat$t)
    
    ## Run for each chain
    chains <- NULL
    for(j in 1:nchains){
      ## Get random starting values
      startTab <- generate_viable_start_pars(parTab,obs_dat,
                                             create_posterior_func,
                                             inc_func_use,
                                             prior_func_use)
      covMat <- diag(nrow(startTab))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
      
      output <- run_MCMC(parTab=startTab,
                         data=obs_dat,
                         INCIDENCE_FUNC=inc_func_use,
                         PRIOR_FUNC = prior_func_use,
                         solve_likelihood=TRUE,
                         mcmcPars=mcmcPars_ct,
                         filename=paste0(chainwd,"/",lab1,"/",lab1,"_",obs_time,"_",run_index,"_chainno_",j),
                         CREATE_POSTERIOR_FUNC=create_posterior_func,
                         mvrPars=mvrPars,
                         OPT_TUNING=0.2,
                         use_pos=use_pos,
                         t_dist=NULL)
      
      ## Read in chain and remove burn in period
      chain <- read.csv(output$file)
      chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
      chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
      chains[[j]] <- chain
    }
    chain <- do.call("bind_rows",chains)
  }
}

res <- NULL
reruns <- NULL
min_ess_list <- NULL
gr_results <- list()
gr_samples <- list()
## For growth, peak and decline phases
for(i in 1:nrow(run_table)){
  obs_time <- run_table$obs_time[i]
  lab1 <- run_table$label[i]
  
  ## Fits
  obs_dat <- read_csv(paste0(plot_wd,"/",lab1,"_",obs_time,"_cts.csv"))
  
  ## Observation times
  if(!is.na(age_max)){
    obs_dat <- obs_dat %>% mutate(t = t - max(t),t = t + age_max)
  }
  ## Only from start to observation time
  ages <- 1:max(obs_dat$t)
  times <- 0:max(obs_dat$t)
  
  chainwd_tmp <- paste0(chainwd,"/",lab1,"/")
  unfixed_list <- load_mcmc_chains(chainwd_tmp, parTab,TRUE,1,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=FALSE)
  
  min_ess <- min(effectiveSize(unfixed_list$list))
  gelmans <- gelman_diagnostics(unfixed_list$list)
  
  chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=TRUE)$chain
  chain <- as.data.frame(chain)
  res[[i]] <- chain
  reruns[[i]] <- gelmans$Rerun
  min_ess_list[[i]] <- min_ess
  
  ## Pull pre-read chain and get this observation time
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  ## Get smoothed growth rates
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    #trajs[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
    trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  }
  
  predictions <- plot_prob_infection(chain_comb, 100, inc_func_use,
                                     times,
                                     obs_dat=obs_dat)
  model_func <- create_posterior_func(parTab,obs_dat,NULL,inc_func_use,"model")
  p2 <- plot_distribution_fits(chain_comb, obs_dat, model_func,100)
  
  ## Growth rates
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  prob_growth <- sum(trajs1[,ncol(trajs1)] > 0)/n_samp
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- times[2:length(times)]
  colnames(trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate samples
  gr_samples[[i]] <- tibble(daily_gr=trajs1[,ncol(trajs1)],runname=runname,label=lab1, index=run_index,obs_time=obs_time)
  
  best_traj <- pmax(inc_func_use(get_best_pars(chain_comb),times),0.0000001)
  best_grs <- log(best_traj[2:length(best_traj)]/best_traj[1:(length(best_traj)-1)])
  best_gr <- best_grs[length(best_grs)]
  
  res_table <-  trajs1_quants %>% filter(t == max(t)) %>%
    mutate(MAP=best_gr) %>%
    mutate(label=lab1) %>%
    mutate(index=run_index) %>%
    mutate(obs_time=obs_time) %>%
    mutate(runname=runname) %>%
    mutate(prob_growth = prob_growth,
           min_ess=min_ess,
           rerun=gelmans$Rerun)
  gr_results[[i]] <- res_table
      
  ## Growth rate plot
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + 
    coord_cartesian(ylim=c(-0.5,0.5))
}
gr_samples_comb <- do.call("bind_rows", gr_samples)
gr_results <- do.call("bind_rows", gr_results)
write_csv(gr_samples_comb, path=paste0(plot_wd,"/gr_samples.csv"))
write_csv(gr_results, path=paste0(plot_wd,"/diagnostics.csv"))
