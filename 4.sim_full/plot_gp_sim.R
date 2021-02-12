########################################
## SCRIPT TO PERFORM SIM-RECOVERY OF MASSACHUSSETTS EPIDEMIC ON THE CLUSTER
## James Hay 6 November 2020
##  - Simulation-recovery of nursing homes with different population sizes and prior strengths
##  - First part of the script simulates an outbreak in a nursing home with specified population size and R0
##  - Then, simulates observed Ct values at the specified observation times. Here, use 3 obs times for growth, peak and decline phases
##  - Finally, re-fit the SEIR model to these simulations
##  - CLUSTER NOTE: this reads in pre-simulted Ct data and SEIR dynamics and fits the Ct 
##    model with specified settings to whatever timepoint was chosen

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
library(magick)
library(gganimate)
HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source(paste0(HOME_WD,"/virosolver_paper/code/priors.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/plot_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/linelist_sim_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/odin_funcs.R"))

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=700000,"popt"=0.44,"opt_freq"=5000,
                 "thin"=200,"adaptive_period"=300000,"save_block"=10000)

## Arguments for this run, controlled by a csv file
#control_table <- read_csv(paste0(HOME_WD,"/virosolver_paper/pars/sim_ma_control_use.csv"))
control_table <- readxl::read_excel(paste0(HOME_WD,"/virosolver_paper/pars/massachusetts/sim_ma_control.xlsx"))
control_table$chainwd <- "~/Documents/GitHub/virosolver_paper/mcmc_chains/3.sim_ma_ct/sim_MA_gp/2/"
ps <- NULL
trajs_quants_all <- NULL
obs_dat_use_all <- NULL
## Get task ID, used to read options from control table
for(index in 4:14){
  simno <- index#as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  ## Set random seed
  set.seed(simno)
  
  ########################################
  ## 2. Read in run control parameters
  ########################################
  ## Name of this run
  runname <- control_table$run_name[simno] 
  ## Which subfolder to save in? Usually different, but should remain unchanged if running multiple chains
  run_index <- control_table$run_index[simno] 
  ## File for the SEIR simulation
  data_wd_seir <- control_table$data_wd_seir[simno] 
  ## File for the Ct values
  data_file_cts <- control_table$data_file_cts[simno] 
  ## Using all preceding times or just the most recent time?
  cumu_times <- control_table$cumu_times[simno] 
  ## Using all Cts or just positives?
  use_pos <- control_table$use_pos[simno] 
  ## Is the GP prior fixed or estimated? (if applicable)
  gp_fixed <- control_table$gp_fixed[simno] 
  ## Which incidence model are we assuming?
  model_version <- control_table$model_version[simno] 
  ## Where to store MCMC chains?
  top_chainwd <- control_table$chainwd[simno] 
  ## Where to save output plots?
  top_plotwd <- control_table$plotwd[simno] 
  ## In the vector of timepoints, which observation time are we at?
  timepoint <- control_table$timepoint[simno] 
  ## If running repeats, what's the MCMC chain index?
  chainno <- control_table$chainno[simno]
  ## Pointer to prior function to use
  prior_func_use <- get(control_table$prior_func[simno])
  ## Pointer to incidence function to use
  inc_func_use <- get(control_table$inc_func[simno])
  ## Maximum allowable age of infection?
  age_max <- control_table$age_max[simno]
  
  ## Most importantly, where is the parameter control table stored?
  parTab <- read.csv("~/Documents/GitHub/virosolver_paper/pars/massachusetts/partab_gp_model.csv", stringsAsFactors=FALSE)
  
  ########################################
  ## 3. Final admin
  ########################################
  ## Fixed control parameters
  n_samp <- 1000
  
  ## CHANGE TO MAIN WD and manage save locations
  #chainwd <- paste0(top_chainwd,runname,"/",run_index,"/timepoint_",timepoint)
  chainwd <- paste0(top_chainwd,"/timepoint_",timepoint)
  plot_wd <- paste0(top_plotwd,runname,"/",run_index)
  setwd(HOME_WD)
  
  ########################################
  ## 4. Model parameters
  ########################################
  ## GP model parameters for fitting
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
  
  ## Read in Ct data
  obs_dat <- read_csv(data_file_cts)
  ## Get possible observation times
  obs_times <- unique(obs_dat$t)
  
  ## Cannot observe after the last observation time
  timepoint <- min(timepoint, length(obs_times))
  
  ## Filter data for either all preceding time points or only this time point
  if(cumu_times) {
    obs_dat_use <- obs_dat %>% filter(t <= obs_times[timepoint])
  } else {
    obs_dat_use <- obs_dat %>% filter(t == obs_times[timepoint])
  }
  ## If only using positive Cts
  if(use_pos) {
    obs_dat_use <- obs_dat_use %>% filter(ct < pars["intercept"])
  }
  
  ## Re-center observation times if we have a maximum age of infection
  if(!is.na(age_max)){
    obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t),t = t + age_max)
  }
  
  ## Vectors of times/infection ages for ismulation
  ages <- 1:max(obs_dat_use$t)
  times <- 0:max(obs_dat_use$t)
  
  ## This is for the GP version
  mat <- matrix(rep(times, each=length(times)),ncol=length(times))
  t_dist <- abs(apply(mat, 2, function(x) x-times))
  if(model_version == "gp"){
    parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
    pars <- parTab$values
    names(pars) <- parTab$names
    
    ## Means for priors
    means <- parTab$values
    names(means) <- parTab$names
  }
  
  ## Epidemic cannot start after first observation time
  parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)
  
  ## Test that model works with these settings
  f <- create_posterior_func(parTab, obs_dat_use, prior_func_use, 
                             inc_func_use,solve_ver="likelihood",
                             use_pos=use_pos,
                             t_dist=t_dist)
  f(pars)
  
  ########################################
  ## 6. Check MCMC run
  ########################################    
  ## Read in chain and remove burn in period
  
  chains_check <- load_mcmc_chains(chainwd,parTab,unfixed=TRUE,multi=FALSE,
                             burnin=mcmcPars_ct["adaptive_period"],
                             chainNo=FALSE)
  print(paste0("Rerun? ", gelman_diagnostics(chains_check$list)$Rerun))
  
  chains <- load_mcmc_chains(chainwd,parTab,unfixed=FALSE,multi=FALSE,
                             burnin=mcmcPars_ct["adaptive_period"],
                             chainNo=FALSE)
  chain <- as.data.frame(chains$chain)
  chain$sampno <- 1:nrow(chain)
  chain1 <- chain
  chain_comb <- chain
  #chain_comb <- chain_comb[,c("sampno",colnames(chain_comb)[!(colnames(chain_comb) %in% c("sampno", "chain"))])]
  #colnames(chain_comb) <- c("sampno",parTab$names,"chain")
  
  ## Get smoothed growth rates and raw incidence projections
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs_smoothed <- matrix(0, nrow=n_samp,ncol=length(times))
  trajs_normal <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    ## Ensure these are all positive
    trajs_smoothed[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
    trajs_normal[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  }
  ## Get daily growth rate from smoothed version, find quantiles
  trajs1 <- t(apply(trajs_smoothed, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + 
    coord_cartesian(ylim=c(-0.5,0.5))
  
  ## Get daily incidence from unsmoothed version, find quantiles
  trajs_quants <- t(apply(trajs_normal, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower_95","lower_50","median","upper_50","upper_95","t")
  
  ## Ct distribution plot
  p_dat <- ggplot(obs_dat_use %>% filter(ct < pars["intercept"])) + 
    geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",
                draw_quantiles=c(0.025,0.5,0.975)) + 
    geom_dotplot(aes(x=t, y=ct,group=t),binaxis="y",
                 binwidth=1,stackdir="center",binpositions="all",dotsize=0.25) +
    scale_y_continuous(trans="reverse") +
    export_theme +
    scale_x_continuous(limits=c(0,200))
  
  seir_outputs <- read_csv("~/Documents/GitHub/virosolver_paper/data/MA_SEIR_sim/sim_MA_gp_1_seir_outputs.csv")
  trajs_quants$frame <- timepoint
  trajs_quants_all <- bind_rows(trajs_quants_all, trajs_quants)
  ## Incidence rate plot
  
  obs_dat_use <- obs_dat_use %>% filter(ct < pars["intercept"])
  
  p_dat <-  ggplot(obs_dat_use) + 
    geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",
                draw_quantiles=c(0.025,0.5,0.975)) + 
    geom_dotplot(aes(x=t, y=ct,group=t),binaxis="y",
                 binwidth=1,stackdir="center",binpositions="all",dotsize=0.25) +
    scale_y_continuous(trans="reverse") +
    export_theme +
    scale_x_continuous(limits=c(0,220))
  
  obs_dat_use$frame <- timepoint
  
  obs_dat_use_all <- bind_rows(obs_dat_use_all, obs_dat_use)
  
  
  p_inc <- ggplot(trajs_quants) +
    #geom_line(data=obs_dat %>% mutate(is_pos=ct < 40) %>% group_by(t) %>% summarize(detect=sum(is_pos)/n()),
    #          aes(x=t,y=detect)) +
    geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25,fill=AAAS_palette["blue1"]) + 
    geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5,fill=AAAS_palette["blue1"]) + 
    geom_line(aes(x=t,y=median)) + 
    #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
    geom_line(data=tibble(t=1:250,y=(seir_outputs$inc/1000000)[1:250]),aes(x=t,y=y),col=AAAS_palette["red1"]) +
    export_theme +
    ylab("Per capita incidence") +
    xlab("Days since start") +
    #scale_x_continuous(limits=c(0,220)) #+
    coord_cartesian(ylim=c(0,0.002))
   ps[[index-3]] <- p_inc
}

p_dat <-  ggplot(obs_dat_use_all) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=t, y=ct,group=t),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.25) +
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_continuous(limits=c(0,210)) +
  labs(x = 'Days since start', y = 'Ct value') +
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10)) +
  transition_time(frame) +
  ease_aes('linear')

p_comb <- ggplot(trajs_quants_all) + 
  geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=t,y=median),col=AAAS_palette["blue1"]) + 
  geom_line(data=tibble(t=1:200,y=(seir_outputs$inc/1000000)[1:200]),aes(x=t,y=y),col=AAAS_palette["red1"]) + 
  geom_vline(data=obs_dat_use_all,aes(xintercept=t),linetype="dashed") +
  export_theme +
  ylab("Per capita incidence") +
  xlab("Days since start") +
  scale_x_continuous(limits=c(0,210)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.0025)) +
  labs(x = 'Days since start', y = 'Per capita incidence') +
  transition_time(frame) +
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10)) +
  ease_aes('linear')

a_gif <- animate(p_dat, width = 8, height = 3,units="in",res=300,duration=15,end_pause=5,start_pause=5)
b_gif <- animate(p_comb, width = 8, height = 3,units="in",res=300,duration=15,end_pause=5,start_pause=5)

a_mgif <- image_read(a_gif)
b_mgif <- image_read(b_gif)

new_gif <- image_append(c(a_mgif[1], b_mgif[1]),stack=TRUE)
for(i in 2:100){
  combined <- image_append(c(a_mgif[i], b_mgif[i]),stack=TRUE)
  new_gif <- c(new_gif, combined)
}

anim_save("~/Documents/GitHub/virosolver_paper/figures/sim_recovery_allcts.gif",new_gif)

