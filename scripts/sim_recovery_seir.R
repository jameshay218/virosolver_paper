########################################
## SCRIPT TO SIMULATE LINE LIST DATA AND CT VALUES 
## FROM AN SEIR PROCESS, THEN RE-ESTIMATE Rt
## Purpose is to show value added beyond the prior
## 31st March 2021
## ------ Structure ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence curve and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
#' 4. Simulate Ct values for the sub-setted line list, giving a Ct value for each observed individual
#' 5. Feed these Ct values to the inference framework and re-estimate the simulated epidemic trajectory

########################################
## 1. Headers
########################################
## Install these packages
library(tidyverse)
library(ggplot2)
library(extraDistr)
library(patchwork)
library(ggthemes)
library(odin) ## install from CRAN
library(doParallel)
library(fitdistrplus)
library(lazymcmc)

HOME_WD <- "~/Documents/GitHub/"

## Where to perform the simulations
main_wd <- paste0(HOME_WD,"/virosolver_paper")
setwd(main_wd)
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
## Probably fine if you just leave these unchanged.
source("code/priors.R")

export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                    axis.text.y=element_text(size=7),
                                    axis.title.x=element_text(size=8),
                                    axis.title.y=element_text(size=8),
                                    strip.text=element_text(face="bold",size=10)
                                    )

## What to call this run
run_name <- paste0("sim_recover_seir")

## Where to save the simulated data
save_wd <- paste0(HOME_WD, "/virosolver_paper/data/sim_recovery_seir/")
save_file <- paste0(save_wd, run_name)
if(!file.exists(save_wd)) dir.create(save_wd,recursive = TRUE)

## Where to save MCMC chains and outputs
chainwd_seir <- paste0(main_wd,"/mcmc_chains/",run_name)
plot_wd <- paste0(main_wd,"/plots/",run_name)
if(!file.exists(chainwd_seir)) dir.create(chainwd_seir,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

## Arguments for this run
set.seed(2)

n_samp <- 1000 ## How many posterior samples to use for plots etc
rerun_seir <- FALSE

## NOTE THIS BIT TO USE FOREACH LOOP
## Manage MCMC runs and parallel runs
n_clusters <- 10
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

########################################
## 2. Model parameters and simulation settings
########################################
## Simulation parameters
population_n <- 100000
## Over the course of 250 days
times <- 0:250

## Sampling procedure - how often do we take samples from the population?
sampling_frequency <- 14
## How many samples do we take on each sample day?
sampling_number <- 1000

## Model parameters
model_pars <- read.csv("pars/nursing_homes/partab_seir_model.csv")
model_pars[model_pars$names == "I0","values"] <- 10/population_n 
#View(model_pars)
pars <- model_pars$values
names(pars) <- model_pars$names

## IMPORTANT CHECK
## - Check that the assumed viral kinetics are in line
##   with your data. This means things like peak Ct value,
##   waning rate, level of the "plateau" phase, and the 
##   limit of detection
test_ages <- seq(1,50,by=1)
cts <- viral_load_func(pars, test_ages)
prop_detect <- prop_detectable(test_ages,pars, cts)
p1 <- ggplot(data.frame(ct=cts,t=test_ages)) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,0)) +
  ylab("Modal Ct value") +
  xlab("Days since infection")
p2 <- ggplot(data.frame(p=prop_detect,t=test_ages)) + 
  geom_line(aes(x=t,y=p)) + 
  ylab("Proportion of infections still detectable") +
  xlab("Days since infection")
p1/p2

########################################
## 3. Full simulated line list
########################################
## Simulate SEIR dynamics, incidence, growth rates, Rt etc
## Use "ode" for deterministic
## Use "odin" for stochastic model
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, ver="odin",switch_model=FALSE)

## Look at simulated dynamics over time
seir_dynamics$plot

## Basic plot, showing the proportion of the population infected each day.
## This is the epidemic incidence curve that you're trying to estimate.
## So you can substitute this object for **any** incidence curve that
## you want to estimate.
plot(seir_dynamics$incidence/population_n)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
## These entries give the **true** infection (or lack of infection) timings for each individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                               population_n=population_n)

########################################
## 4. Simulate observation process
########################################
## This bit of code generates linelist data for some specified observation process.
## Here, we simulate sampling 2000 people at random from the population every 14 days
## arguments changed above
sample_probs <- c(rep(0, sampling_frequency-1),sampling_number/population_n)
sample_probs <- rep(sample_probs, length(times)/sampling_frequency +1)
sample_probs <- sample_probs[1:length(times)]
frac_report <- tibble(t=times,prob=sample_probs)
frac_report <- frac_report %>% filter(t >= 50 & t <= 160)

## frac_report is a table, giving the proportion (prob) of the population
## sampled on day t
head(frac_report)

## This function takes the complete linelist and sub-samples a proportion
## of the entire population on the days specified by `frac_report`
observed_linelist <- simulate_reporting(complete_linelist, 
                                        frac_report=NULL,
                                        timevarying_prob=frac_report,
                                        solve_times=times, 
                                        symptomatic=FALSE)

## Simulate viral load/Ct value for each person in the line list
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,
                                                      kinetics_pars=pars)

## Clean data to form expected by virosolver
obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)

## Save simulated line list, Cts and SEIR dynamics
write_csv(seir_dynamics$seir_outputs, path=paste0(save_file,"_seir_outputs.csv"))
write_csv(complete_linelist, path=paste0(save_file,"_full_linelist.csv"))
write_csv(obs_dat, path=paste0(save_file,"_cts.csv"))

## Save SEIR plots
## Ct distribution plot
p_dat_all <- ggplot(obs_dat %>% filter(t < 125) %>%  filter(ct < pars["intercept"])) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_continuous(limits=c(min(times),max(times)+50)) +
  ylab("Ct value") +
  xlab("Observation time")

ggsave(paste0(save_file,"_ct_dist.png"),p_dat_all,width=7,height=4,dpi=150)
ggsave(paste0(save_file,"_seir_model.png"),seir_dynamics$plot,width=7,height=7,dpi=150)
ggsave(paste0(save_file,"_growth_rates.png"),seir_dynamics$growth_rate_p,width=7,height=6,dpi=150)

head(obs_dat)
obs_dat <- obs_dat %>% filter(t < 150)
ages <- 1:max(obs_dat$t)
times <- 0:max(obs_dat$t)


########################################
## 5. Fit the SEIR model to the simulated data
########################################
## REMEMBER, obs_dat is your simulated dataset.
## Your aim is to use this dataset to re-estimate the
## contents of seir_dynamics.
## Manage MCMC runs and parallel runs
n_clusters <- 10
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

nchains <- 3
n_temperatures <- 5
mcmcPars_ct <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                    "thin"=10,"adaptive_period"=20000,"save_block"=1000,"temperature" = seq(1,101,length.out=n_temperatures),
                    "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                    "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)
inc_func_use <- solveSEIRModel_rlsoda_wrapper
prior_func_use <- prior_func_hinge_seir

all_grs <- NULL
all_rts <- NULL
all_incs <- NULL
all_dat <- NULL

for(run_version in c("prior_only","pos_only","all_cts")){
  print(run_version)
  ## Where to save MCMC chains and outputs
  chainwd_tmp <- paste0(main_wd,"/mcmc_chains/",run_name,"/",run_version)
  plot_wd_tmp <- paste0(main_wd,"/plots/",run_name,"/",run_version)
  if(!file.exists(chainwd_tmp)) dir.create(chainwd_tmp,recursive = TRUE)
  if(!file.exists(plot_wd_tmp)) dir.create(plot_wd_tmp,recursive = TRUE)
  
  solve_likelihood <- TRUE
  use_pos <- FALSE
  if(run_version == "pos_only"){
    use_pos <- TRUE
  }
  if(run_version == "prior_only"){
    solve_likelihood <- FALSE
  }
  
  ## GP model parameters for fitting
  parTab <- read.csv("pars/nursing_homes/partab_seir_model.csv")
  if(!use_pos) {
    parTab[parTab$names == "overall_prob","fixed"] <- 0
  }
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
  
  obs_times <- unique(obs_dat$t)
  ages <- 1:max(obs_dat$t)
  times <- 0:max(obs_dat$t)
  
  ## Check that posterior function solves correctly
  f <- create_posterior_func(parTab, obs_dat, prior_func_use, 
                             inc_func_use,solve_ver="likelihood",
                             use_pos=use_pos,
                             t_dist=t_dist)
  f(pars)
  
  
  if(rerun_seir){
    res <- foreach(i=seq_along(obs_times),.packages = c("extraDistr","tidyverse","patchwork","virosolver")) %dopar% {
      ## Need to read in R package
      devtools::load_all(paste0(HOME_WD,"/lazymcmc"))
      
      timepoint <- obs_times[i]
      runname_use <- paste0(run_name,"_time_",timepoint)
      dir.create(paste0(chainwd_tmp,"/",timepoint),recursive = TRUE)
      obs_dat_use <- obs_dat %>% filter(t == timepoint)
      
      if(use_pos){
        obs_dat_use <- obs_dat_use %>% filter(ct < 40)
      }
      
      ages <- 1:max(obs_dat_use$t)
      times <- 0:max(obs_dat_use$t)
      parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)
      
      chains <- NULL
      for(j in 1:nchains){
        if(n_temperatures > 1){
          startTab <- rep(list(parTab),n_temperatures)
          for(k in 1:length(startTab)){
            startTab[[k]] <- generate_viable_start_pars(parTab,obs_dat_use,
                                                        create_posterior_func,
                                                        inc_func_use,
                                                        prior_func_use,
                                                        t_dist=NULL,
                                                        use_pos=use_pos)
            }
        }
        
        covMat <- diag(nrow(startTab[[1]]))
        mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
        mvrPars <- rep(list(mvrPars), n_temperatures)
        output <- run_MCMC(parTab=startTab,
                           data=obs_dat_use,
                           INCIDENCE_FUNC=inc_func_use,
                           PRIOR_FUNC = prior_func_use,
                           solve_likelihood=solve_likelihood,
                           mcmcPars=mcmcPars_ct,
                           filename=paste0(chainwd_tmp,"/",timepoint,"/",runname_use,"_chainno_",j),
                           CREATE_POSTERIOR_FUNC=create_posterior_func,
                           mvrPars=mvrPars,
                           OPT_TUNING=0.2,
                           use_pos=use_pos,
                           t_dist=NULL)
        
        ## Read in chain and remove burn in period
        chain <- read.csv(output$file)
        chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
        chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
        chain$chain <- j
        chains[[j]] <- chain
      }
      chain <- do.call("bind_rows",chains)
    }
  } 
  
  ########################################
  ## 6. Generate plots from the estimated posterior
  ########################################
  res <- NULL
  res_trajs <- NULL
  res_rt <- NULL
  res_gr <- NULL
  
  for(i in seq_along(obs_times)){
    timepoint <- obs_times[i]
    chainwd_tmp2 <- paste0(chainwd_tmp,"/",timepoint)
    chain <- lazymcmc::load_mcmc_chains(chainwd_tmp2, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                        multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
    chain <- as.data.frame(chain)
    res[[i]] <- chain
  }
  
  ## Go through each time point that we estimated dynamics from
  ## Get the MCMC chains ie. the posteriors and generate some plots
  ## and analyses
  for(i in seq_along(obs_times)){
    timepoint <- obs_times[i]
    runname_use <- runname_use <- paste0(run_name,"_time_",timepoint)
    
    obs_dat_tmp <- obs_dat_use <- obs_dat %>% filter(t == timepoint)
    
    ages <- 1:max(obs_dat_use$t)
    times <- 0:max(obs_dat_use$t)
    
    chain <- res[[i]]
    chain_comb <- chain
    chain_comb$sampno <- 1:nrow(chain_comb)
    chain1 <- chain
    chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
    
    
    ## Check convergence of MCMC chains. 
    ## We want to see lots of hairy caterpillars.
    ## Some of the parameters will be flat lines - these
    ## were the parameters that we fixed.
    p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
      mutate(chain = as.factor(chain)) %>%
      pivot_longer(-c(sampno,chain)) %>%
      ggplot() +
      geom_line(aes(x=sampno,y=value,col=chain)) +
      facet_wrap(~name,scales="free_y")+
      scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),length.out=5)) +
      export_theme
    
    ## Get smoothed growth rates
    samps <- sample(unique(chain_comb$sampno),n_samp)
    dat_inc <- NULL 
    
    for(ii in seq_along(samps)){
      pars <- get_index_pars(chain_comb, samps[ii])
      inc_est <- inc_func_use(pars,times)[1:length(times)]
      #inc_est <- inc_est/sum(inc_est)
      tmp_inc <- tibble(t=times,inc=inc_est)
      x <- tmp_inc$inc
      
      ## Get Rt
      seir_pars <- c("beta"=pars["R0"]*(1/pars["infectious"]),1/pars["incubation"],1/pars["infectious"])
      names(seir_pars) <- c("beta","sigma","gamma")
      init <- c(1-pars["I0"],0,pars["I0"],0)
      sol <- solveSEIRModel_rlsoda(times, init, seir_pars,compatible=TRUE)
      S_compartment <- c(rep(1, pars["t0"]),sol[,2])[1:length(times)]
      tmp_inc$Rt <- S_compartment*unname(pars["R0"])
      tmp_inc$gr <- c(NaN,log(x[2:length(x)]/x[1:(length(x)-1)]))
      
      dat_inc[[ii]] <- tmp_inc
      dat_inc[[ii]]$samp <- ii
    }
    trajs <- do.call("bind_rows",dat_inc)
    trajs1_quants <- trajs %>% group_by(t) %>%
      summarize(lower_95=quantile(gr, 0.025,na.rm=TRUE),
                lower_50=quantile(gr,0.25,na.rm=TRUE),
                median=median(gr,na.rm=TRUE),
                upper_50=quantile(gr,0.75,na.rm=TRUE),
                upper_95=quantile(gr, 0.975,na.rm=TRUE)) %>%
      mutate(obs_time = timepoint,
             run_version=run_version)
    res_gr[[i]] <- trajs1_quants
    
    
    Rt_quants <- trajs %>% group_by(t) %>%
      summarize(lower_95=quantile(Rt, 0.025,na.rm=TRUE),
                lower_50=quantile(Rt,0.25,na.rm=TRUE),
                median=median(Rt,na.rm=TRUE),
                upper_50=quantile(Rt,0.75,na.rm=TRUE),
                upper_95=quantile(Rt, 0.975,na.rm=TRUE)) %>%
      mutate(obs_time = timepoint,
             run_version=run_version)
    res_rt[[i]] <- Rt_quants
    
    ## Growth rate plot
    p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) + 
      geom_line(aes(x=t+model_pars[model_pars$names=="t0","values"],y=median)) + 
      coord_cartesian(ylim=c(-0.5,0.5))
    ## Rt plot
    p_rt <- ggplot(Rt_quants) + 
      geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) + 
      geom_line(aes(x=t,y=median)) +
      geom_line(data=seir_dynamics$seir_outputs %>% filter(step <= max(trajs1_quants$t)),
                aes(x=step+model_pars[model_pars$names=="t0","values"],y=Rt),col="green") +
      scale_x_continuous(limits=c(min(times),max(times)+50)) +
      coord_cartesian(ylim=c(0,10))
    
    inc_quants <- trajs %>% group_by(t) %>%
      summarize(lower_95=quantile(inc, 0.025,na.rm=TRUE),
                lower_50=quantile(inc,0.25,na.rm=TRUE),
                median=median(inc,na.rm=TRUE),
                upper_50=quantile(inc,0.75,na.rm=TRUE),
                upper_95=quantile(inc, 0.975,na.rm=TRUE)) %>%
      mutate(obs_time = timepoint,
             run_version=run_version)
    res_trajs[[i]] <- trajs %>% mutate(obs_time = timepoint,
                                       run_version=run_version)
    
    tmp_inc_real <- seir_dynamics$incidence[seq_along(times)]/population_n
    all_dat[[i]] <- tibble(t=times,y=seir_dynamics$incidence[seq_along(times)]/population_n) %>% filter(t <= timepoint) %>% mutate(obs_time=timepoint)
    
    #tmp_inc_real <- tmp_inc_real/sum(tmp_inc_real)
    
    ## Incidence plot
    p_inc <- ggplot(inc_quants) + 
      geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) + 
      geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5) + 
      geom_line(aes(x=t,y=median)) + 
      geom_line(data=tibble(x=times,y=tmp_inc_real),
                aes(x=x,y=y),col="darkgreen") +
      scale_x_continuous(limits=c(min(times),max(times)+50)) +
      export_theme +
      ylab("Per capita incidence") +
      xlab("Days since start") +
      coord_cartesian(ylim=c(0,0.1))
    
    
    p_dat <- ggplot(obs_dat_tmp %>% filter(ct < pars["intercept"])) + 
      geom_violin(aes(x=t,group=t,y=ct),width=10,fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
      geom_jitter(aes(x=t,y=ct),size=0.25,width=10,height=0) + 
      scale_y_continuous(trans="reverse") +
      export_theme +
      scale_x_continuous(limits=c(min(times),max(times)+50)) +
      ylab("Ct value") +
      xlab("Observation time")
    
    
    vl_trajs <-  matrix(0, nrow=n_samp,ncol=length(ages))
    for(ii in 1:n_samp){
      tmp_pars <- get_index_pars(chain_comb, samps[ii])
      tmp <- viral_load_func(tmp_pars,ages,FALSE)
      tmp1 <- extraDistr::rgumbel(length(tmp),tmp, tmp_pars["obs_sd"])
      vl_trajs[ii,] <- tmp1
    }
    vl_trajs1_quants <- t(apply(vl_trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
    vl_trajs1_quants <- as.data.frame(vl_trajs1_quants)
    vl_trajs1_quants$t <- 1:nrow(vl_trajs1_quants)
    colnames(vl_trajs1_quants) <- c("lower","median","upper","t")
    
    ## Growth rate plot
    p_vl <- ggplot(vl_trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
      geom_line(aes(x=t,y=median))
    
    dir.create(paste0(plot_wd_tmp,"/traces/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/predictions/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/Rt/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/grs/"),recursive = TRUE)
    
    ggsave(paste0(plot_wd_tmp,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
    ggsave(paste0(plot_wd_tmp,"/predictions/",runname_use,"_predictions.png"),p_dat/p_inc,width=7,height=7)
    ggsave(paste0(plot_wd_tmp,"/Rt/",runname_use,"_grs.png"),p_rt,width=7,height=4)
    ggsave(paste0(plot_wd_tmp,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  }
  
  
  all_grs[[run_version]] <- do.call("bind_rows", res_gr)
  all_rts[[run_version]] <- do.call("bind_rows", res_rt)
  all_incs[[run_version]] <- do.call("bind_rows", res_trajs)
}
all_grs <- do.call("bind_rows",all_grs)
all_rts <- do.call("bind_rows",all_rts)
all_incs <- do.call("bind_rows",all_incs)
all_dat <- do.call("bind_rows",all_dat)

run_version_key <- c("all_cts"="All Ct values",
                     "prior_only"="Prior (no data)",
                     "pos_only"="Only detectable Ct values")
all_incs$run_version <- run_version_key[all_incs$run_version]
all_grs$run_version <- run_version_key[all_grs$run_version]
all_rts$run_version <- run_version_key[all_rts$run_version]

all_incs$obs_time <- paste0("t=",all_incs$obs_time)
all_grs$obs_time <- paste0("t=",all_grs$obs_time)
all_rts$obs_time <- paste0("t=",all_rts$obs_time)
all_dat$obs_time <- paste0("t=",all_dat$obs_time)

all_incs$obs_time <- factor(all_incs$obs_time, levels=paste0("t=",obs_times))
all_grs$obs_time <- factor(all_grs$obs_time, levels=paste0("t=",obs_times))
all_rts$obs_time <- factor(all_rts$obs_time, levels=paste0("t=",obs_times))
all_dat$obs_time <- factor(all_dat$obs_time, levels=paste0("t=",obs_times))

p_all_cts <- all_incs %>% 
  filter(!(obs_time %in% c("t=125","t=139"))) %>% 
  filter(run_version != "Only detectable Ct values") %>% 
  group_by(t, run_version, obs_time) %>%
  summarize(lower_95=quantile(inc, 0.025,na.rm=TRUE),
            lower_50=quantile(inc,0.25,na.rm=TRUE),
            median=median(inc,na.rm=TRUE),
            upper_50=quantile(inc,0.75,na.rm=TRUE),
            upper_95=quantile(inc, 0.975,na.rm=TRUE)) %>%
  ggplot() + 
  geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25,fill=AAAS_palette["blue1"]) +
  geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5,fill=AAAS_palette["blue1"]) +
  geom_line(aes(x=t,y=median),col=AAAS_palette["blue1"]) +
  geom_line(data=all_dat %>% group_by(obs_time) %>% filter(!(obs_time %in% c("t=125","t=139"))),aes(x=t,y=y),col=AAAS_palette["red1"],size=0.5) +
  export_theme +
  theme(panel.grid.major = element_line(color="grey40",size=0.1),
        panel.border = element_rect(color="black",fill=NA,size=0.2)) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.1,by=0.02)) +
  coord_cartesian(ylim=c(0,0.1)) +
  scale_x_continuous(limits=c(min(times),125),expand=c(0,0),breaks=seq(0,125,by=25)) +
  xlab("Time") +
  ylab("Incidence") +
  facet_grid(obs_time~run_version) +
  labs(tag="B")

p_pos_cts <- all_incs %>% 
  filter(!(obs_time %in% c("t=125","t=139"))) %>% 
  filter(run_version != "All Ct values") %>%
  group_by(run_version, obs_time, samp) %>%
  mutate(inc=inc/sum(inc)) %>%
  group_by(t, run_version, obs_time) %>%
  summarize(lower_95=quantile(inc, 0.025,na.rm=TRUE),
            lower_50=quantile(inc,0.25,na.rm=TRUE),
            median=median(inc,na.rm=TRUE),
            upper_50=quantile(inc,0.75,na.rm=TRUE),
            upper_95=quantile(inc, 0.975,na.rm=TRUE)) %>%
  ggplot() + 
  geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25,fill=AAAS_palette["blue1"]) +
  geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5,fill=AAAS_palette["blue1"]) +
  geom_line(aes(x=t,y=median),col=AAAS_palette["blue1"]) +
  geom_line(data=all_dat %>% group_by(obs_time) %>%mutate(y=y/sum(y)) %>% filter(!(obs_time %in% c("t=125","t=139"))),aes(x=t,y=y),col=AAAS_palette["red1"]) +
  facet_grid(obs_time~run_version) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.2,by=0.05)) +
  scale_x_continuous(limits=c(min(times),125),expand=c(0,0),breaks=seq(0,125,by=25)) +
  export_theme +
  theme(panel.grid.major = element_line(color="grey40",size=0.1),
        panel.border = element_rect(color="black",fill=NA,size=0.2)) +
  xlab("Time") +
  ylab("Relative probability of infection") +
  coord_cartesian(ylim=c(0,0.2))+
  labs(tag="C")

all_real_seir <- NULL
for(obs_time in obs_times){
  tmp_seir <- seir_dynamics$seir_outputs %>% mutate(obs_time=obs_time) %>% filter(step <= obs_time)
  all_real_seir[[obs_time]] <- tmp_seir
}
all_real_seir <- do.call("bind_rows", all_real_seir)
all_rts %>% 
  filter(!(obs_time %in% c("t=125","t=139"))) %>% 
  ggplot() + 
  geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  geom_line(data=all_real_seir%>% filter(obs_time < 125),aes(x=step,y=Rt),col="darkgreen") +
  facet_grid(obs_time~run_version)

all_grs %>% 
  filter(!(obs_time %in% c("t=125","t=139"))) %>% 
  ggplot() + 
  geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_50,ymax=upper_50),alpha=0.5) +
  facet_grid(obs_time~run_version)


percent_pos <- obs_dat %>% filter(t <= 111) %>% 
  group_by(t) %>% 
  mutate(is_pos=ct < pars["intercept"]) %>% 
  summarize(percent_pos=100*sum(is_pos)/n()) %>%
  mutate(percent_pos =paste0(sprintf(percent_pos,fmt='%#.1f'),"%"))

obs_dat_lhs <- obs_dat %>% filter(t <= 111) %>%  filter(ct < pars["intercept"]) %>% mutate(run_version="All Ct values")
obs_dat_rhs <- obs_dat %>% filter(t <= 111) %>%  filter(ct < pars["intercept"]) %>% mutate(run_version="Prior (no data)")
obs_dat_lhs <- bind_rows(obs_dat_lhs,obs_dat_rhs)

coeff1 <- 400
p_dat_lhs <- ggplot(seir_dynamics$seir_outputs%>% filter(step <= 111)) + 
  geom_violin(data=obs_dat_lhs,aes(x=t,group=t,y=(40-ct)/coeff1),scale="width",fill="grey70",
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(data=obs_dat_lhs,aes(x=t,y=(40-ct)/coeff1),size=0.01,width=2,height=0) + 
  geom_line(aes(x=step,y=inc/population_n),col=AAAS_palette["red1"],size=1) +
  geom_text(data=percent_pos, aes(x=t,y=0.09,label=percent_pos),size=3) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.1,by=0.02),
                     sec.axis=sec_axis(~.*-coeff1 + 40, name="Ct value")) +
  coord_cartesian(ylim=c(0,0.1)) +
  export_theme +
  theme(panel.grid.major = element_line(color="grey40",size=0.1),
        panel.border = element_rect(color="black",fill=NA,size=0.2),
        strip.text=element_blank()) +
  scale_x_continuous(limits=c(min(times),125),expand=c(0,0),breaks=seq(0,125,by=25)) +
  ylab("Incidence") +
  xlab("") + 
  facet_wrap(~run_version) +
  labs(tag="A")

obs_dat_lhs1 <- obs_dat %>% filter(t <= 111) %>%  filter(ct < pars["intercept"]) %>% mutate(run_version="Only detectable Ct values")
obs_dat_rhs1 <- obs_dat %>% filter(t <= 111) %>%  filter(ct < pars["intercept"]) %>% mutate(run_version="Prior (no data)")
obs_dat_lhs1 <- bind_rows(obs_dat_lhs1,obs_dat_rhs1)

coeff2 <- 200
p_dat_rhs <- ggplot(seir_dynamics$seir_outputs%>% filter(step <= 111)%>%
                      mutate(inc=inc/sum(inc))) + 
  geom_violin(data=obs_dat_lhs1,aes(x=t,group=t,y=(40-ct)/coeff2),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(data=obs_dat_lhs1,aes(x=t,y=(40-ct)/coeff2),size=0.01,width=2,height=0) + 
  geom_line(aes(x=step,y=inc),col=AAAS_palette["red1"],size=1) +
  geom_text(data=percent_pos, aes(x=t,y=0.18,label=percent_pos),size=3) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.2,by=0.05),
                     sec.axis=sec_axis(~.*-coeff2 + 40, name="Ct value")) +
  coord_cartesian(ylim=c(0,0.2)) +
  export_theme +
  theme(panel.grid.major = element_line(color="grey40",size=0.1),
        panel.border = element_rect(color="black",fill=NA,size=0.2),
        strip.text=element_blank()) +
  scale_x_continuous(limits=c(min(times),125),expand=c(0,0),breaks=seq(0,125,by=25)) +
  ylab("Relative probability of infection") +
  xlab("") + 
  facet_wrap(~run_version)

bot_p <- p_all_cts | p_pos_cts
top_p <- p_dat_lhs | p_dat_rhs

main_p <- (top_p / bot_p) + plot_layout(heights=c(1,4))
ggsave("figures/supplement/cross_section_fits_raw.png",plot=main_p,height=8,width=10,units="in",dpi=300)
ggsave("figures/supplement/cross_section_fits_raw.pdf",plot=main_p,height=8,width=10)
