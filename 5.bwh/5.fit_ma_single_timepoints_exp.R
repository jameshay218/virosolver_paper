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
library(doParallel)
#library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")

HOME_WD <- "~/Documents/GitHub/"

devtools::load_all(paste0(HOME_WD,"virosolver"))
devtools::load_all(paste0(HOME_WD,"lazymcmc"))

## Arguments for this run
set.seed(1)
n_samp <- 1000
runname <- "ma_exp"
run_version <- "exp" ##gp, seir or exp##

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc <- TRUE

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- paste0(HOME_WD,"/virosolver_paper/")
chainwd <- paste0(HOME_WD, "/virosolver_paper/mcmc_chains/5.real_ma_single_timepoint/",runname,"/")
plot_wd <- paste0(HOME_WD, "/virosolver_paper/plots/5.real_ma_single_timepoint/",runname,"/")
setwd(main_wd)

## Manage MCMC runs and parallel runs
nchains <- 3
n_clusters <- 11
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

## MCMC parameters for Ct model fits
## MCMC control parameters
n_temperatures <- 10
mcmcPars_ct <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                 "thin"=10,"adaptive_period"=30000,"save_block"=1000,"temperature" = seq(1,101,length.out=n_temperatures),
                 "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                 "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)

## Code for plotting
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)


########################################
## 2. Model parameters and simulation settings
########################################
#max_age <- NA
max_age <- 35
inc_func_use <- exponential_growth_model
prior_func_use <- prior_func_hinge_exp
#inc_func_use <- solveSEIRModel_rlsoda_wrapper
#prior_func_use <- prior_func_hinge_seir

## GP model parameters for fitting
#parTab <- read.csv("pars/massachusetts/partab_seir_model.csv")
parTab <- read.csv("pars/massachusetts/partab_exp_pos_model.csv")
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

## Simulation parameters
times <- 0:200
## Extend to account for delays
times_extended <-c(times,max(times):(max(times)+50)) 


########################################
## 3. Read in MA data
########################################
#obs_dat_all <- read_csv("~/Documents/GitHub/ct_inference_preprint/data/BWH_COVID_Cts_deid_20200403-20200831.csv") %>%
#  mutate(id=1:n())
obs_dat_all <- read_csv("data/panther_Ct_20200403-20201110.csv") %>% rename(panther_Ct=ORF1ab_Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n())

obs_dat1 <- obs_dat_all

obs_dat1 <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct)

obs_times <- unique(obs_dat1$t)

obs_dat_all <- obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  arrange(first_day) %>%
  rename(date = first_day, ct=panther_Ct)

comb_dat <- left_join(obs_dat1, obs_dat_all)
date_key <- distinct(comb_dat %>% dplyr::select(t, date))

date_min_date <- min(date_key$date)
date_min_t <- min(date_key$t)

date_max_date <- max(date_key$date)
date_max_t <- max(date_key$t)

integer_seq_times <- seq(0, date_max_t)
date_seq_times <- seq(date_min_date-date_min_t, date_max_date,by="1 day")
date_key <- tibble(t=integer_seq_times,date=date_seq_times)

p_dat <- ggplot(obs_dat_all %>% filter(ct < 40)) + 
  geom_violin(aes(x=date,group=date,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_date(limits=as.Date(c("2020-01-01","2020-10-01")),breaks="1 month") +
  xlab("Date of sample") +
  ylab("Detectable Ct")

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
#parTab[!(parTab$names %in% c("prob","rho","nu","obs_sd")),"fixed"] <- 1

f <- create_posterior_func(parTab, obs_dat1, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=TRUE,
                           t_dist=t_dist)
f(pars)
## Run for each chain
if(rerun_mcmc){
  res <- foreach(i=seq_along(obs_times),.packages = c("extraDistr","tidyverse","patchwork")) %dopar% {
    ## Need to read in R package
    devtools::load_all(paste0(HOME_WD,"/virosolver"))
    devtools::load_all(paste0(HOME_WD,"/lazymcmc"))
    timepoint <- obs_times[i]
    runname_use <- paste0(runname,"_time_",timepoint)
    dir.create(paste0(chainwd,"/",timepoint),recursive = TRUE)
    
    obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
      
    ## Observation times
    if(!is.na(max_age)){
      obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
    }
    
    ages <- 1:max(obs_dat_use$t)
    times <- 0:max(obs_dat_use$t)
    
    parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)
    
    chains <- NULL
    for(j in 1:nchains){
      ## Get random starting values
      #startTab <- generate_viable_start_pars(parTab,obs_dat_use,
      #                                       create_posterior_func,
      #                                       inc_func_use,
      #                                       prior_func_use,
      #                                       t_dist=NULL,
      #                                       use_pos=TRUE)
      if(n_temperatures > 1){
        startTab <- rep(list(parTab),n_temperatures)
        for(k in 1:length(startTab)){
          startTab[[k]] <- generate_viable_start_pars(parTab,obs_dat_use,
                                                      create_posterior_func,
                                                      inc_func_use,
                                                      prior_func_use,
                                                      t_dist=NULL,
                                                      use_pos=TRUE)
          startTab[[k]][startTab[[k]]$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)
        }
      }
      
      covMat <- diag(nrow(startTab[[1]]))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
      mvrPars <- rep(list(mvrPars), n_temperatures)
      output <- run_MCMC(parTab=startTab,
                         data=obs_dat_use,
                         INCIDENCE_FUNC=inc_func_use,
                         PRIOR_FUNC = prior_func_use,
                         solve_likelihood=TRUE,
                         mcmcPars=mcmcPars_ct,
                         filename=paste0(chainwd,"/",timepoint,"/",runname_use,"_chainno_",j),
                         CREATE_POSTERIOR_FUNC=create_posterior_func,
                         mvrPars=mvrPars,
                         OPT_TUNING=0.2,
                         use_pos=TRUE,
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

res <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  chainwd_tmp <- paste0(chainwd,timepoint)
  chain <- lazymcmc::load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                      multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain <- as.data.frame(chain)
  #chain$sampno <- 1:nrow(chain)
  res[[i]] <- chain
}

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use <- runname_use <- paste0(runname,"_time_",timepoint)
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  
  ## Observation times
  if(!is.na(max_age)){
    obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
  }
  
  ages <- 1:max(obs_dat_use$t)
  times <- 0:max(obs_dat_use$t)
  
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  
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
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    #trajs[ii,] <- pmax(smooth.spline(inc_func_use(get_index_pars(chain_comb, samps[ii]),times))$y,0.0000001)
    trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
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
  
  trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")
  
  ## Growth rate plot
  p_inc <- ggplot(trajs_quants %>% left_join(date_key)) + 
    geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5) + 
    geom_line(aes(x=date,y=median)) + 
    #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
    #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
    export_theme +
    ylab("Per capita incidence") +
    xlab("Days since start") +
    scale_x_date(limits=as.Date(c("2020-01-01","2020-10-01")),breaks="1 month") +
    coord_cartesian(ylim=c(0,0.03))
  
  vl_trajs <-  matrix(0, nrow=n_samp,ncol=length(ages))
  for(ii in 1:n_samp){
    tmp_pars <- get_index_pars(chain_comb, samps[ii])
    #tmp_pars <- get_best_pars(chain_comb)
    tmp <- viral_load_func(tmp_pars,ages,FALSE)
    tmp1 <- extraDistr::rgumbel(length(tmp),tmp, tmp_pars["obs_sd"])
    vl_trajs[ii,] <- tmp1
  }
  #vl_trajs[vl_trajs < -3] <- -3
  vl_trajs1_quants <- t(apply(vl_trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  vl_trajs1_quants <- as.data.frame(vl_trajs1_quants)
  vl_trajs1_quants$t <- 1:nrow(vl_trajs1_quants)
  colnames(vl_trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_vl <- ggplot(vl_trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median))
  
  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  #dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  #dir.create(paste0(plot_wd,"/posteriors/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/grs/"),recursive = TRUE)
  
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",runname_use,"_predictions.png"),p_dat/p_inc,width=7,height=7)
  #ggsave(paste0(plot_wd_tmp,"/distributions/",runname_tmp,"_distributions.png"),p2,
  #       width=(7/5) * length(unique(dat_tmp$t)),height=6)
  #ggsave(paste0(plot_wd_tmp,"/posteriors/",runname_tmp,"_densities.png"),p_densities,width=7,height=4)
  ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  
  
}
