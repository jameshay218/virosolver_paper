########################################
## SCRIPT TO FIT MODEL TO SIMULATED LINE LIST DATA
## 2nd November 2020
## ------ Aims ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame and the desired line list data, a tibble with variables t (sample time) and ct (observed ct value)
#' 3. Options for the Ct model fit
#' 4. Fit model to Ct values, MCMC run
#' 5. Read in MCMC chains, pre-computed or new
#' 6. Plot various fits, posterior diagnostics and predictions from the MCMC fit
########################################

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

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- "~/Documents/GitHub/virosolver_paper/mcmc_chains/TEST.linelists"
plot_wd <- "~/Documents/GitHub/virosolver_paper/plots/TEST.linelists"
run_name <- "testing_linelists"
savewd <- paste0(chainwd, "/", run_name)

setwd(main_wd)

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(savewd)) dir.create(savewd,recursive = TRUE)

## Code for plotting
source("code/plot_funcs.R")
export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                     axis.text.y=element_text(size=7),
                                     axis.title.x=element_text(size=8),
                                     axis.title.y=element_text(size=8))

## Manage parallel runs
n_clusters <- 8
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

########################################
## 2. Read in line list data produced by simulate_linelists.R
########################################
obs_dat <- read_csv("data/sim_linelist.csv")
obs_dat <- obs_dat %>% filter(t == 140)
obs_times <- unique(obs_dat$t)

## For SEIR model, use pars/partab_seir_model.csv
## For exp model all Cts, use pars/partab_exp_model.csv
## For exp model only +ve Cts, use pars/partab_exp_pos_model.csv
## For the GP model, use pars/partab_gp_model.csv
parTab <- read.csv("pars/partab_seir_model.csv",stringsAsFactors=FALSE)
#parTab <- read.csv("pars/partab_exp_model.csv",stringsAsFactors=FALSE)
#parTab <- read.csv("pars/partab_exp_pos_model.csv",stringsAsFactors=FALSE)
#parTab <- read.csv("pars/partab_gp_model.csv",stringsAsFactors=FALSE)

pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

########################################
## 3. Choose models, priors and MCMC control parameters
########################################
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

## Can run this script without rerunning the MCMC
rerun_mcmc_ct <- TRUE

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=10000,"popt"=0.234,"opt_freq"=1000,
                 "thin"=10,"adaptive_period"=5000,"save_block"=1000)
nchains <- 3
n_samp <- 1000 ## Number of posterior samples for plots

## CHOOSE INCIDENCE FUNCTION
is_gp_run <- FALSE
## For exponential growth model
#inc_func_use <- exponential_growth_model
## For SEIR model
inc_func_use <- solveSEIRModel_rlsoda_wrapper
## For GP model
#inc_func_use <- gaussian_process_model
#is_gp_run <- TRUE

## CHOOSE PRIOR FUNCTION
## ** SEE code/priors.R **
prior_func_use <- prior_func_hinge_seir

## Independent cross sections or combined?
cumu_data <- TRUE

## Only detectable or all Cts?
only_pos <- FALSE

## How old can infections be? Try 35 for single cross section models.
## For full SEIR model, set this to NA and the code will use the start time of the epidemic
max_age <- NA

########################################
## 4. Fit model to Ct values
########################################
if(rerun_mcmc_ct){
  ## For each observation time, fit nchain chains
  res <- foreach(i=seq_along(obs_times),.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
    devtools::load_all("~/Documents/GitHub/virosolver")
    obs_time <- obs_times[i]
    
    runname_use <- paste0(run_name,"_",obs_time)
    dir.create(paste0(savewd,"/",runname_use),recursive = TRUE)
    
    ## If using independent cross-sections or cumulative data
    if(cumu_data) {
      obs_dat_tmp <- obs_dat %>% filter(t <= obs_time)
    } else {
      obs_dat_tmp <- obs_dat %>% filter(t == obs_time)
    }
    
    ## If only using detectable Cts
    if(only_pos){
      obs_dat_tmp <- obs_dat_tmp %>% filter(ct < pars["intercept"])
    }
    
    dat_tmp_used_in_run <- obs_dat_tmp
    
    ## Observation times
    if(!is.na(max_age)){
      dat_tmp_used_in_run <- dat_tmp_used_in_run %>% mutate(t = t - max(t),
                                                            t = t + max_age)
    }
    ## Only from start to observation time
    ages <- 1:max(dat_tmp_used_in_run$t)
    times <- 0:max(dat_tmp_used_in_run$t)
    
    ## This is for the GP version
    mat <- matrix(rep(times, each=length(times)),ncol=length(times))
    t_dist <- abs(apply(mat, 2, function(x) x-times))
    
    if(is_gp_run){
      parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
      pars <- parTab$values
      names(pars) <- parTab$names
      
      ## Means for priors
      means <- parTab$values
      names(means) <- parTab$names
    }
    
    ## Epidemic cannot start after first observation time
    parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_tmp$t)
    
    ## Run for each chain
    chains <- NULL
    for(j in 1:nchains){
      ## Get random starting values
      startTab <- generate_viable_start_pars(parTab,obs_dat_tmp,
                                             create_posterior_func,
                                             inc_func_use,
                                             prior_func_use)
      covMat <- diag(nrow(startTab))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
      
      output <- run_MCMC(parTab=startTab,
                         data=obs_dat_tmp,
                         INCIDENCE_FUNC=inc_func_use,
                         PRIOR_FUNC = prior_func_use,
                         solve_likelihood=TRUE,
                         mcmcPars=mcmcPars_ct,
                         filename=paste0(savewd,"/", runname_use,"/",runname_use,"_chainno_",j),
                         CREATE_POSTERIOR_FUNC=create_posterior_func,
                         mvrPars=mvrPars,
                         OPT_TUNING=0.2,
                         use_pos=only_pos,
                         t_dist=t_dist)
      
      ## Read in chain and remove burn in period
      chain <- read.csv(output$file)
      chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
      chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
      chains[[j]] <- chain
    }
    chain <- do.call("bind_rows",chains)
  }
}

########################################
## 5. Load in MCMC chains
########################################
res <- NULL
for(i in seq_along(obs_times)){
  runname_use <- paste0(run_name,"_",obs_time)
  chainwd_tmp <- paste0(savewd,"/",runname_use)
  chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=TRUE)$chain
  chain <- as.data.frame(chain)
  res[[i]] <- chain
}

########################################
## 6. All plots from this run
########################################
## Plot useful summary plots
for(i in seq_along(obs_times)){
  runname_use <- paste0(run_name,"_",obs_time)
  obs_time <- obs_times[i]
  
  ## Set up data as if doing the full fitting
  ## If using independent cross-sections or cumulative data
  if(cumu_data) {
    obs_dat_tmp <- obs_dat %>% filter(t <= obs_time)
  } else {
    obs_dat_tmp <- obs_dat %>% filter(t == obs_time)
  }
  
  ## If only using detectable Cts
  if(only_pos){
    obs_dat_tmp <- obs_dat_tmp %>% filter(ct < pars["intercept"])
  }
  
  dat_tmp_used_in_run <- obs_dat_tmp
  
  ## Observation times
  if(!is.na(max_age)){
    dat_tmp_used_in_run <- dat_tmp_used_in_run %>% mutate(t = t - max(t),
                                                          t = t + max_age)
  }
  ## Only from start to observation time
  ages <- 1:max(dat_tmp_used_in_run$t)
  times <- 0:max(dat_tmp_used_in_run$t)
  
  ## This is for the GP version
  mat <- matrix(rep(times, each=length(times)),ncol=length(times))
  t_dist <- abs(apply(mat, 2, function(x) x-times))
  
  if(is_gp_run){
    parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
    pars <- parTab$values
    names(pars) <- parTab$names
    
    ## Means for priors
    means <- parTab$values
    names(means) <- parTab$names
  }
  
  
  ## Pull pre-read chain and get this observation time
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  
  ## Special subsetting if GP model
  if("prob" %in% colnames(chain)){
    use_cols <- c(which(colnames(chain) != "prob"), which(colnames(chain) == "prob")[1])
    chain1 <- chain[,use_cols]
  }
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  ## Trace plots
  p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_line(aes(x=sampno,y=value,col=chain)) +
    facet_wrap(~name,scales="free_y")+
    scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),by=2000)) +
    export_theme
  
  ## Posterior density plots
  p_densities <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_density(aes(x=value,fill=chain),alpha=0.25) +
    facet_wrap(~name,scales="free") +
    export_theme
  
  ## Get predicted incidence trends
  predictions <- plot_prob_infection(chain_comb, 100, inc_func_use,
                                     times,
                                     obs_dat=dat_tmp_used_in_run)
  p1 <- predictions$plot
  
  ## Get fits to observed Ct distribution
  model_func <- create_posterior_func(parTab,dat_tmp_used_in_run,NULL,inc_func_use,"model")
  p2 <- plot_distribution_fits(chain_comb, dat_tmp_used_in_run, model_func,100)
  
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
  
  trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_inc <- ggplot(trajs_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + geom_line(data=tibble(t=1:349,y=seir_dynamics$incidence/population_n),aes(x=t,y=y),col="red")
  
  #p3 <- plot_posterior_density(chain, "beta",parTab, 0, 1000, real_data=TRUE)
  #p4 <- plot_posterior_density(chain, "overall_prob",parTab, 0, 1000, real_data=TRUE)
  
  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/posteriors/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/grs/"),recursive = TRUE)
  
  #dir.create(paste0(plot_wd,"/exp/"),recursive = TRUE)
  #dir.create(paste0(plot_wd,"/overall_prob/"),recursive = TRUE)
  
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",runname_use,"_predictions.png"),p1,width=7,height=4)
  ggsave(paste0(plot_wd,"/distributions/",runname_use,"_distributions.png"),p2,
         width=(7/5) * length(unique(obs_dat_tmp$t)),height=6)
  ggsave(paste0(plot_wd,"/posteriors/",runname_use,"_densities.png"),p_densities,width=7,height=4)
  ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  
}


