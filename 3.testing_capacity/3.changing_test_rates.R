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
library(doParallel)
library(lazymcmc)

## If fitting case counts:
# library(EpiNow2)

HOME_WD <- "~"
#HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))
#devtools::load_all(paste0(HOME_WD,"/lazymcmc"))


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
chainwd <- paste0(HOME_WD, "/virosolver_paper/mcmc_chains/ReportSims/ChangingTests")
plot_wd <- paste0(HOME_WD, "/virosolver_paper/plots/ReportSims/ChangingTests")
data_wd <- paste0(HOME_WD, "/virosolver_paper/data/ReportSims/Symptom")
results_wd <- paste0(HOME_WD, "/virosolver_paper/results/ReportSims/ChangingTests/virosolver/")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(data_wd)) dir.create(data_wd,recursive = TRUE)
if(!file.exists(results_wd)) dir.create(results_wd,recursive = TRUE)


dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/distributions/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/posteriors/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/grs/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/exp/"),recursive = TRUE)
dir.create(paste0(plot_wd,"/overall_prob/"),recursive = TRUE)

setwd(main_wd)

########################################
## 2. Simulation settings
########################################
## Arguments for this run. Can generat automatically here
control_table <- expand.grid(pop_no=1:100,strategy=c("E","F","G","H"),testing_day=c("4","5","6")) %>% arrange(pop_no, testing_day, strategy)
control_table$run_name <- 1:nrow(control_table)

#n_temperatures <- 2
#mcmcPars_ct <- list("iterations"=30000,"popt"=0.44,"opt_freq"=1000,
#                    "thin"=50,"adaptive_period"=30000,"save_block"=10000,"temperature" = seq(1,101,length.out=n_temperatures),
#                    "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
#                    "adaptiveLeeway" = 0.2, "max_total_iterations" = 30000)

mcmcPars_ct <- c("iterations"=50000,"popt"=0.44,"opt_freq"=2000,
                 "thin"=10,"adaptive_period"=30000,"save_block"=1000)

nchains <- 3
n_samp <- 1000 ## Number of posterior samples for plots

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


complete_linelist <- read.csv(paste0(data_wd,"/CompleteLineList_",pop_no,".csv"))
infection_days <- read.csv(paste0(data_wd,"/TrueGRs_",pop_no,".csv"))

True_SEIR <- tibble(t=rep(times,3), variable=rep(c("infections","growth_rate","R"), each=length(times)),
                    median=c(seir_dynamics$seir_outputs$I, c(NA,seir_dynamics$growth_rates[seir_dynamics$growth_rates[,"ver"]=="daily","GR"]), 
                             seir_dynamics$seir_outputs$Rt))
trueVals <- True_SEIR %>% filter(variable=="infections")
trueVals$logN <- ifelse(trueVals$median>0,log(trueVals$median),NA)
trueVals$GR14 <- NA
trueVals$GR28 <- NA
trueVals$GR35 <- NA
d14 <- 1:14
d28 <- 1:28
d35 <- 1:35
for (i in 2:length(trueVals$logN)) {
  if (i-13 > 1) {
    logNs <- trueVals$logN[(i-13):i]
    res <- lm(logNs~d14)
    trueVals$GR14[i] <- unname(coef(res)["d14"])
    if (i-27 > 1) {
      logNs <- trueVals$logN[(i-27):i]
      res <- lm(logNs~d28)
      trueVals$GR28[i] <- unname(coef(res)["d28"])
      if (i-34 > 1) {
        logNs <- trueVals$logN[(i-34):i]
        res <- lm(logNs~d35)
        trueVals$GR35[i] <- unname(coef(res)["d35"])
      }
    }
  }
}
trueVals <- pivot_longer(trueVals %>% dplyr::select(t,GR14,GR28,GR35) %>% filter(!is.na(GR14)),
                         cols=c(GR14,GR28,GR35), names_to="variable", values_to="median")
True_SEIR <- bind_rows(True_SEIR, trueVals %>% filter(!is.na(median)))


########################################
## 5. Simulate observation process and analyze symptom-based case report data
########################################
days <- 36
t_4 <- seq(from=98-35, by=1, length.out=days)
t_5 <- seq(from=126-35, by=1, length.out=days)
t_6 <- seq(from=154-35, by=1, length.out=days)

p_E <- rep(0.002, days)/5
p_F <- seq(from=0.001, to=0.003, length.out=days)/5
p_G <- exp(seq(from=log(0.0005), to=log(0.0041), length.out=days))/5
p_H <- seq(from=0.003, to=0.001, length.out=days)/5

# for (Sim in 1:SimNo) {
print(paste0("Starting Simulation Number: ",Sim))

Cts_Full <- NULL
Pos_Counts <- NULL
Res_Comb <- NULL

## This is where top level for loop was previously
testdays <- get(paste0("t_",Days))
lastday <- max(testdays)

TrueInf <- infection_days %>% dplyr::select(infection_time, N) %>% 
  rename(sampled_time=infection_time) %>% filter(sampled_time %in% testdays)
TrueInf$Tests <- NA
TrueInf$Positive <- TrueInf$N/20
Pos_Counts <- bind_rows(Pos_Counts, bind_cols(Sim=Sim, TestDay=Days, Probs=Probs, TestProbs="True Incidence", TrueInf %>% dplyr::select(!N)))

Res_Comb <- bind_rows(Res_Comb,
                      bind_cols(Sim=Sim, TestDay=Days, TestProbs="SEIR", 
                                True_SEIR %>% filter(t==lastday)),
                      bind_cols(Sim=Sim, TestDay=Days, Probs=Probs, TestProbs="Cases",
                                infection_days %>% filter(infection_time==lastday) %>% 
                                  dplyr::select(infection_time, N, GR, GR14, GR28, GR35, Rt) %>% 
                                  rename(t=infection_time, infections=N, R=Rt) %>% 
                                  pivot_longer(cols=c(infections, GR, GR14, GR28, GR35, R),
                                               names_to="variable", values_to="median")))

## This is where second level for loop was previously
print(paste0("Starting Scenario ",Days,Probs))

########################################
## 2. Generate Reported Data
########################################
tvp <- tibble(t=testdays, prob=get(paste0("p_",Probs)))
LL <- simulate_reporting(complete_linelist,
                         timevarying_prob=tvp,
                         solve_times=testdays,
                         symptomatic=FALSE)
SVL <- simulate_viral_loads_wrapper(LL$sampled_individuals,
                                    kinetics_pars=pars)
Cts_Sim <- SVL %>% dplyr::select(sampled_time, ct_obs)
Cts_Full <- bind_rows(Cts_Full, bind_cols(Sim=Sim, TestDay=Days, Probs=Probs, TestProbs=Probs, Cts_Sim))

Cts_Sim <- Cts_Sim %>% group_by(sampled_time) %>% summarize(N=sum(ct_obs<pars["intercept"]))
Cts_Sim$logN <- log(Cts_Sim$N)
Cts_Sim$GR <- c(NA,Cts_Sim$logN[2:length(Cts_Sim$logN)]/Cts_Sim$logN[1:(length(Cts_Sim$logN)-1)])
Cts_Sim$GR14 <- NA
Cts_Sim$GR28 <- NA
Cts_Sim$GR35 <- NA
d14 <- 1:14
d28 <- 1:28
d35 <- 1:35
for (i in 2:length(Cts_Sim$logN)) {
  if (i-13 > 1) {
    logNs <- Cts_Sim$logN[(i-13):i]
    res <- lm(logNs~d14)
    Cts_Sim$GR14[i] <- unname(coef(res)["d14"])
    if (i-27 > 1) {
      logNs <- Cts_Sim$logN[(i-27):i]
      res <- lm(logNs~d28)
      Cts_Sim$GR28[i] <- unname(coef(res)["d28"])
      if (i-34 > 1) {
        logNs <- Cts_Sim$logN[(i-34):i]
        res <- lm(logNs~d35)
        Cts_Sim$GR35[i] <- unname(coef(res)["d35"])
      }
    }
  }
}
Cts_Sim <- Cts_Sim %>% filter(sampled_time==lastday) %>% 
  dplyr::select(sampled_time, N, GR, GR14, GR28, GR35) %>% 
  rename(t=sampled_time, infections=N) %>% 
  pivot_longer(cols=c(infections, GR, GR14, GR28, GR35),
               names_to="variable", values_to="median")
Cts_Sim$Sim <- Sim
Cts_Sim$TestDay <- Days
Cts_Sim$TestProbs <- Probs
Cts_Sim$variable <- paste0(Cts_Sim$variable,"_ObsCases")

Res_summ <- SVL %>% dplyr::select(sampled_time, ct_obs) %>% group_by(sampled_time) %>% 
  summarize(Tests=n(), Positive=sum(ct_obs<pars["intercept"]))
Pos_Counts <- bind_rows(Pos_Counts, bind_cols(Sim=Sim, TestDay=Days, Probs=Probs, TestProbs=Probs, Res_summ))

for (Type in c("Pos","All")) {
  run_name <- paste0("Rates_",Sim,"_",Days,Probs,"_",Type)
  savewd <- paste0(chainwd, "/", run_name)
  setwd(main_wd)
  
  if(!file.exists(savewd)) dir.create(savewd,recursive = TRUE)
  
  obs_dat <- SVL %>% dplyr::select(sampled_time, ct_obs) %>%
    rename(t=sampled_time, ct=ct_obs)
  obs_times <- lastday
  
  ## For SEIR model, use pars/partab_seir_model.csv
  ## For exp model all Cts, use pars/partab_exp_model.csv
  ## For exp model only +ve Cts, use pars/partab_exp_pos_model.csv
  ## For the GP model, use pars/partab_gp_model.csv
  # parTab <- read.csv("pars/partab_seir_model.csv",stringsAsFactors=FALSE)
  if (Type=="Pos") {
    parTab <- read.csv("pars/massachusetts/partab_exp_pos_model.csv",stringsAsFactors=FALSE)
    ## Only detectable or all Cts?
    only_pos <- TRUE
  } else if (Type=="All") {
    parTab <- read.csv("pars/massachusetts/partab_exp_model.csv",stringsAsFactors=FALSE)
    ## Only detectable or all Cts?
    only_pos <- FALSE
  }
  #parTab <- read.csv("pars/partab_gp_model.csv",stringsAsFactors=FALSE)
  
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
  
  ########################################
  ## 3. Choose models, priors and MCMC control parameters
  ########################################
  ## CHOOSE INCIDENCE FUNCTION
  is_gp_run <- FALSE
  ## For exponential growth model
  inc_func_use <- exponential_growth_model
  ## For SEIR model
  # inc_func_use <- solveSEIRModel_rlsoda_wrapper
  ## For GP model
  #inc_func_use <- gaussian_process_model
  #is_gp_run <- TRUE
  
  ## CHOOSE PRIOR FUNCTION
  ## ** SEE code/priors.R **
  ## For SEIR model:
  # prior_func_use <- prior_func_hinge_seir
  ## For exp model:
  prior_func_use <- prior_func_hinge_exp
  ## For GP model:
  # prior_func_use <- prior_func_hinge_gp
  
  ## Independent cross sections or combined?
  cumu_data <- TRUE
  
  ## How old can infections be? Try 35 for single cross section models.
  ## For full SEIR model, set this to NA and the code will use the start time of the epidemic
  max_age <- 35
  
  ########################################
  ## 4. Fit model to Ct values
  ########################################
    ## For each observation time, fit nchain chains
    # res <- foreach(i=seq_along(obs_times),.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
    # devtools::load_all("../virosolver")
  for(i in 1:length(obs_times)) {
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
      # ## Create posterior function. Need parTab (parameter control table), the data being used (obs_dat), the prior function, the incidence function being used (eg. the SEIR model), whether positives only (use_pos), and any other arguments to the inc and prior functions with ...
      # f <- create_posterior_func(parTab, obs_dat, PRIOR_FUNC=prior_func_use, INCIDENCE_FUNC=inc_func_use,solve_ver = "likelihood",use_pos=FALSE)
      # ## Test if returns a finite likelihood with the default parameters
      # f(parTab$values)
      ## Get random starting values
      startTab <- generate_viable_start_pars(parTab,obs_dat_tmp,
                                             create_posterior_func,
                                             inc_func_use,
                                             prior_func_use,
                                             t_dist=NULL,
                                             use_pos=only_pos)
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
  
  
  ########################################
  ## 5. Load in MCMC chains
  ########################################
  res <- NULL
  for(i in seq_along(obs_times)){
    runname_use <- paste0(run_name,"_",obs_times[i])
    chainwd_tmp <- paste0(savewd,"/",runname_use)
    chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
    chain <- as.data.frame(chain)
    res[[i]] <- chain
  }
  
  ########################################
  ## 6. All plots from this run
  ########################################
  ## Plot useful summary plots
  for(i in seq_along(obs_times)){
    obs_time <- obs_times[i]
    runname_use <- paste0(run_name,"_",obs_time)
    
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
    times_obs <- 0:max(dat_tmp_used_in_run$t)
    
    ## This is for the GP version
    mat <- matrix(rep(times_obs, each=length(times_obs)),ncol=length(times_obs))
    t_dist <- abs(apply(mat, 2, function(x) x-times_obs))
    
    if(is_gp_run){
      parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times_obs),])
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
    # predictions <- plot_prob_infection(chain_comb, 100, inc_func_use,
                                       # times_obs,
                                       # obs_dat=dat_tmp_used_in_run)
    # p1 <- predictions$plot
    
    ## Get fits to observed Ct distribution
    # model_func <- create_posterior_func(parTab,dat_tmp_used_in_run,NULL,inc_func_use,"model")
    # p2 <- plot_distribution_fits(chain_comb, dat_tmp_used_in_run, model_func,100)
    
    ## Get estimated growth rates
    samps <- sample(unique(chain_comb$sampno),n_samp)
    trajs <- matrix(0, nrow=n_samp,ncol=length(times_obs))
    for(ii in seq_along(samps)){
      trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times_obs), 0.0000001)
    }
    trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
    trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
    trajs1_quants <- as.data.frame(trajs1_quants)
    trajs1_quants$t <- 1:nrow(trajs1_quants)
    colnames(trajs1_quants) <- c("lower","median","upper","t")
    
    t1q <- trajs1_quants %>% filter(t==1)
    t1q$t <- lastday
    Res_Comb <- bind_rows(Res_Comb,
                          Cts_Sim,
                          bind_cols(Sim=Sim, TestDay=Days, Probs=Probs,
                                    TestProbs=Probs, 
                                    variable=paste0("Est_",Type), t1q))
    
    
    
    ## Growth rate plot
    p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
      geom_line(aes(x=t,y=median)) + 
      coord_cartesian(ylim=c(-0.5,0.5))
    
    ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
    ggsave(paste0(plot_wd,"/posteriors/",runname_use,"_densities.png"),p_densities,width=7,height=4)
    ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  }
}
write_csv(Res_Comb, paste0(results_wd,"/Res_Comb_Sim_popno",pop_no,"_",Days,"",Probs,".csv"))
write_csv(Cts_Full, paste0(results_wd,"/Cts_Full_Sim_popno",pop_no,"_",Days,"",Probs,".csv"))
write_csv(Pos_Counts, paste0(results_wd,"/Pos_Counts_Sim_popno",pop_no,"_",Days,"",Probs,".csv"))
