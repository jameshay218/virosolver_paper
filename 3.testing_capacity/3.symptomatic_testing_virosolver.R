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
library(doParallel)
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
chainwd <- paste0(HOME_WD, "/virosolver_paper/mcmc_chains/ReportSims/Symptom_rerun")
plot_wd <- paste0(HOME_WD, "/virosolver_paper/plots/ReportSims/Symptom_rerun")
data_wd <- paste0(HOME_WD, "/virosolver_paper/data/ReportSims/Symptom")
results_wd <- paste0(HOME_WD, "/virosolver_paper/results/ReportSims/Symptom_rerun/virosolver/")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(data_wd)) dir.create(data_wd,recursive = TRUE)
if(!file.exists(results_wd)) dir.create(results_wd,recursive = TRUE)


## Where to perform the simulations
setwd(main_wd)

n_clusters <- 3
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

########################################
## 2. Simulation settings
########################################
## Arguments for this run. Can generat automatically here
control_table <- expand.grid(pop_no=1:100,strategy=c("R1","R2","R3"),testing_day=as.character(c(7,8,9))) %>% arrange(pop_no, testing_day, strategy)
control_table$run_name <- 1:nrow(control_table)
control_table <- control_table[c(5L, 6L, 9L, 14L, 15L, 18L, 23L, 24L, 27L, 32L, 33L, 36L, 41L, 
                                 42L, 45L, 50L, 51L, 54L, 59L, 60L, 63L, 68L, 69L, 72L, 77L, 78L, 
                                 81L, 86L, 87L, 90L, 95L, 96L, 99L, 104L, 105L, 108L, 113L, 114L, 
                                 117L, 122L, 123L, 126L, 131L, 132L, 135L, 140L, 141L, 144L, 149L, 
                                 150L, 153L, 158L, 159L, 167L, 168L, 171L, 176L, 177L, 180L, 185L, 
                                 186L, 189L, 194L, 195L, 198L, 203L, 204L, 207L, 212L, 213L, 216L, 
                                 221L, 222L, 225L, 230L, 231L, 234L, 239L, 240L, 243L, 248L, 249L, 
                                 252L, 257L, 258L, 261L, 266L, 267L, 270L, 275L, 276L, 279L, 284L, 
                                 285L, 288L, 293L, 294L, 297L, 302L, 303L, 306L, 311L, 312L, 315L, 
                                 320L, 321L, 324L, 329L, 330L, 333L, 338L, 339L, 342L, 347L, 348L, 
                                 351L, 356L, 357L, 360L, 365L, 366L, 368L, 369L, 374L, 375L, 378L, 
                                 383L, 384L, 387L, 392L, 393L, 396L, 401L, 402L, 405L, 410L, 411L, 
                                 414L, 419L, 420L, 423L, 428L, 429L, 432L, 437L, 438L, 441L, 446L, 
                                 447L, 450L, 455L, 456L, 459L, 464L, 465L, 468L, 473L, 474L, 477L, 
                                 482L, 483L, 486L, 491L, 492L, 495L, 500L, 501L, 504L, 509L, 510L, 
                                 513L, 518L, 519L, 522L, 527L, 528L, 531L, 536L, 537L, 540L, 545L, 
                                 546L, 549L, 554L, 555L, 558L, 563L, 564L, 567L, 572L, 573L, 576L, 
                                 581L, 582L, 585L, 590L, 591L, 594L, 600L, 608L, 609L, 612L, 617L, 
                                 618L, 621L, 626L, 627L, 630L, 635L, 636L, 639L, 644L, 645L, 648L, 
                                 653L, 654L, 657L, 662L, 663L, 666L, 671L, 672L, 675L, 680L, 681L, 
                                 684L, 689L, 690L, 693L, 698L, 699L, 702L, 707L, 708L, 711L, 716L, 
                                 717L, 720L, 725L, 726L, 729L, 735L, 738L, 743L, 744L, 747L, 752L, 
                                 753L, 756L, 761L, 762L, 765L, 770L, 771L, 774L, 779L, 780L, 783L, 
                                 788L, 789L, 792L, 796L, 797L, 798L, 801L, 805L, 806L, 807L, 810L, 
                                 815L, 816L, 819L, 824L, 825L, 828L, 833L, 834L, 837L, 841L, 842L, 
                                 843L, 846L, 851L, 852L, 855L, 860L, 861L, 864L, 869L, 870L, 873L, 
                                 877L, 878L, 879L, 882L, 886L, 887L, 888L, 891L, 895L, 896L, 897L, 
                                 900L),]


n_temperatures <- 10
mcmcPars_ct <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                    "thin"=20,"adaptive_period"=30000,"save_block"=10000,"temperature" = seq(1,101,length.out=n_temperatures),
                    "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                    "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)
nchains <- 3
n_samp <- 1000 ## Number of posterior samples for plots

## Set Simulation Number and get sim settings
Sim <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#Sim <- 1
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
st <- proc.time()

## Symptom-Based Testing Schemes:
days <- 36
t_7 <- 1:(which.max(seir_dynamics$incidence)-2*7)
t_8 <- 1:(which.max(seir_dynamics$incidence)+2*7)
t_9 <- 1:(which.max(seir_dynamics$incidence)+6*7)

lastdays <- c(max(t_7),max(t_8),max(t_9))
names(lastdays) <- c("t_7","t_8","t_9")
p_A <- rep(.1, length.out=days)
p_B <- seq(from=.1,to=.2, length.out=days)
p_C <- exp(seq(from=log(.1),to=log(.25),length.out=days))
p_D <- seq(from=.1,to=.01, length.out=days)
p_E <- 0.2-exp(seq(from=log(.1),to=log(.195),length.out=days))
p_R <- c(0.003,0.003,0.003)
days_R1 <- c(0)
days_R2 <- c(7,0)
days_R3 <- c(14,7,0)

########################################
## 6. Fit virosolver model
########################################
Cases_Full <- NULL
Cts_Full <- NULL
Ests_Full <- NULL
R_Ests_Full <- NULL
True_SEIR_sim <- NULL

testdays <- get(paste0("t_",Days))
lastday <- max(testdays)

True_SEIR_sim <- bind_cols(True_SEIR %>% filter(t==lastday),Sim=Sim, TestDay=Days, TestProbs="T")

## Sampling for random sample with viral loads calculation for comparison:
tvp_R <- tibble(t=lastday-c(14,7,0), prob=p_R)
LL_R <- simulate_reporting(complete_linelist,
                           timevarying_prob=tvp_R,
                           solve_times=testdays,
                           symptomatic=FALSE)
SVL <- simulate_viral_loads_wrapper(LL_R$sampled_individuals,
                                    kinetics_pars=pars)
Cts_Full <- bind_cols(Sim=Sim, TestDay=Days, TestProbs="R", SVL %>% dplyr::select(sampled_time, ct_obs))


run_name <- paste0("RepScen_Pos_",Sim,"_",Days,Probs)
savewd <- paste0(chainwd, "/", run_name)
if(!file.exists(savewd)) dir.create(savewd,recursive = TRUE)

print(paste0("Analyzing Scenario: ",Days,Probs))
obs_times <- lastday-get(paste0("days_",Probs))
obs_dat <- Cts_Full %>% dplyr::select(sampled_time,ct_obs) %>% 
  dplyr::filter(sampled_time %in% obs_times) %>% 
  rename(t=sampled_time, ct=ct_obs)
obs_times <- lastday

## For SEIR model, use pars/partab_seir_model.csv
## For exp model all Cts, use pars/partab_exp_model.csv
## For exp model only +ve Cts, use pars/partab_exp_pos_model.csv
## For the GP model, use pars/partab_gp_model.csv
parTab <- read.csv(paste0(main_wd,"pars/massachusetts/partab_seir_model.csv"),stringsAsFactors=FALSE)
# parTab <- read.csv("pars/partab_exp_model.csv",stringsAsFactors=FALSE)
# parTab <- read.csv("pars/partab_exp_pos_model.csv",stringsAsFactors=FALSE)
#parTab <- read.csv("pars/partab_gp_model.csv",stringsAsFactors=FALSE)
pars <- parTab$values
names(pars) <- parTab$names
## Means for priors
means <- parTab$values
names(means) <- parTab$names

########################################
## Choose models, priors and MCMC control parameters
########################################
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!

## CHOOSE INCIDENCE FUNCTION
is_gp_run <- FALSE
## For exponential growth model
# inc_func_use <- exponential_growth_model
## For SEIR model
inc_func_use <- solveSEIRModel_rlsoda_wrapper
## For GP model
#inc_func_use <- gaussian_process_model
#is_gp_run <- TRUE

## CHOOSE PRIOR FUNCTION
## ** SEE code/priors.R **
## For SEIR model:
prior_func_use <- prior_func_hinge_seir
## For exp model:
# prior_func_use <- prior_func_hinge_exp
## For GP model:
# prior_func_use <- prior_func_hinge_gp

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
## For each observation time, fit nchain chains
rerun_mcmc <- FALSE
if(rerun_mcmc){
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
    res <- foreach(j=1:nchains,.packages = c("extraDistr","tidyverse","patchwork")) %dopar% {
      devtools::load_all(paste0(HOME_WD,"/virosolver"))
      devtools::load_all(paste0(HOME_WD,"/lazymcmc"))
    #for(j in 1:nchains){
      # ## Create posterior function. Need parTab (parameter control table), the data being used (obs_dat), the prior function, the incidence function being used (eg. the SEIR model), whether positives only (use_pos), and any other arguments to the inc and prior functions with ...
      # f <- create_posterior_func(parTab, obs_dat, PRIOR_FUNC=prior_func_use, INCIDENCE_FUNC=inc_func_use,solve_ver = "likelihood",use_pos=FALSE)
      # ## Test if returns a finite likelihood with the default parameters
      # f(parTab$values)
      ## Get random starting values
      if(n_temperatures > 1){
        startTab <- rep(list(parTab),n_temperatures)
        for(k in 1:length(startTab)){
          startTab[[k]] <- generate_viable_start_pars(parTab,obs_dat_tmp,
                                                      create_posterior_func,
                                                      inc_func_use,
                                                      prior_func_use,
                                                      t_dist=NULL,
                                                      use_pos=TRUE)
          startTab[[k]][startTab[[k]]$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_tmp$t)
        }
      }
      covMat <- diag(nrow(startTab[[1]]))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
      mvrPars <- rep(list(mvrPars), n_temperatures)
      
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
  runname_use <- paste0(run_name,"_",obs_times[i])
  chainwd_tmp <- paste0(savewd,"/",runname_use)
  chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=TRUE,PTchain=TRUE)$chain
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
  
  ## Get smoothed growth rates and Rt estimates
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  Rts <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    samp_pars <- get_index_pars(chain_comb, samps[ii])
    trajs[ii,] <- pmax(inc_func_use(samp_pars,times), 0.0000001)
    Rts[ii,] <- (1-cumsum(trajs[ii,]))*unname(samp_pars["R0"])
  }
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.5,0.025,0.25,0.75,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  trajs_quants <- as.data.frame(t(apply(trajs, 2, function(x) quantile(x,c(0.5,0.025,0.25,0.75,0.975)))))
  trajs_quants$t <- 0:(nrow(trajs_quants)-1)
  Rts_quants <- as.data.frame(t(apply(Rts, 2, function(x) quantile(x,c(0.5,0.025,0.25,0.75,0.975)))))
  Rts_quants$t <- 0:(nrow(Rts_quants)-1)
  colnames(trajs1_quants) <- c("median","lower_95","lower_50","upper_50","upper_95","t")
  colnames(Rts_quants) <- c("median","lower_95","lower_50","upper_50","upper_95","t")
  colnames(trajs_quants) <- c("median","lower_95","lower_50","upper_50","upper_95","t")
  R_Ests_Full <- bind_rows(R_Ests_Full,
                           bind_cols(Sim=Sim, TestDay=Days, TestProbs=Probs, variable="infections", trajs_quants),
                           bind_cols(Sim=Sim, TestDay=Days, TestProbs=Probs, variable="growth_rate", trajs1_quants),
                           bind_cols(Sim=Sim, TestDay=Days, TestProbs=Probs, variable="R", Rts_quants))
  
  ## Growth rate plot
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower_95,ymax=upper_95),alpha=0.25) +
    geom_line(aes(x=t,y=median)) +
    coord_cartesian(ylim=c(-0.5,0.5))
  
  if(!file.exists(paste0(plot_wd,"/traces/"))) dir.create(paste0(plot_wd,"/traces/"),recursive=TRUE)
  if(!file.exists(paste0(plot_wd,"/predictions/"))) dir.create(paste0(plot_wd,"/predictions/"),recursive=TRUE)
  if(!file.exists(paste0(plot_wd,"/distributions/"))) dir.create(paste0(plot_wd,"/distributions/"),recursive=TRUE)
  if(!file.exists(paste0(plot_wd,"/posteriors/"))) dir.create(paste0(plot_wd,"/posteriors/"),recursive=TRUE)
  if(!file.exists(paste0(plot_wd,"/grs/"))) dir.create(paste0(plot_wd,"/grs/"),recursive=TRUE)
  
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",runname_use,"_predictions.png"),p1,width=7,height=4)
  ggsave(paste0(plot_wd,"/distributions/",runname_use,"_distributions.png"),p2,
         width=(7/5) * length(unique(obs_dat_tmp$t)),height=6)
  ggsave(paste0(plot_wd,"/posteriors/",runname_use,"_densities.png"),p_densities,width=7,height=4)
  ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
  print(proc.time() - st)
}

save(list=c("Cts_Full","R_Ests_Full","True_SEIR_sim"), file=paste0(results_wd,"/Sim_viro_pop",pop_no,"_",Days,"", Probs,".Rda"))
