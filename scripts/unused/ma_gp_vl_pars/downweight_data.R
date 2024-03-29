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
n_samp <- 1000
runname <- "downweight"
run_version <- "gp" ##gp, seir or exp##
rerun_mcmc <- TRUE

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/TESTS/",runname)
plot_wd <- paste0("~/Documents/GitHub/virosolver_paper/plots/TESTS/",runname)
setwd(main_wd)

## Manage MCMC runs and parallel runs
nchains <- 2
n_clusters <- 2
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=700000,"popt"=0.44,"opt_freq"=2000,
                 "thin"=100,"adaptive_period"=300000,"save_block"=1000)

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
inc_func_use <- gaussian_process_model
prior_func_use <- prior_func_hinge_gp

## GP model parameters for fitting
parTab <- read.csv("pars/partab_gp_model.csv")
#parTab <- read.csv("pars/partab_gp_model_start.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("nu","rho"), "fixed"] <- 0
parTab[parTab$names =="overall_prob", c("values","fixed")] <- c(1,1)
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

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
  filter(coll_date > "2020-04-15" & coll_date < "2020-11-10") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct)# %>%
  #filter(!(t %in% c(63,77)))

obs_dat_all <- obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15" & coll_date < "2020-11-10") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  arrange(first_day) %>%
  rename(date = first_day, ct=panther_Ct) #%>%
  #filter(!(date %in% as.Date(c("2020-05-10","2020-05-24"))))

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
  scale_x_date(limits=as.Date(c("2020-01-01","2020-12-01")),breaks="1 month") +
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
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat1$t)
#parTab[!(parTab$names %in% c("prob","rho","nu","obs_sd")),"fixed"] <- 1

parTab[parTab$names == "prob","values"] <- reverse_gp_model(rep(0.001,length(times)),pars, times)

f <- create_posterior_func(parTab, obs_dat1, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=TRUE,
                           t_dist=t_dist)
pars <- parTab$values
names(pars) <- parTab$names
prior_func_use(pars,t_dist)
f(parTab$values)

## Run for each chain
chains <- NULL
if(rerun_mcmc){
  res <- foreach(j=1:nchains,.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
      devtools::load_all("~/Documents/GitHub/virosolver")
      
    ## Get random starting values
    startTab <- generate_viable_start_pars(parTab,obs_dat1,
                                           create_posterior_func,
                                           inc_func_use,
                                           prior_func_use,
                                           t_dist=t_dist,
                                           use_pos=TRUE)
    
    startTab <- parTab
    covMat <- diag(nrow(startTab))
    mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
    
    output <- run_MCMC(parTab=startTab,
                       data=obs_dat1,
                       INCIDENCE_FUNC=inc_func_use,
                       PRIOR_FUNC = prior_func_use,
                       solve_likelihood=TRUE,
                       mcmcPars=mcmcPars_ct,
                       filename=paste0(chainwd,"/",runname,"_chainno_",j),
                       CREATE_POSTERIOR_FUNC=create_posterior_func,
                       mvrPars=NULL,
                       OPT_TUNING=0.2,
                       use_pos=TRUE,
                       t_dist=t_dist)
    
    ## Read in chain and remove burn in period
    chain <- read.csv(output$file)
    chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
    chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
    chains[[j]] <- chain
    chain <- do.call("bind_rows",chains)
  } 
}

res <- NULL

chains_diag <- lazymcmc::load_mcmc_chains(chainwd, parTab,TRUE,1,mcmcPars_ct["adaptive_period"],
                                     multi=FALSE,chainNo=FALSE,PTchain = FALSE)

gelman_diagnostics(chains_diag$list,threshold=1.1)
chains <- lazymcmc::load_mcmc_chains(chainwd, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                      multi=FALSE,chainNo=TRUE,PTchain = FALSE)
chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
chain_comb <- chain
chain_comb$sampno <- 1:nrow(chain_comb)

## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
trajs <- matrix(0, nrow=n_samp,ncol=length(times))
vl_trajs <- matrix(0, nrow=n_samp, ncol=length(ages))
detect_trajs <- matrix(0, nrow=n_samp, ncol=length(ages))

for(ii in seq_along(samps)){
  tmp_pars <- get_index_pars(chain_comb, samps[ii])
  vl <- viral_load_func(tmp_pars,ages,TRUE)
  vl_trajs[ii,]  <- extraDistr::rgumbel(length(vl),vl,tmp_pars["obs_sd"])
  detect_trajs[ii,] <- sapply(ages, function(a) prop_detectable(a, tmp_pars,vl))
  trajs[ii,] <- pmax(smooth.spline(inc_func_use(tmp_pars,times))$y,0.0000001)
  #trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
}

trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs1_quants <- as.data.frame(trajs1_quants)
trajs1_quants$t <- 1:nrow(trajs1_quants)
colnames(trajs1_quants) <- c("lower","median","upper","t")
trajs1_quants <- trajs1_quants %>% left_join(date_key)
## Growth rate plot
p_gr <- ggplot(trajs1_quants) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=date,y=median),col=AAAS_palette["blue1"]) + 
  export_theme +
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Daily growth rate") +
  xlab("Date") +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5))

detect_trajs1 <- t(apply(detect_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
detect_trajs1 <- as.data.frame(detect_trajs1)
colnames(detect_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
detect_trajs1$t <- ages

p_detect <- ggplot(detect_trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme

vl_trajs1 <- t(apply(vl_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
vl_trajs1 <- as.data.frame(vl_trajs1)
colnames(vl_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
vl_trajs1$t <- ages

p_vl <- ggplot(vl_trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Days since infection") +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  coord_cartesian(ylim=c(40,0)) +
  ylab("Ct value") +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme


trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
trajs_quants <- as.data.frame(trajs_quants)
trajs_quants$t <- 1:nrow(trajs_quants)
colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")

obs_dat_all <- obs_dat_all %>% left_join(obs_dat_all %>% group_by(date) %>% tally())

p_dat <- ggplot(obs_dat_all) + 
  geom_violin(aes(x=date,group=date,y=ct,fill=n),scale="width",
              alpha=0.5,color="black",size=0.1,
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=date, y=ct,group=date),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data=obs_dat_all %>% group_by(date) %>% summarize(median_ct=median(ct)),
              aes(x=date,y=median_ct),col=AAAS_palette["blue1"],se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(42, 5),expand=c(0,0)) +
  scale_fill_gradient(low=AAAS_palette["blue1"],high=AAAS_palette["red1"]) +
  geom_hline(yintercept=40,linetype="dashed") +
  export_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = c(0.05,0.5)) +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) +
  xlab("Date of sample") +
  ylab("Ct value") +
  labs(tag="B")
p_dat

trajs_melted <- reshape2::melt(trajs)
colnames(trajs_melted) <- c("samp","t","inc")
trajs_melted$t <- times[trajs_melted$t]
trajs_quants <- trajs_quants %>% left_join(date_key)
## Growth rate plot
p_inc <- ggplot(trajs_quants) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=date,y=median),col=AAAS_palette["blue1"]) + 
  #geom_line(data=trajs_melted[trajs_melted$samp %in% 1:10,], aes(x=t,y=inc,group=samp)) +
  #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
  #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
  export_theme +
  ylab("Relative probability of infection") +
  xlab("Date") +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) +
  #coord_cartesian(ylim=c(-0.0001,0.005)) +
  scale_y_continuous(expand=c(0,0))  +
  labs(tag="C")


nyt_dat <- read_csv("~/Documents/GitHub/covid-19-data/us-states.csv")
ma_dat <- nyt_dat %>% filter(state == "Massachusetts") %>%
  mutate(new_cases = lead(cases,1) - cases)


p_nyt <- ggplot(ma_dat) + 
  geom_bar(aes(x=date,y=new_cases),stat="identity") +
  export_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5000)) +
  xlab("Date") +
  ylab("Number of cases") +
  labs(tag="A")


