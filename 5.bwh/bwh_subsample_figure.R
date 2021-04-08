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

#HOME_WD <- "~/"
HOME_WD <- "~/Documents/GitHub"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Arguments for this run
index <- 1994
set.seed(index)
n_samp <- 1000
runnames <- c("ma_gp_subsample5perc","ma_gp_subsample10perc","ma_gp_subsample25perc",
              "ma_gp_subsample25","ma_gp_subsample50","ma_gp")
run_version <- "gp" ##gp, seir or exp##

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- paste0(HOME_WD,"/virosolver_paper/")
setwd(main_wd)

## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=500000,"popt"=0.44,"opt_freq"=2000,
                 "thin"=350,"adaptive_period"=200000,"save_block"=100)

## Code for plotting
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors_tighter.R")

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

########################################
## 2. Model parameters and simulation settings
########################################
inc_func_use <- gaussian_process_model
prior_func_use <- prior_func_hinge_gp

## GP model parameters for fitting
parTab <- read.csv(paste0(main_wd,"/pars/massachusetts/partab_gp_model.csv"))
#parTab <- read.csv("pars/partab_gp_model_start.csv",stringsAsFactors=FALSE)
parTab[parTab$names %in% c("nu","rho"), "values"] <- c(1.5,0.03)
#parTab[parTab$names %in% c("nu","rho"), "fixed"] <- 0
parTab[parTab$names %in% c("nu","rho"), "fixed"] <- 1
parTab[parTab$names =="overall_prob", c("values","fixed")] <- c(1,1)
parTab[parTab$names =="level_switch", "fixed"] <- 1
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
obs_dat_all <- read_csv(paste0(main_wd,"/data/panther_Ct_20200403-20201110.csv")) %>% rename(panther_Ct=ORF1ab_Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n()) %>%
  filter(panther_Ct < 40)

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

## Subsample data to use at most subsamp_n samples per timepoint
obs_dat_old <- obs_dat1

f <- create_posterior_func(parTab, obs_dat1, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=TRUE,
                           t_dist=t_dist)
f(parTab$values)
store_trajs <- NULL
index <- 1
for(runname in runnames){
  setwd(main_wd)
  chainwd <- paste0(HOME_WD,"/virosolver_paper/mcmc_chains/4.real_ma_ct/",runname)
  plot_wd <- paste0(HOME_WD,"/virosolver_paper/plots/4.real_ma_ct/",runname)
  
  chains <- lazymcmc::load_mcmc_chains(chainwd, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                       multi=FALSE,chainNo=TRUE,PTchain = FALSE)
  chain <- as.data.frame(chains$chain)
  chain$sampno <- 1:nrow(chain)
  chain_comb <- chain[chain$chain == 1,]
  chain_comb$sampno <- 1:nrow(chain_comb)
  
  model_func <- create_posterior_func(parTab,obs_dat1,NULL,inc_func_use,"model")
  p_dist_fits <- plot_distribution_fits_bwh(chain_comb, obs_dat1, model_func,1000,TRUE,date_key)
  
  ## Get smoothed growth rates
  test_ages <- 0:50
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  
  for(ii in seq_along(samps)){
    tmp_pars <- get_index_pars(chain_comb, samps[ii])
    trajs[ii,] <- pmax(smooth.spline(inc_func_use(tmp_pars,times))$y,0.0000001)
  }
  
  trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")
  
  trajs_melted <- reshape2::melt(trajs)
  colnames(trajs_melted) <- c("samp","t","inc")
  trajs_melted$t <- times[trajs_melted$t]
  trajs_quants <- trajs_quants %>% left_join(date_key)
  trajs_quants$runname <- runname
  store_trajs[[index]] <- trajs_quants
  index <- index + 1
}
samp_size_5perc <- obs_dat1 %>% sample_frac(0.05) %>% tally()
samp_size_10perc <- obs_dat1 %>% sample_frac(0.10) %>% tally()
samp_size_25perc <- obs_dat1 %>% sample_frac(0.25) %>% tally()
samp_size_25 <-  obs_dat1 %>% 
  group_by(t) %>% 
  tally() %>% 
  mutate(n_use = pmin(n, 25)) %>% 
  left_join(obs_dat1) %>%
  group_by(t) %>%
  sample_n(n_use) %>%
  dplyr::select(t, ct)%>%
  ungroup() %>%
  tally()
samp_size_50 <-  obs_dat1 %>% 
  group_by(t) %>% 
  tally() %>% 
  mutate(n_use = pmin(n, 50)) %>% 
  left_join(obs_dat1) %>%
  group_by(t) %>%
  sample_n(n_use) %>%
  dplyr::select(t, ct) %>%
  ungroup() %>%
  tally()

runname_key <- c("ma_gp_subsample5perc"=paste0("Subsample 5% (N=",samp_size_5perc,")"),
                 "ma_gp_subsample10perc"=paste0("Subsample 10% (N=",samp_size_10perc,")"),
                 "ma_gp_subsample25perc"=paste0("Subsample 25% (N=",samp_size_25perc,")"),
                 "ma_gp_subsample25"=paste0("Subsample 25 Cts per week (N=",samp_size_25,")"),
                 "ma_gp_subsample50"=paste0("Subsample 50 Cts per week (N=",samp_size_50,")"),
                 "ma_gp"=paste0("All data (N=",nrow(obs_dat1),")"))


trajs_quants <- do.call("bind_rows",store_trajs)
trajs_quants$runname <- runname_key[trajs_quants$runname]
trajs_quants$runname <- factor(trajs_quants$runname, levels=runname_key[c(6,5,4,3,2,1)])
## Growth rate plot
p_inc <- ggplot(trajs_quants) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=date,y=median),col=AAAS_palette["blue1"]) + 
  geom_hline(yintercept=0,size=0.75) +
  export_theme +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1)) +
  ylab("Relative probability of infection") +
  xlab("Date") +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-11-25")),breaks="1 month",expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.02),breaks=seq(0,0.02,by=0.005)) +
  facet_wrap(~runname,ncol=2,dir="v")
ggsave("figures/supplement/bwh_subsampled.png",height=6,width=8,dpi=300,units="in")
ggsave("figures/supplement/bwh_subsampled.pdf",height=6,width=8)
