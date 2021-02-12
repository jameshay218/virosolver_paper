## James Hay 16 October 2020
##  - This script runs the analyses to investigate transmission dynamics in 4 Massachussetts care homes
##  - First, we look at the prevalence of positive samples at 3 sampling rounds
##  - Next, we fit a modified SEIR model (SEEIRR) to these point prevalence estimates
##  - Assuming that these fits demonstrate "ground truth", we then fit the Ct model to each individual cross section
##  - Finally, we fit the overall SEIR model to all 3 cross sections comparing Ct use and pos/neg only

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

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc <- TRUE

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
#main_wd <- setwd("~/Documents/Harvard/Research/Inectious Diseases/COVID-19/PCR Cts/Code & Results/ct_inference")
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- "~/Documents/GitHub/virosolver_paper/mcmc_chains/1.nursing_home" ## Where to save MCMC chains
setwd(main_wd)
source("code/plot_funcs.R")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)


nchains <- 1
n_clusters <- 8
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)


## Change these MCMC control parameters to run longer/shorter chains
mcmcPars1 <- c("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)
mcmcPars2 <- c("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=50000,"save_block"=1000)


## Note that data are grouped by week of collection
dat_wide <- read_csv(file = "data/nursing_home_data.csv")
dat_wide$collection_date <- as.Date(dat_wide$collection_date,format="%m/%d/%y")
dat_wide$mean_week_date <- as.Date(dat_wide$mean_week_date,format="%m/%d/%y")

##########################################
# 1. Look at raw data for prevalence
##########################################
## Look at prevalence over time
prev_time <- dat_wide %>%
  filter(is_nursinghome == 1) %>%
  filter(result %in% c("POS","NEG")) %>%
  group_by(location, result, week) %>%
  tally() %>%
  pivot_wider(values_from=n, names_from=result,values_fill=list(n=0)) %>%
  group_by(location, week) %>%
  mutate(prev=POS/(POS+NEG),
         lower_confint=prop.test(POS, POS+NEG)$conf.int[1],
         upper_confint=prop.test(POS,POS+NEG)$conf.int[2])

p_prev <- prev_time %>% ggplot() +
  geom_point(aes(x=week,y=prev),size=0.25) +
  geom_errorbar(aes(x=week,ymin=lower_confint,ymax=upper_confint),size=0.25,width=0.1) +
  facet_wrap(~location,ncol=2) +
  ylab("Prevalence (95% binomial confidence intervals)") +
  xlab("Week") +
  ggtitle("Prevalence over time by location") +
  export_theme +
  theme(legend.position=c(0.8,0.1),
        axis.text.x=element_text(angle=0,hjust=0.5),
        panel.grid.major = element_line(color="grey70",size=0.1)) +
  scale_y_continuous(limits=c(0,1))
p_prev


##########################################
# 2. Fit SEIR model
##########################################
parTab <- read.csv("pars/partab_seeirr_combined.csv",stringsAsFactors=FALSE)

prior_func_seeirr <- function(pars){
  sds <- c(2,0.5,0.5,2,1)*2
  names(pars) <- parTab$names
  p1 <- dnorm(pars["R0_1"], 0, 100,log=TRUE)
  p1 <- dnorm(pars["R0_2"], 0, 100,log=TRUE)
  p1 <- dnorm(pars["R0_3"], 0, 100,log=TRUE)
  p1 <- dnorm(pars["R0_4"], 0, 100,log=TRUE)
  p2 <- dnorm(pars["latent"],2,sds[2],log=TRUE)
  p3 <- dnorm(pars["incubation"],2,sds[3],log=TRUE)
  p4 <- dnorm(pars["infectious"],10,sds[4],log=TRUE)
  p5 <- dnorm(pars["recovery"],5, sds[5], log=TRUE)
  return(sum(p1, p2,p3,p4,p5))
}

test_institutions <- c("NH 1","NH 2", "NH 3","NH 4")

## Duration of epidemic in days
times <- 0:365
t<-times
min_date <- as.Date("2020-01-01")
model_func <- create_post_func_seeirr(parTab, data=NULL, ts=times,PRIOR_FUNC=prior_func_seeirr,ver="model")
dat <- model_func(parTab$values)

## Make Jan 1st day 0
nh_prev1 <- dat_wide %>%
  filter(result %in% c("POS","NEG")) %>%
  filter(location %in% test_institutions) %>%
  group_by(location, result, week, mean_week_date) %>%
  tally() %>%
  group_by(location, result, mean_week_date, week) %>%
  summarize(n=sum(n)) %>%
  pivot_wider(values_from=n, names_from=result,values_fill=list(n=0)) %>%
  group_by(location, mean_week_date) %>%
  mutate(prev=POS/(POS+NEG),
         lower_confint=prop.test(POS, POS+NEG)$conf.int[1],
         upper_confint=prop.test(POS,POS+NEG)$conf.int[2]) %>%
  ungroup()  %>%
  rename(date=mean_week_date) %>%
  mutate(date = as.numeric(date - min_date))


#posterior <- create_post_func_seeirr(parTab, data=nh_prev1 %>% filter(location == "NH 1"),
#                                     ts=times,PRIOR_FUNC=prior_func_seeirr,ver="likelihood")
#dat <- posterior(parTab$values)


posterior <- create_post_func_seeirr_combined(parTab, data=nh_prev1,
                                     ts=times,PRIOR_FUNC=prior_func_seeirr,ver="likelihood")
dat <- posterior(parTab$values)


#res <- foreach(i=seq_along(test_institutions),.packages = c("lazymcmc","rethinking","extraDistr","ggthemes","tidyverse","deSolve")) %dopar% {
#
#  devtools::load_all("~/Documents/GitHub/virosolver")
#  data_tmp <- nh_prev1 %>% filter(location == test_institutions[i])
startTab <- generate_start_tab(parTab)
output <- run_MCMC(parTab=startTab, data=nh_prev1, mcmcPars=mcmcPars1,
                     filename=paste0(chainwd,"/",test_institutions[i]),
                     CREATE_POSTERIOR_FUNC = create_post_func_seeirr_combined, mvrPars = NULL,
                   PRIOR_FUNC=prior_func_seeirr, ts=times)
chain <- read.csv(output$file)
best_pars <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,0.5,w=0.8)
## Start from best location of previous chain
parTab$values <- best_pars
## Run second chain
output <- run_MCMC(parTab=parTab, data=nh_prev1, mcmcPars=mcmcPars1,
                   filename=paste0(chainwd,"/",test_institutions[i]),
                   CREATE_POSTERIOR_FUNC = create_post_func_seeirr_combined, mvrPars = mvrPars,
                   PRIOR_FUNC=prior_func_seeirr, ts=times)
chain <- read.csv(output$file)


quants_all <- NULL
random_traj_all <- NULL
t0_all <- NULL
best_traj_all <- NULL
n_samp <- 100
times1 <- 1:max(nh_prev1$date)

for(i in seq_along(test_institutions)){
  chain <- res[[i]]
  samps <- sample(unique(chain$sampno), n_samp)
  data_tmp <- nh_prev1 %>% filter(location == test_institutions[i])

  dat <- matrix(0, nrow=n_samp, ncol=length(times1))

  model_func <- create_post_func_seeirr(parTab, data=data_tmp,
                                        ts=times1,PRIOR_FUNC=prior_func_seeirr,ver="model")

  best_pars <- get_best_pars(chain)
  best_traj <- tibble(time=times1, y=model_func(best_pars)[1:length(times1)]) %>%
    mutate(time=as.Date(time + min_date)) %>%
    mutate(location=test_institutions[i])

  for(j in seq_along(samps)){
    pars <- get_index_pars(chain, samps[j])
    dat[j,] <- model_func(pars)[seq_along(times1)]
  }
  quants <- t(apply(dat, 2, function(x) quantile(x, c(0.025,0.5,0.975))))
  quants <- as.data.frame(quants)
  colnames(quants) <- c("lower","median","upper")
  quants$time <- times1
  quants <- quants %>%
    mutate(time=as.Date(time + min_date)) %>%
    mutate(location=test_institutions[i])

  ## Get the random trajectories
  dat <- as_tibble(dat)
  colnames(dat) <- times1
  dat$samp <- 1:nrow(dat)
  random_traj <- dat %>%
    pivot_longer(-samp) %>%
    mutate(time=as.numeric(name)) %>%
    filter(samp %in% sample(unique(dat$samp), 25)) %>%
    mutate(time = as.Date(time + min_date)) %>%
    mutate(location=test_institutions[i])

  ## Get point-range for t0
  t0_pointrange <- as.Date(quantile(chain$t0,c(0.025,0.5,0.975)) + min_date)
  t0_pointrange <- data.frame(t(data.frame(t0_pointrange)))
  colnames(t0_pointrange) <- c("lower","median","upper")
  t0_pointrange$lower <- as.Date(t0_pointrange$lower)
  t0_pointrange$median <- as.Date(t0_pointrange$median)
  t0_pointrange$upper <- as.Date(t0_pointrange$upper)
  t0_pointrange$location <- test_institutions[i]

  quants_all <- bind_rows(quants_all, quants)
  random_traj_all <- bind_rows(random_traj_all, random_traj)
  best_traj_all <- bind_rows(best_traj_all, best_traj)
  t0_all <- bind_rows(t0_all, t0_pointrange)
}

ggplot() +
  geom_line(data=random_traj_all,aes(x=time,y=value,group=samp),size=0.05) +
  geom_ribbon(data=quants_all,
              aes(ymin=lower,ymax=upper,x=time),alpha=0.1,fill="blue") +
  geom_point(data=nh_prev1 %>%
               mutate(date=as.Date(date + min_date)), aes(x=date, y=prev),
             col="darkgreen") +
  geom_errorbar(data=nh_prev1 %>%
                  mutate(date=as.Date(date + min_date)),
                aes(x=date,ymin=lower_confint, ymax=upper_confint),
                width=2, col="darkgreen") +
  geom_line(data=best_traj_all, aes(x=time,y=y)) +
  geom_pointrange(data=t0_all,aes(xmin=lower,xmax=upper,x=median),y=0.1, col="red",size=0.25) +
  geom_label(x=t0_pointrange$median[1] + 5, y=0.2,label="Seed time estimate",size=2) +
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
  scale_x_date(breaks="7 days") +
  ylab("PCR detectable prevalence") +
  xlab("Date") +
  export_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~location,nrow=1)

