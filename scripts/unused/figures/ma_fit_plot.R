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
set.seed(1)
n_samp <- 1000

runname_seir <- "real_ma_seir_pt"
runname_exp <- "real_ma_exp"

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/5.real_ma_single_timepoint/")
plot_wd <- paste0("~/Documents/GitHub/virosolver_paper/plots/5.real_ma_single_timepoint/")
setwd(main_wd)

## MCMC parameters for Ct model fits
mcmcPars_ct_seir <- c("iterations"=100000,"popt"=0.44,"opt_freq"=5000,
                      "thin"=100,"adaptive_period"=50000,"save_block"=1000)
mcmcPars_ct_exp <- c("iterations"=200000,"popt"=0.44,"opt_freq"=5000,
                     "thin"=100,"adaptive_period"=100000,"save_block"=1000)

## Code for plotting
source("code/plot_funcs.R")

########################################
## 2. Model parameters and simulation settings
########################################
## GP model parameters for fitting
parTab_exp <- read.csv("pars/partab_exp_pos_model.csv")
parTab_seir <- read.csv("pars/partab_seir_model.csv")

########################################
## 3. Read in MA data
########################################
#obs_dat_all <- read_csv("~/Documents/GitHub/ct_inference_preprint/data/BWH_COVID_Cts_deid_20200403-20200831.csv") %>%
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

res_seir <- NULL
res_exp <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  
  ## Read in SEIR chains
  chainwd_tmp <- paste0(chainwd,"/",runname_seir,"/",timepoint)
  chain_seir <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,
                                 #mcmcPars_ct_seir["adaptive_period"],
                                 30000,
                                 multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir <- as.data.frame(chain_seir)
  chain_seir$sampno <- 1:nrow(chain_seir)
  res_seir[[i]] <- chain_seir
  
  ## Read in exp chains
  chainwd_tmp <- paste0(chainwd,"/",runname_exp,"/",timepoint)
  chain_exp <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct_exp["adaptive_period"],multi=TRUE,chainNo=TRUE)$chain
  chain_exp <- as.data.frame(chain_exp)
  chain_exp$sampno <- 1:nrow(chain_exp)
  res_exp[[i]] <- chain_exp
}

grs_daily <- NULL
grs_average <- NULL

prop_grow_daily <- NULL
prop_grow_average <- NULL

ps_daily <- NULL
ps_average <- NULL

beta_max <- 0.5
overall_circle_top <- data.frame(y=c(0,seq(0,beta_max,by=0.001)),x=c(0,sqrt(beta_max^2-seq(0,beta_max,by=0.001)^2)))
overall_circle_bot <- data.frame(y=c(0,seq(-beta_max,0,by=0.001)),x=c(0,sqrt(beta_max^2-seq(-beta_max,0,by=0.001)^2)))
beta_max1 <- beta_max*0.7
small_overall_circle_top <- data.frame(y=c(0,seq(0,beta_max1,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(0,beta_max1,by=0.001)^2)))
small_overall_circle_bot <- data.frame(y=c(0,seq(-beta_max1,0,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(-beta_max1,0,by=0.001)^2)))


for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use_seir <- runname_seir
  runname_use_exp <- runname_exp
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  obs_dat_seir <- obs_dat_exp <- obs_dat_use
  
  ## Observation times
  obs_dat_exp <- obs_dat_exp %>% mutate(t = t - min(t), t = t + 35)
  
  ages_exp <- 1:max(obs_dat_exp$t)
  times_exp <- 0:max(obs_dat_exp$t)
  ages_seir <- 1:max(obs_dat_seir$t)
  times_seir <- 0:max(obs_dat_seir$t)
  
  chain_seir <- res_seir[[i]]
  chain_comb_seir <- chain_seir
  chain_comb_seir$sampno <- 1:nrow(chain_comb_seir)
  chain1_seir <- chain_seir
  chain_comb_seir <- chain_comb_seir[,colnames(chain_comb_seir) != "chain"]
  
  ## Get daily growth rate
  samps <- sample(unique(chain_comb_seir$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  for(j in seq_along(samps)){
    trajs[j,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_seir, samps[j]),times_seir),0.0000001)
  }
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  tmp_beta <- trajs1[,ncol(trajs1)]
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2)))))
  
  dats_all <- NULL
  
  p <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill=AAAS_palette["red1"],alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill=AAAS_palette["green1"],alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=y,x=sqrt(beta_max^2 - y^2),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all[[i]] <- tmp
    p <- p + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.03,fill=AAAS_palette["blue1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5)
  ps_daily[[i]] <- p + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=med_segment,x=0,xend=sqrt((beta_max^2) - med_segment^2)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col=AAAS_palette["grey1"]) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    geom_hline(yintercept = 0,linetype="dashed",size=0.25,color="grey40") +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
  chain_exp <- res_exp[[i]]
  chain_comb_exp <- chain_exp
  chain_comb_exp$sampno <- 1:nrow(chain_comb_exp)
  chain1_exp <- chain_exp
  chain_comb_exp <- chain_comb_exp[,colnames(chain_comb_exp) != "chain"]
  
  
  grs_daily[[i]] <- trajs1_quants[nrow(trajs1_quants),]
  grs_average[[i]] <- c(quantile(chain_comb_exp$beta, c(0.025,0.25,0.5,0.75,0.975)),timepoint)
  
  prop_grow_daily[[i]] <- sum(trajs1[,ncol(trajs1)] > 0)/nrow(trajs1)
  prop_grow_average[[i]] <- sum(chain_comb_exp$beta > 0)/nrow(chain_comb_exp)
}

