## IMPORTANT NOTE: the single cross section SEIR models are indexed from 1st February, so t[0] = 1st February there.
## The exponential and GP models are indexed from 35 days prior to the first sampling time. So all SEIR model times are 36
## days after the exp and GP models. HOWEVER, we index obs_dat based on the exponential and GP models. This is why
## you will see some times as 35, and others as 71 referencing the same day.

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

## Set flag to write plots to disk or not
save_plots <- FALSE

## CHANGE TO MAIN WD
main_wd <- "~/Documents/GitHub/virosolver_paper/"
setwd(main_wd)
source("code/plot_funcs.R")

export_theme <- export_theme + theme(plot.margin=unit(c(0,0,0,0),units="cm"),plot.tag=element_text(size=10,face="bold"))


## Arguments for this run
set.seed(1)
n_samp <- 1000
## How dates translate to observation times, GP model
obs_times <- c(35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119, 126, 
               133, 140, 147, 154, 161, 168, 175, 182, 189, 196, 203, 210, 217, 
               224, 231, 238, 245)
## In SEIR cross section fits, times are shifted differently
seir_timeshift <- 36

## Run names for the various MCMC chains
runname_seir <- "ma_seir_constrained" ## Constrained prior on seed time
runname_seir_prior <- "ma_seir_prior_constrained"
runname_seir_prior_broad <- "ma_seir_rerun_prior"
runname_seir_alt <- "ma_seir_rerun" ## Less constrained prior on seed time
runname_exp <- "ma_exp"

## Where the MCMC chains are stored
chainwd <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/5.real_ma_single_timepoint/")
chainwd_gp <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/4.real_ma_ct/ma_gp/")

## MCMC parameters for Ct model fits
mcmcPars_ct_seir <- c("adaptive_period"=30000)
mcmcPars_ct_exp <- c("adaptive_period"=30000)

parTab_exp <- read.csv("pars/massachusetts/partab_exp_pos_model.csv")
parTab_seir <- read.csv("pars/massachusetts/partab_seir_model.csv")

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

########################################
## 2. Read in MCMC chain
########################################
res_seir <- NULL
res_exp <- NULL
res_seir_alt <- NULL
res_seir_prior <- NULL
res_seir_prior_broad <- NULL

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  ## Read in SEIR chains constrained t0
  chainwd_tmp <- paste0(chainwd,"/",runname_seir,"/",timepoint + seir_timeshift)
  chain_seir <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1, mcmcPars_ct_seir["adaptive_period"],
                                 multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir <- as.data.frame(chain_seir)
  chain_seir$sampno <- 1:nrow(chain_seir)
  res_seir[[i]] <- chain_seir
  
  ## Read in SEIR chains prior
  chainwd_tmp <- paste0(chainwd,"/",runname_seir_prior,"/",timepoint+seir_timeshift)
  chain_seir_prior <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1, mcmcPars_ct_seir["adaptive_period"],
                                 multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir_prior <- as.data.frame(chain_seir_prior)
  chain_seir_prior$sampno <- 1:nrow(chain_seir_prior)
  res_seir_prior[[i]] <- chain_seir_prior
  
  ## Read in SEIR chains from loose t0
  chainwd_tmp <- paste0(chainwd,"/",runname_seir_alt,"/",timepoint+seir_timeshift)
  chain_seir_alt <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1, mcmcPars_ct_seir["adaptive_period"],
                                     multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir_alt <- as.data.frame(chain_seir_alt)
  chain_seir_alt$sampno <- 1:nrow(chain_seir_alt)
  res_seir_alt[[i]] <- chain_seir_alt
  
  ## Read in SEIR prior chains from loose t0
  chainwd_tmp <- paste0(chainwd,"/",runname_seir_prior_broad,"/",timepoint+seir_timeshift)
  chain_seir_prior_broad <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1, mcmcPars_ct_seir["adaptive_period"],
                                       multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir_prior_broad <- as.data.frame(chain_seir_prior_broad)
  chain_seir_prior_broad$sampno <- 1:nrow(chain_seir_prior_broad)
  res_seir_prior_broad[[i]] <- chain_seir_prior_broad
  
  ## Read in exp chains
  chainwd_tmp <- paste0(chainwd,"/",runname_exp,"/",timepoint)
  chain_exp <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct_exp["adaptive_period"],multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_exp <- as.data.frame(chain_exp)
  chain_exp$sampno <- 1:nrow(chain_exp)
  res_exp[[i]] <- chain_exp
}

########################################
## 3. MA incidence plot and Rt
########################################
## NYT data
#nyt_dat <- read_csv("data/us-states.csv")
nyt_dat <- read_csv("~/Documents/GitHub/covid-19-data/us-states.csv")
nyt_dat <- read_csv("~/Documents/GitHub/covid-19-data/us-counties.csv")
nyt_dat <- nyt_dat %>% 
  filter(state=="Massachusetts",
         county == "Suffolk") %>%
  group_by(state) %>%
  mutate(new_cases=cases-lag(cases, 1))
## Rt fit
#estimates <- readRDS("~/Documents/GitHub/ct_dynamics_preprint/results/ma_rt_fit.RData")
estimates <- readRDS("results/ma_county_rt_fit.RData")
rt_dat <- estimates$estimates$summarised %>% filter(variable %in% c("growth_rate") & type=="estimate")

## Split based on above or below Rt = 1
rt_dat_top <- rt_dat %>%
  mutate(bottom=pmax(0,bottom),
         top=pmax(0,top))
rt_dat_bot <- rt_dat %>%
  mutate(bottom=pmin(0,bottom),
         top=pmin(0,top))
## Same thing for median line
rt_median <- rt_dat %>% mutate(is_grow=ifelse(median>0,"Growing","Declining"))
index <- 1
rt_median$index <- index
for(i in 2:nrow(rt_median)){
  if(rt_median$is_grow[i] != rt_median$is_grow[i-1]) {
    index <- index + 1
  }
  rt_median$index[i] <- index
}
coeff <- 2500
yshift <- 500
p1 <-  ggplot(nyt_dat)+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.5,fill=AAAS_palette["grey1"]) +
  #geom_ribbon(data=dat_infections,aes(x=date,ymin=bottom,ymax=top),
  #            fill=AAAS_palette["grey1"],alpha=0.5) +
  #geom_line(data=dat_infections,aes(x=date,y=mean),col=AAAS_palette["grey1"]) +
  geom_hline(yintercept=yshift,linetype="dashed",col=AAAS_palette["grey1"]) +
  geom_ribbon(data=rt_dat_top,aes(x=date,ymin=bottom*coeff + yshift,ymax=top*coeff + yshift),
              fill="darkorange",alpha=0.25) +
  geom_ribbon(data=rt_dat_bot,aes(x=date,ymin=bottom*coeff + yshift,ymax=top*coeff + yshift),
              fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=rt_median,aes(x=date,y=median*coeff + yshift,col=is_grow,group=index)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,1000),
                     sec.axis=sec_axis(~.*1/coeff - yshift/coeff, name="Growth rate")) +
  scale_color_manual(values=c("Growing"="darkorange","Declining"=as.character(AAAS_palette["green1"])))+
  scale_x_date(limits=as.Date(c("2020-04-01", "2020-11-15"), "%Y-%m-%d"), breaks="7 days",
               expand=c(0,0)) +
  export_theme+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New cases") +
  labs(tag="A")




########################################
## 4. BWH Ct data
########################################
obs_dat_all <- read_csv("data/panther_Ct_20200403-20201110.csv") %>% rename(panther_Ct=ORF1ab_Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n()) %>%
  filter(panther_Ct < 40)

obs_dat1 <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  #mutate(earliest_date=as.Date("2020-02-01")) %>%
  #mutate(first_day = as.numeric(first_day - earliest_date)) %>% ## Assume 1st February is start date
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

p_dat <- ggplot(obs_dat_all) + 
  geom_violin(aes(x=date,group=date,y=ct),scale="width",fill=AAAS_palette["grey1"],
              alpha=0.5,color="black",size=0.1,
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=date, y=ct,group=date),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data=obs_dat_all %>% group_by(date) %>% summarize(median_ct=median(ct)),
              aes(x=date,y=median_ct),col=AAAS_palette["blue1"],se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(42, 10),expand=c(0,0)) +
  geom_hline(yintercept=40,linetype="dashed") +
  export_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_date(limits=as.Date(c("2020-04-01","2020-11-15")),breaks="1 month",expand=c(0,0)) +
  xlab("Date of sample") +
  ylab("Ct value") +
  labs(tag="C")

########################################
## 5. Rt scatterplot
########################################
## Plot Ct data
bwh_data <- read_csv("data/panther_Ct_20200403-20201110.csv") %>% rename(panther_Ct=ORF1ab_Ct) %>%
  mutate(platform="Panther",first_pos=1) %>%
  mutate(id=1:n())

bwh_data_use <-  bwh_data %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  rename(date=coll_date) %>%
  left_join(epi_calendar)

bwh_data_week <- bwh_data_use %>% 
  group_by(date) %>% 
  summarize(median_ct=median(panther_Ct),
            skew_ct=moments::skewness(panther_Ct),
            n=n())

combined_dat <- bwh_data_week %>% 
  full_join(estimates$estimates$summarised %>% as_tibble() %>% 
              filter(variable == "R") %>%
              mutate(date = date)
  ) %>%
  arrange(date) %>% 
  rename(mean_rt=mean) %>%
  dplyr::select(date, median_ct,skew_ct, n, mean_rt) %>%
  filter(n >= 10)

combined_dat1 <- combined_dat %>% 
  pivot_longer(-c(date,n)) %>%
  left_join(epi_calendar) %>%
  group_by(first_day, name) %>% 
  summarize(mean=mean(value)) %>%
  pivot_wider(values_from=mean,names_from=name)
combined_dat1 <- combined_dat

median_cts <- smooth(combined_dat$median_ct)
skew_cts <- combined_dat$skew_ct
Rt <- combined_dat$mean_rt
ccf(median_cts, Rt,lag.max = 25)
ccf(skew_cts, Rt,lag.max = 25)

p_rt <- ggplot(combined_dat1 %>% rename(Rt=mean_rt)) +
  geom_point(aes(x=skew_ct,y=median_ct,col=Rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1)+
  scale_y_continuous(trans="reverse") +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  export_theme + 
  theme(legend.position=c(0.2,0.7)) + 
  labs(tag="B")

########################################
## 6. GP dials
########################################
grs_daily <- NULL
grs_average <- NULL
prop_grow_daily <- NULL
prop_grow_average <- NULL
ps_daily <- NULL
ps_average <- NULL

trajs_dat_all <- NULL

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
  obs_dat_seir <- obs_dat_seir %>% mutate(t=t+seir_timeshift)
  
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
  
  chain_exp <- res_exp[[i]]
  chain_comb_exp <- chain_exp
  chain_comb_exp$sampno <- 1:nrow(chain_comb_exp)
  chain1_exp <- chain_exp
  chain_comb_exp <- chain_comb_exp[,colnames(chain_comb_exp) != "chain"]
  
  
  ## Get daily growth rate
  samps <- sample(unique(chain_comb_seir$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  t0s <- NULL
  for(k in seq_along(samps)){
    trajs[k,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_seir, samps[k]),times_seir),0)
    t0s[k] <- get_index_pars(chain_comb_seir, samps[k])["t0"]
  }
  #trajs <- t(apply(trajs,1,function(x) x/sum(x)))
  trajs_dat <- data.frame(trajs)
  colnames(trajs_dat) <- times_seir
  trajs_dat$samp <- 1:nrow(trajs_dat)
  trajs_dat$t0 <- t0s
  trajs_dat <- reshape2::melt(trajs_dat,id.vars=c("samp","t0"))
  
  colnames(trajs_dat) <- c("samp","t0","t","inc")
  trajs_dat$t <- as.numeric(trajs_dat$t)
  trajs_dat$obs_time <- timepoint
  
  trajs_dat_all[[i]] <- trajs_dat
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  tmp_beta <- trajs1[,ncol(trajs1)]
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), 
                                                        function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2),na.rm=TRUE))))
  
  dats_all <- NULL
  p_daily <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="darkorange",alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill=AAAS_palette["green1"],alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=beta_max*sin(pi*y),x=beta_max*cos(pi*y),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all[[i]] <- tmp
    p_daily <- p_daily + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.05,fill=AAAS_palette["blue1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5,na.rm=TRUE)
  ps_daily[[i]] <-  p_daily + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=beta_max*sin(pi*med_segment),x=0,xend=beta_max*cos(pi*med_segment)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="yellow",size=0.5) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
  tmp_beta <- chain_comb_exp$beta
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), 
                                                        function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2),na.rm=TRUE))))
  ## Average growth rates
  p_average <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="darkorange",alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill="#008B45FF",alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=beta_max*sin(pi*y),x=beta_max*cos(pi*y),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all[[i]] <- tmp
    p_average <- p_average + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.05,fill=AAAS_palette["teal1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5,na.rm=TRUE)
  ps_average[[i]] <-  p_average + 
    scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=beta_max*sin(pi*med_segment),x=0,xend=beta_max*cos(pi*med_segment)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="yellow",size=0.5) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
  chain_exp <- res_exp[[i]]
  chain_comb_exp <- chain_exp
  chain_comb_exp$sampno <- 1:nrow(chain_comb_exp)
  chain1_exp <- chain_exp
  chain_comb_exp <- chain_comb_exp[,colnames(chain_comb_exp) != "chain"]
  
  grs_daily[[i]] <- trajs1_quants[nrow(trajs1_quants),] %>% mutate(t=t-seir_timeshift)
  grs_average[[i]] <- c(quantile(chain_comb_exp$beta, c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE),timepoint)
  
  prop_grow_daily[[i]] <- sum(trajs1[,ncol(trajs1)] > 0)/nrow(trajs1)
  prop_grow_average[[i]] <- sum(chain_comb_exp$beta > 0)/nrow(chain_comb_exp)
}

trajs_dat_all_comb <- do.call("bind_rows", trajs_dat_all)


ps_all_daily <- ps_daily[[1]] + 
  #geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40") +
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="E")
ps_all_average <- ps_average[[1]] + 
  #geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40") +
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="E")
for(i in seq(4,length(obs_times),by=3)){
  ps_all_average <- ps_all_average | (ps_average[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
  ps_all_daily <- ps_all_daily | (ps_daily[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
}

beta_max2 <- beta_max*1.2
text_circle_top <- data.frame(y=c(0,beta_max2*cos(pi*seq(0,beta_max,by=beta_max/5))),
                              x=c(0,beta_max2*sin(pi*seq(0,beta_max,by=beta_max/5))),
                              text=c(NA,rev(seq(0,beta_max,by=0.1))))

text_circle_bot <- data.frame(y=c(0,-beta_max2*cos(pi*seq(0,beta_max,by=beta_max/5))),
                              x=c(0,beta_max2*sin(pi*seq(0,beta_max,by=beta_max/5))),
                              text=c(NA,seq(-beta_max,0,by=0.1)))

beta_lines <- data.frame(x=0,y=0,group=seq_along(seq(0,beta_max,by=beta_max/5)),
                         yend=beta_max2*cos(pi*seq(0,beta_max,by=beta_max/5)),
                         xend=beta_max2*sin(pi*seq(0,beta_max,by=beta_max/5)))


p_use_dial <-ps_daily[[10]] +
  scale_x_continuous(expand=c(0,0),limits=c(-0.1,beta_max+0.3)) +
  scale_y_continuous(expand=c(0,0),limits=c(-beta_max-0.2,beta_max+0.2)) +
  geom_vline(xintercept=0,size=0.1)+
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=-yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_label(data=text_circle_top,aes(x=x,y=y,label=text),size=2) +
  geom_label(data=text_circle_bot,aes(x=x,y=y,label=text),size=2) +
  geom_text(data=data.frame(x=0.6,y=0.5,label="Growth"),aes(x=x,y=y,label=label),size=3) +
  geom_text(data=data.frame(x=0.6,y=-0.5,label="Decline"),aes(x=x,y=y,label=label),size=3) +
  theme_void() +
  theme(plot.margin=unit(c(0,-4,0,-4),unit="cm"),
        plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="D")

########################################
## 7. GP fit
########################################
chains <- lazymcmc::load_mcmc_chains(chainwd_gp, parTab,FALSE,1,mcmcPars_ct_exp["adaptive_period"],
                                     multi=FALSE,chainNo=TRUE,PTchain = FALSE)
chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
chain_comb <- chain
chain_comb$sampno <- 1:nrow(chain_comb)

times <- 0:max(obs_dat1$t)
## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
gp_trajs <- matrix(0, nrow=n_samp,ncol=length(times))
for(ii in seq_along(samps)){
  tmp <- pmax(smooth.spline(gaussian_process_model(get_index_pars(chain_comb, samps[ii]),times))$y,0)
  gp_trajs[ii,] <- tmp/sum(tmp)
  #trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
}

gp_trajs1 <- t(apply(gp_trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
gp_trajs1_quants <- t(apply(gp_trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975),na.rm=TRUE)))
gp_trajs1_quants <- as.data.frame(gp_trajs1_quants)
gp_trajs1_quants$t <- 1:nrow(gp_trajs1_quants)
colnames(gp_trajs1_quants) <- c("lower","median","upper","t")

## Growth rate plot
p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
  geom_line(aes(x=t,y=median)) + 
  coord_cartesian(ylim=c(-0.5,0.5))


gp_trajs_quants <- t(apply(gp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
gp_trajs_quants <- as.data.frame(gp_trajs_quants)
gp_trajs_quants$t <- 1:nrow(gp_trajs_quants)
colnames(gp_trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")

## Growth rate plot
p_inc <- ggplot(gp_trajs_quants %>% left_join(date_key)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["red1"]) + 
  geom_ribbon(aes(x=date,ymin=mid_lower,ymax=mid_upper),alpha=0.5,fill=AAAS_palette["red1"]) + 
  geom_line(aes(x=date,y=median),col=AAAS_palette["red1"]) + 
  #geom_line(data=tibble(t=times,y=inc_func_use(get_best_pars(chain_comb),times)),aes(x=t,y=y),col="green") +
  #geom_line(data=tibble(t=1:200,y=(seir_dynamics$incidence/population_n)[1:200]),aes(x=t,y=y),col="red") +
  export_theme +
  ylab("Scaled probability of infection") +
  xlab("Date") +
  scale_x_date(limits=as.Date(c("2020-04-01","2020-11-15")),breaks="1 month",expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.0001,0.025)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(tag="E")


########################################
## 8. Growth rate comparison
########################################
gp_trajs1_quants_joined <- left_join(gp_trajs1_quants, date_key)
gr_dat <- estimates$estimates$summarised %>% filter(variable == "growth_rate") %>% 
  filter(type=="estimate") %>%
  dplyr::select(date, mean,top,bottom) %>%
  rename(gr_mean=mean,
         gr_lower=bottom,
         gr_upper=top)

gr_dat_top <- gr_dat %>%
  mutate(gr_lower=pmax(0,gr_lower),
         gr_upper=pmax(0,gr_upper))
gr_dat_bot <- gr_dat %>%
  mutate(gr_lower=pmin(0,gr_lower),
         gr_upper=pmin(0,gr_upper))

gr_dat_median <- gr_dat %>% mutate(is_grow=ifelse(gr_mean>0,"R(t), growing","R(t), declining"))
index <- 1
gr_dat_median$index <- index
for(i in 2:nrow(gr_dat_median)){
  if(gr_dat_median$is_grow[i] != gr_dat_median$is_grow[i-1]) {
    index <- index + 1
  }
  gr_dat_median$index[i] <- index
}


trajs_dat_top <- gp_trajs1_quants_joined %>%
  mutate(lower=pmax(0,lower),
         upper=pmax(0,upper))
trajs_dat_bot <- gp_trajs1_quants_joined %>%
  mutate(lower=pmin(0,lower),
         upper=pmin(0,upper))

trajs_dat_median <- gp_trajs1_quants_joined %>% mutate(is_grow=ifelse(median>0,"Ct estimate","Ct estimate"))
index <- 1
trajs_dat_median$index <- index
for(i in 2:nrow(trajs_dat_median)){
  if(trajs_dat_median$is_grow[i] != trajs_dat_median$is_grow[i-1]) {
    index <- index + 1
  }
  trajs_dat_median$index[i] <- index
}

p_grs <- ggplot() + 
  geom_ribbon(data=gr_dat_top,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill="darkorange") +
  geom_ribbon(data=gr_dat_bot,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill=AAAS_palette["green1"]) +
  geom_line(data=gr_dat_median,aes(x=date,y=gr_mean,col=is_grow,group=index)) +
  geom_ribbon(data=trajs_dat_top,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_ribbon(data=trajs_dat_bot,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_line(data=trajs_dat_median,aes(x=date,y=median,col=is_grow,group=index)) +
  coord_cartesian(ylim=c(-0.125,0.25)) +
  scale_x_date(expand=c(0,0),limits=as.Date(c("2020-03-09", "2020-12-01"), "%Y-%m-%d"), 
               breaks=as.Date(c("2020-04-01","2020-07-01","2020-10-01"))) +
  geom_hline(yintercept=0,linetype="dashed",col=AAAS_palette["grey1"]) +
  scale_color_manual(values=c("R(t), growing"="darkorange",
                              "R(t), declining"=as.character(AAAS_palette["green1"]),
                              "Ct estimate"=as.character(AAAS_palette["purple1"])))+
  export_theme +
  theme(legend.position=c(0.5,0.9),
        legend.title = element_blank(),
        axis.text.x=element_text(size=6),
        legend.text=element_text(size=6),
        plot.margin=unit(c(0,0,0,0),units="cm")) +
  xlab("Date") +
  ylab("Growth rate") +
  labs(tag="F")

########################################
## 9. Pull together
########################################
if(FALSE){
  pdf("figures/Figure4.pdf",height=6,width=9)
  (((p1/p_dat/((ps_all_daily | plot_spacer()))/p_inc) + plot_layout(heights=c(4,3.5,0.5,4))) | (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
  
  png("figures/Figure4.png",height=6,width=9,units="in",res=300)
  (((p1/p_dat/((ps_all_daily | plot_spacer()))/p_inc) + plot_layout(heights=c(4,3.5,0.5,4))) | (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
}
if(save_plots){
  pdf("figures/Figure4.pdf",height=6,width=9)
  ((p1/p_dat/p_inc)| (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
  
  png("figures/Figure4.png",height=6,width=9,units="in",res=300)
  ((p1/p_dat/p_inc)| (((plot_spacer()|p_rt)+plot_layout(widths=c(1,50)))/p_use_dial/p_grs)) + plot_layout(widths=c(2,1))
  dev.off()
}


########################################
## 10. Alternative dials and example timepoint
########################################
grs_daily_alt <- NULL ## Loose prior on seed time
prop_grow_daily_alt <- NULL
ps_daily_alt <- NULL
ps_daily_prior <- NULL ## Prior on constrained seed time
ps_daily_prior_broad <- NULL ## Prior loose seed time
trajs_dat_all_alt <- NULL

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use_seir <- runname_seir_alt
  runname_use_prior <- runname_seir_prior
  runname_use_prior_broad <- runname_seir_prior_broad
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  obs_dat_seir <- obs_dat_exp <- obs_dat_use
  obs_dat_seir <- obs_dat_seir %>% mutate(t=t+seir_timeshift)
  
  ages_seir <- 1:max(obs_dat_seir$t)
  times_seir <- 0:max(obs_dat_seir$t)
  
  chain_seir <- res_seir_alt[[i]]
  chain_comb_seir <- chain_seir
  chain_comb_seir$sampno <- 1:nrow(chain_comb_seir)
  chain1_seir <- chain_seir
  chain_comb_seir <- chain_comb_seir[,colnames(chain_comb_seir) != "chain"]
  
  chain_prior <- res_seir_prior[[i]]
  chain_comb_prior <- chain_prior
  chain_comb_prior$sampno <- 1:nrow(chain_comb_prior)
  chain1_prior <- chain_prior
  chain_comb_prior <- chain_comb_prior[,colnames(chain_comb_prior) != "chain"]
  
  chain_prior_broad <- res_seir_prior_broad[[i]]
  chain_comb_prior_broad <- chain_prior_broad
  chain_comb_prior_broad$sampno <- 1:nrow(chain_comb_prior_broad)
  chain1_prior_broad <- chain_prior_broad
  chain_comb_prior_broad <- chain_comb_prior_broad[,colnames(chain_comb_prior_broad) != "chain"]
  
  ###################################################
  ## A. Get daily growth rate with alternative seed prior
  ###################################################
  samps <- sample(unique(chain_comb_seir$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  t0s <- NULL
  for(k in seq_along(samps)){
    trajs[k,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_seir, samps[k]),times_seir),0)
    t0s[k] <- get_index_pars(chain_comb_seir, samps[k])["t0"]
  }
  #trajs <- t(apply(trajs,1,function(x) x/sum(x)))
  trajs_dat <- data.frame(trajs)
  colnames(trajs_dat) <- times_seir
  trajs_dat$samp <- 1:nrow(trajs_dat)
  trajs_dat$t0 <- t0s
  trajs_dat <- reshape2::melt(trajs_dat,id.vars=c("samp","t0"))
  
  colnames(trajs_dat) <- c("samp","t0","t","inc")
  trajs_dat$t <- as.numeric(trajs_dat$t)
  trajs_dat$obs_time <- timepoint
  
  trajs_dat_all_alt[[i]] <- trajs_dat
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  tmp_beta <- trajs1[,ncol(trajs1)]
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), 
                                                        function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2),na.rm=TRUE))))
  
  dats_all_alt <- NULL
  p_daily <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="darkorange",alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill="#5F559BFF",alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=beta_max*sin(pi*y),x=beta_max*cos(pi*y),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all_alt[[i]] <- tmp
    p_daily <- p_daily + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.05,fill=AAAS_palette["red1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5,na.rm=TRUE)
  ps_daily_alt[[i]] <-  p_daily + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=beta_max*sin(pi*med_segment),x=0,xend=beta_max*cos(pi*med_segment)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="yellow",size=0.5) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
  grs_daily_alt[[i]] <- trajs1_quants[nrow(trajs1_quants),]
  prop_grow_daily_alt[[i]] <- sum(trajs1[,ncol(trajs1)] > 0)/nrow(trajs1)
  
  ###################################################
  ## B. Priors
  ###################################################
  samps <- sample(unique(chain_comb_prior$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  t0s <- NULL
  for(k in seq_along(samps)){
    trajs[k,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_prior, samps[k]),times_seir),0)
    t0s[k] <- get_index_pars(chain_comb_prior, samps[k])["t0"]
  }
  #trajs <- t(apply(trajs,1,function(x) x/sum(x)))
  trajs_dat <- data.frame(trajs)
  colnames(trajs_dat) <- times_seir
  trajs_dat$samp <- 1:nrow(trajs_dat)
  trajs_dat$t0 <- t0s
  trajs_dat <- reshape2::melt(trajs_dat,id.vars=c("samp","t0"))
  
  colnames(trajs_dat) <- c("samp","t0","t","inc")
  trajs_dat$t <- as.numeric(trajs_dat$t)
  trajs_dat$obs_time <- timepoint
  
  #trajs_dat_all_alt[[i]] <- trajs_dat
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  tmp_beta <- trajs1[,ncol(trajs1)]
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), 
                                                        function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2),na.rm=TRUE))))
  
  dats_all_alt <- NULL
  p_daily <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="darkorange",alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill=AAAS_palette["green1"],alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=beta_max*sin(pi*y),x=beta_max*cos(pi*y),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all_alt[[i]] <- tmp
    p_daily <- p_daily + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.05,fill=AAAS_palette["grey1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5,na.rm=TRUE)
  ps_daily_prior[[i]] <-  p_daily + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=beta_max*sin(pi*med_segment),x=0,xend=beta_max*cos(pi*med_segment)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="black",size=0.5) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
  
  ###################################################
  ## C. Prior broad
  ###################################################
  samps <- sample(unique(chain_comb_prior_broad$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  t0s <- NULL
  for(k in seq_along(samps)){
    trajs[k,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_prior_broad, samps[k]),times_seir),0)
    t0s[k] <- get_index_pars(chain_comb_prior_broad, samps[k])["t0"]
  }
  #trajs <- t(apply(trajs,1,function(x) x/sum(x)))
  trajs_dat <- data.frame(trajs)
  colnames(trajs_dat) <- times_seir
  trajs_dat$samp <- 1:nrow(trajs_dat)
  trajs_dat$t0 <- t0s
  trajs_dat <- reshape2::melt(trajs_dat,id.vars=c("samp","t0"))
  
  colnames(trajs_dat) <- c("samp","t0","t","inc")
  trajs_dat$t <- as.numeric(trajs_dat$t)
  trajs_dat$obs_time <- timepoint
  
  #trajs_dat_all_alt[[i]] <- trajs_dat
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  tmp_beta <- trajs1[,ncol(trajs1)]
  beta_quants <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), 
                                                        function(y) quantile(tmp_beta, c(0.5-y/2,0.5+y/2),na.rm=TRUE))))
  
  dats_all_alt <- NULL
  p_daily <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="darkorange",alpha=0.25) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill=AAAS_palette["green1"],alpha=0.25) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  
  for(j in 1:nrow(beta_quants)){
    alpha <- beta_quants[j,1]
    lower <- beta_quants[j,2]
    upper <- beta_quants[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=beta_max*sin(pi*y),x=beta_max*cos(pi*y),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all_alt[[i]] <- tmp
    p_daily <- p_daily + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.05,fill=AAAS_palette["grey1"])
  }
  
  med_segment <- quantile(tmp_beta, 0.5,na.rm=TRUE)
  ps_daily_prior_broad[[i]] <-  p_daily + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=beta_max*sin(pi*med_segment),x=0,xend=beta_max*cos(pi*med_segment)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="black",size=0.5) +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank(),plot.margin= unit(c(0,0,0,0),"in"))
  
}

trajs_dat_all_comb_alt <- do.call("bind_rows", trajs_dat_all_alt)

p2 <- ggplot(trajs_dat_all_comb_alt %>% filter(samp <= 1000, obs_time == 56)) + 
  geom_line(aes(x=t,y=inc,col=t0 <= (29+seir_timeshift),group=samp))

supp_ps_all_daily_prior <- ps_daily_prior[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="A")
supp_ps_all_daily <- ps_daily[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="B")
supp_ps_all_average <- ps_average[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="C")
supp_ps_all_daily_prior_broad <- ps_daily_prior_broad[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="D")
supp_ps_all_daily_alt <- ps_daily_alt[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="E")

for(i in seq(2,length(obs_times),by=1)){
  supp_ps_all_daily_prior <- supp_ps_all_daily_prior / (ps_daily_prior[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
  supp_ps_all_daily_prior_broad <- supp_ps_all_daily_prior_broad / (ps_daily_prior_broad[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
  supp_ps_all_daily_alt <- supp_ps_all_daily_alt / (ps_daily_alt[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
  supp_ps_all_daily <- supp_ps_all_daily / (ps_daily[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
  supp_ps_all_average <- supp_ps_all_average / (ps_average[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
}

test_p <- (supp_ps_all_daily_prior) | (supp_ps_all_daily)  |  (supp_ps_all_average) | (supp_ps_all_daily_prior_broad) | (supp_ps_all_daily_alt)
if(save_plots){
  ggsave("figures/supplement/all_dials.png",test_p,width=3.2,height=8,dpi=300,units="in")
  ggsave("figures/supplement/all_dials.pdf",test_p,width=3.2,height=8,units="in")
}

supp_ps_all_daily_alt_horiz <- ps_daily_alt[[1]] + 
  theme(plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="A")
for(i in seq(4,length(obs_times),by=3)){
  supp_ps_all_daily_alt_horiz <- supp_ps_all_daily_alt_horiz | (ps_daily_alt[[i]] + geom_hline(yintercept = 0,linetype="dashed",size=0.5,color="grey40"))
}

supp_p_use_dial_right <- ps_daily[[4]] +
  scale_x_continuous(expand=c(0,0),limits=c(-0.1,beta_max+0.3)) +
  scale_y_continuous(expand=c(0,0),limits=c(-beta_max-0.2,beta_max+0.2)) +
  geom_vline(xintercept=0,size=0.1)+
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=-yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_label(data=text_circle_top,aes(x=x,y=y,label=text),size=2) +
  geom_label(data=text_circle_bot,aes(x=x,y=y,label=text),size=2) +
  geom_text(data=data.frame(x=0.6,y=0.5,label="Growth"),aes(x=x,y=y,label=label),size=3) +
  geom_text(data=data.frame(x=0.6,y=-0.5,label="Decline"),aes(x=x,y=y,label=label),size=3) +
  theme_void() +
  theme(plot.margin=unit(c(0,-4,0,-4),unit="cm"),
        plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="C")

supp_p_use_dial_left <- ps_daily_alt[[4]] +
  scale_x_continuous(expand=c(0,0),limits=c(-0.1,beta_max+0.3)) +
  scale_y_continuous(expand=c(0,0),limits=c(-beta_max-0.2,beta_max+0.2)) +
  geom_vline(xintercept=0,size=0.1)+
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_segment(data=beta_lines,aes(x=x,y=y,xend=xend,yend=-yend,group=as.factor(group)),
               linetype="dashed",col="black",size=0.1) +
  geom_label(data=text_circle_top,aes(x=x,y=y,label=text),size=2) +
  geom_label(data=text_circle_bot,aes(x=x,y=y,label=text),size=2) +
  geom_text(data=data.frame(x=0.6,y=0.5,label="Growth"),aes(x=x,y=y,label=label),size=3) +
  geom_text(data=data.frame(x=0.6,y=-0.5,label="Decline"),aes(x=x,y=y,label=label),size=3) +
  theme_void() +
  theme(plot.margin=unit(c(0,-4,0,-4),unit="cm"),
        plot.tag = element_text(family="sans",size=10,face="bold")) +
  labs(tag="B")


plot_dat_supp  <- trajs_dat_all_comb_alt %>% filter(obs_time==56) %>%
  mutate(t = as.Date(t, origin="2020-02-01")) %>%
  mutate("Epidemic seed date"=ifelse(t0 < 60, "2020-02-01 to 2020-04-01","2020-04-01 to 2020-04-25"))

plot_dat_supp_t0 <- plot_dat_supp %>% group_by(`Epidemic seed date`) %>%
  summarize(t0_mean=mean(t0),t0_lower=quantile(t0, 0.025,na.rm=TRUE),t0_upper=quantile(t0, 0.975,na.rm=TRUE)) %>%
  mutate(t0_mean = as.Date(t0_mean, origin="2020-02-01"),
         t0_lower = as.Date(t0_lower, origin="2020-02-01"),
         t0_upper = as.Date(t0_upper, origin="2020-02-01"))
  
  
supp_p_traj <- plot_dat_supp %>%
  #filter(t >= as.Date("2020-05-03")-35) %>%
  ggplot() + 
  scale_color_manual(values=c("#3B4992FF","#EE0000FF")) +
  geom_line(aes(x=t,y=inc,color=`Epidemic seed date`,group=samp),size=0.1) +
  geom_pointrange(data=plot_dat_supp_t0,y=0.15,aes(xmin=t0_lower,xmax=t0_upper,x=t0_mean,col=`Epidemic seed date`)) +
  export_theme +
  theme(legend.position=c(0.3,0.9),
        panel.grid.major=element_line(color="grey40",size=0.1),
        plot.tag=element_text(size=10,face="bold")) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Incidence") +
  xlab("Date") +
  labs(tag="C")


chainwd_tmp <- paste0(chainwd,"/",runname_seir_alt,"/",56+seir_timeshift)
chain_seir_alt <- load_mcmc_chains(chainwd_tmp, parTab_seir,TRUE,1, mcmcPars_ct_seir["adaptive_period"],
                                   multi=TRUE,chainNo=FALSE,PTchain = TRUE)$list
coda::gelman.diag(chain_seir_alt)
coda::effectiveSize(chain_seir_alt)

chain1 <- res_seir_alt[[4]] %>%
  group_by(chain) %>%
  mutate(sampno=1:n())

supp_traces <- chain1[,c("sampno",unique(parTab_seir[which(parTab_seir$fixed == 0),"names"]),"chain")] %>%
  mutate(chain = as.factor(chain)) %>%
  pivot_longer(-c(sampno,chain)) %>%
  ggplot() +
  geom_line(aes(x=sampno,y=value,col=chain)) +
  scale_color_viridis_d() +
  facet_wrap(~name,scales="free_y")+
  scale_x_continuous(seq(1,5000,by=1000)) +
  export_theme +
  theme(plot.tag=element_text(size=10,face="bold")) +
  labs(tag="D")

figS12 <- (supp_ps_all_daily_alt_horiz) / ((supp_p_use_dial_left) | (supp_p_use_dial_right)) / (supp_p_traj) / (supp_traces) + plot_layout(heights=c(0.5,3,4,3))

if(save_plots){
  ggsave("figures/supplement/figS12_raw.png",figS12,width=7,height=8,dpi=300,units="in")
  ggsave("figures/supplement/figS12_raw.pdf",figS12,width=7,height=8,units="in")
}

########################################
## 11. Biobot data
########################################
biobot_dat <- read_csv("data/biobot_data_melted.csv")
p_inc_biobot <- p_inc + 
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(xintercept=as.Date("2020-04-15"),linetype="dashed") +
  geom_vline(xintercept=as.Date("2020-09-01"),linetype="dashed") +
theme(panel.grid.major=element_line(size=0.1,color="grey40"))  +
  labs(tag="B")

p_biobot <- ggplot(biobot_dat %>% filter(date <= "2020-12-01") %>% mutate(system2=paste0(system, ", 7 day average"))) +
  geom_point(aes(x=date,y=copies,color=system),size=0.5) +
  geom_errorbar(aes(x=date,ymin=copies-lower,ymax=copies+upper,color=system),width=0.1)+
  geom_line(aes(x=date,y=avg,color=system,linetype=system2))+
  #geom_point(aes(x=date,y=southern_cp),color="blue",size=0.5)+
  #geom_errorbar(aes(x=date,ymin=southern_cp-southern_lower,ymax=southern_cp+southern_upper,color="blue",width=0.1)+
  #geom_point(aes(x=date,y=northern_cp),color="darkorange",size=0.5)+
  #geom_errorbar(aes(x=date,ymin=northern_cp-northern_lower,ymax=northern_cp+northern_upper),color="darkorange",width=0.1)+
  #geom_line(aes(x=date,y=southern_7day),color="blue")+
  #geom_line(aes(x=date,y=northern_7day),color="darkorange")+
  geom_vline(xintercept=as.Date("2020-04-15"),linetype="dashed") +
  geom_vline(xintercept=as.Date("2020-09-01"),linetype="dashed") +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-12-01")),breaks="1 month",expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("Northern System"="darkorange","Southern System"="blue")) +
  scale_linetype_manual(values=c("Northern System, 7 day average"=1,"Southern System, 7 day average"=1))+
  guides(linetype=guide_legend(title=element_blank(),override.aes=list(color=c("Northern System, 7 day average"="darkorange","Southern System, 7 day average"="blue"))),
         color=guide_legend(title=element_blank())) +
  export_theme +
  theme(plot.margin=margin(0,0.5,0,0,"cm"),
        legend.position=c(0.5,0.8),
        legend.direction="horizontal") +
  labs(tag="A") + 
  ylab("copies/mL") +
  theme(panel.grid.major=element_line(size=0.1,color="grey40"),
        legend.background = element_rect(fill="white",color=NA)) +
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
p_biobot  
  
fig_S_biobot <-p_biobot/p_inc_biobot
if(save_plots){
  pdf("figures/supplement/biobot.pdf",height=5,width=7)
  fig_S_biobot
  dev.off()
  png("figures/supplement/biobot.png",height=5,width=7,units="in",res=300)
  fig_S_biobot
  dev.off()
}



## Pulling out growth rates for comparison
grs_daily_comb <- do.call("bind_rows",grs_daily)
grs_daily_comb$date <- as.Date(grs_daily_comb$t,origin="2020-03-08")


p_grs_compare <- ggplot() + 
  #geom_ribbon(data=gr_dat_top,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill="darkorange") +
  #geom_ribbon(data=gr_dat_bot,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill=AAAS_palette["green1"]) +
  #geom_line(data=gr_dat_median,aes(x=date,y=gr_mean,col=is_grow,group=index)) +
  geom_ribbon(data=trajs_dat_top,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_ribbon(data=trajs_dat_bot,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["purple1"]) +
  geom_line(data=trajs_dat_median %>% mutate(is_grow="Gaussian process estimate"),aes(x=date,y=median,col=is_grow,group=index)) +
  geom_errorbar(data=grs_daily_comb,aes(x=date,ymin=lower95,ymax=upper95,col="Cross-sectional Ct estimate"),width=3,alpha=0.25) +
  geom_errorbar(data=grs_daily_comb,aes(x=date,ymin=lower50,ymax=upper50,col="Cross-sectional Ct estimate"),width=1.5) +
  geom_point(data=grs_daily_comb,aes(x=date,y=median,col="Cross-sectional Ct estimate"),size=1) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  scale_x_date(expand=c(0,0),limits=as.Date(c("2020-03-09", "2020-12-01"), "%Y-%m-%d"), 
               breaks="1 month") +
  geom_hline(yintercept=0,linetype="dashed",col=AAAS_palette["grey1"]) +
  scale_color_manual(values=c(#"R(t), growing"="darkorange",
                         #"R(t), declining"=as.character(AAAS_palette["green1"]),
                         "Gaussian process estimate"=as.character(AAAS_palette["purple1"]),
                         "Cross-sectional Ct estimate"=as.character(AAAS_palette["green1"])))+
  export_theme +
  theme(legend.position=c(0.2,0.9),
        legend.title = element_blank(),
        axis.text.x=element_text(size=6),
        legend.text=element_text(size=6),
        plot.margin=unit(c(0.5,1,0.5,0.5),units="cm")) +
  xlab("Date") +
  ylab("Growth rate") 
if(save_plots){
  ggsave(plot=p_grs_compare,filename="figures/supplement/bwh_growths_compare.png",height=4,width=7,units="in",dpi=300)
  ggsave(plot=p_grs_compare,filename="figures/supplement/bwh_growths_compare.pdf",height=4,width=7,units="in",dpi=300)
}
