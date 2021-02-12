## Pre-script 1: fit kinetics model based on expected parameters
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(lazymcmc)
library(patchwork)
library(lhs)
library(optimx)
library(ggthemes)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop must move to the correct working directory to source the model functions
#main_wd <- "~/Documents/Harvard/Research/Infectious Diseases/COVID-19/PCR Cts/Code & Results/ct_inference"
main_wd <- "~/Documents/GitHub/virosolver_paper/"
setwd(main_wd)

source("code/plot_funcs.R")
source("code/priors.R")

devtools::load_all("~/Documents/GitHub/virosolver/")
SAVE_PLOTS <- TRUE
#run_ver <- "bwh"
run_ver <- "nh"

parTab <- read.csv("pars/partab_for_optim.csv")

if(run_ver == "bwh"){
  bwh_data <- read_csv("data/panther_Ct_20200403-20201110.csv")
  ## Ct value plateaus at 38 in BWH data, 33 in NH data
  parTab[parTab$names == "level_switch","values"] <- 38
  ct_lower <- range(bwh_data$ORF1ab_Ct,na.rm=TRUE)[1]
  ct_upper <- 32
} else {
  nh_data <- read_csv("data/nursing_home_data.csv") %>% filter(is_nursinghome==1)
  parTab[parTab$names == "level_switch","values"] <- 33
  ct_lower <- range(nh_data$N2, na.rm=TRUE)[1]
  ct_upper <- 30
}
######################################################
## 1. FIT 2-STAGE HINGE MODEL SUBJECTIVELY
######################################################
## Dealing with log10 per ul
parTab[parTab$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
to_fit_probs <- read_csv("data/borremands_urt_pos.csv")
ages_observed <- c(-assumed_incu_period, to_fit_probs$time_since_onset) + assumed_incu_period
ages <- ages_observed

desired_probs <- (c(0,to_fit_probs$pos)/100)
desired_probs <- pmax(desired_probs, 0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_hinge2 <- parTab$values
names(pars_hinge2) <- parTab$names


######################################################
## COST FUNCTION
cost_function_hinge <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1, 
                                ct_lower=12.1,ct_upper=32){

  pars1 <- pars_hinge2
  
  ## Pull out the parameters being optimized
  pars1["t_switch"] <- explore_pars[1]
  pars1["obs_sd"] <- explore_pars[2]
  pars1["viral_peak"] <- explore_pars[3]
  pars1["prob_detect"] <- explore_pars[6]
  pars1["sd_mod"] <- explore_pars[7]
  
  #pars1["level_switch"] <- explore_pars[4]
  #pars1["wane_rate2"] <- explore_pars[5]

  ## Find Cts at peak and at day 30
  peak <- viral_load_func(pars1,pars1["desired_mode"]+pars1["tshift"],convert_vl=FALSE)
  day30 <- viral_load_func(pars1,30,convert_vl=FALSE)
  
  ## Calculate cost function for peak and day 30 Cts
  cost1 <- 0
  cost3 <- 0
  ## Also include highest recorded Ct
  ## Lower 1% of distribution should be ct_upper at peak and
  ## upper 1% should be ct_lower at day 30
  if(use_ct_cost){
    cost1 <- (ct_lower - qgumbel(0.01,peak,pars1["obs_sd"]))^2
    cost3 <- (ct_upper - qgumbel(0.01,day30,pars1["obs_sd"]*pars1["sd_mod"]))^2
  }
  ## Calculate proportion detectable curve
  ## Get proportion detectable over time and fit to observations
  vls <- viral_load_func(pars1, 1:max(ages), convert_vl=FALSE)
  obs <- c(0,prop_detectable(ages[ages > 0], pars1, vls))
  costs2 <- (obs-desired_probs)^2
  sum(costs2 + cost1*ct_cost_weight + cost3)
}

######################################################
## A) RUN OPTIM
######################################################
fit_hinge2 <- optimx(c(10,5.5,30,40,40,0.2,0.8), cost_function_hinge,
                     lower=c(0,5,0,30,15,.000001,0.25),
                     upper=c(25,10,40,40,100,1,1),
                     method="L-BFGS-B",
                   use_ct_cost=TRUE, ct_cost_weight=1,
                   ct_lower=ct_lower,ct_upper=ct_upper)

pars1 <- pars_hinge2

pars1["t_switch"] <- as.numeric(fit_hinge2[1,1])
pars1["obs_sd"] <- as.numeric(fit_hinge2[1,2])
pars1["viral_peak"] <- as.numeric(fit_hinge2[1,3])
pars1["prob_detect"] <- as.numeric(fit_hinge2[1,6])
pars1["sd_mod"] <- as.numeric(fit_hinge2[1,7])
pars1["t_unit"] <- 1
#pars1["wane_rate2"] <- 20
#pars1["level_switch"] <- as.numeric(fit_hinge2[1,4])
#pars1["wane_rate2"] <- as.numeric(fit_hinge2[1,5])

######################################################
## B) PLOT AGAINST PROPORTION DETECTABLE
######################################################
test_ages <- seq(0,55,by=1)
vls <- viral_load_func(pars1, test_ages, FALSE)
probs <- c(0,prop_detectable(test_ages[test_ages > 0], pars1, vls))
fitted_dat <- tibble(age=test_ages,prob=probs)

## Generate draws from our priors and get quantiles
sample_pars_from_prior <- function(pars){
  sds <- sds_exp
  ## Draw parameters from normal priors
  tmp_pars <- pars
  tmp_pars["obs_sd"] <- max(rnorm(1, pars["obs_sd"], sds["obs_sd"]), 1)
  tmp_pars["viral_peak"] <- rnorm(1, pars["viral_peak"], sds["viral_peak"])
  tmp_pars["t_switch"] <- rnorm(1,pars["t_switch"],sds["t_switch"])
  
  #tmp_pars["wane_rate2"] <- rnorm(1,pars["wane_rate2"],sds["wane_rate2"])
  #tmp_pars["level_switch"] <- min(rnorm(1,pars["level_switch"],sds["level_switch"]), pars["intercept"])
  
  beta1_mean <- pars["prob_detect"]
  beta1_sd <- sds["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  tmp_pars["prob_detect"] <- rbeta(1,beta_alpha,beta_beta)
  tmp_pars
}

## Draw 10000 samples from prior
n_samp1 <- 10000
tmp_trajs <- matrix(0, nrow=n_samp1, ncol=length(test_ages))
tmp_trajs2 <- matrix(0, nrow=n_samp1, ncol=length(test_ages))
all_pars <- matrix(nrow=n_samp1,ncol=length(pars1))

## For 10000 samples
for(i in 1:n_samp1){
  ## Get parameters
  tmp_pars <- sample_pars_from_prior(pars1)
  ## Get trajectory Ct scale
  tmp_vls <- viral_load_func(tmp_pars, test_ages, FALSE)
  ## Get trajectory viral load scale
  tmp_vls2 <- viral_load_func(tmp_pars, test_ages, TRUE)
  ## Get proportion detectable, but 0 on day 0
  probs <- c(0,prop_detectable(test_ages[test_ages > 0], tmp_pars, tmp_vls))
  if(any(is.na(probs))) print(tmp_pars)
  ## Save
  tmp_trajs[i,] <- probs
  tmp_trajs2[i,] <- tmp_vls
  all_pars[i,] <- tmp_pars
}

## Get median, 95 and 50% quantiles
prob_quants <- t(apply(tmp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
prob_quants <- as.data.frame(prob_quants)
colnames(prob_quants) <- c("lower","mid_low","median","mid_high","upper")
prob_quants$t <- test_ages

tmp_trajs_melted <- reshape2::melt(tmp_trajs)
colnames(tmp_trajs_melted) <- c("samp","t","prob")
tmp_trajs_melted$t <- test_ages[tmp_trajs_melted$t]

vl_quants <- t(apply(tmp_trajs2, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
vl_quants <- as.data.frame(vl_quants)
colnames(vl_quants) <- c("lower","mid_low","median","mid_high","upper")
vl_quants$t <- test_ages

tmp_vl_melted <- reshape2::melt(tmp_trajs2)
colnames(tmp_vl_melted) <- c("samp","t","vl")
tmp_vl_melted$t <- test_ages[tmp_vl_melted$t]

## Simple plots to check
p_detect_tmp <- ggplot(prob_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=mid_low,ymax=mid_high),alpha=0.5) +
  geom_line(aes(x=t,y=median))

p_vl_tmp <- ggplot(vl_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper), alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=mid_low,ymax=mid_high),alpha=0.5) +
  scale_y_continuous(trans="reverse") +
  coord_cartesian(ylim=c(40,10)) +
  geom_line(aes(x=t,y=median))

######################################################
## C) PLOT DISTRIBUTION OF VIRAL LOADS
######################################################
## Lots of samples from prior
n_samps <- 100000
## Get point estimate trajectory
vls <- viral_load_func(pars1, test_ages, FALSE)
line_dat <- tibble(mean_load=vls,t=test_ages)

## Will be simulating observations at each of these times
violin_times <- c(1,seq(5,50,by=5))
hists <- matrix(ncol=n_samps,nrow=length(violin_times))

## Standard deviation changes over time, need to choose right one
t_switch <-  pars1["t_switch"] + pars1["desired_mode"] + pars1["tshift"]
sd_mod <- rep(pars1["sd_mod"], max(test_ages))
unmod_vec <- 1:min(t_switch,max(ages))
sd_mod[unmod_vec] <- 1
decrease_vec <- (t_switch+1):(t_switch+pars1["sd_mod_wane"])
sd_mod[decrease_vec] <- 1 - ((1-pars1["sd_mod"])/pars1["sd_mod_wane"])*seq_len(pars1["sd_mod_wane"])

for(i in 1:length(violin_times)){
  vl_tmp <- viral_load_func(pars1, violin_times[i],FALSE)
  hists[i,] <- rgumbel(n_samps,vl_tmp, pars1["obs_sd"]*sd_mod[test_ages[violin_times[i]+1]+1])

}
hists[hists > pars1["intercept"]] <- NA# pars1["intercept"]
hists <- reshape2::melt(hists)
colnames(hists) <- c("t","i","value")
hists$t <- violin_times[hists$t]
hists <- as_tibble(hists)

line_dat$lower <- qgumbel(0.975, mu=line_dat$mean_load, sigma=pars1["obs_sd"],lower.tail = TRUE, log.p=FALSE)
line_dat$upper <- qgumbel(0.025, mu=line_dat$mean_load, sigma=pars1["obs_sd"],lower.tail = TRUE, log.p=FALSE)

line_dat <- line_dat %>% left_join(fitted_dat %>% rename(t=age))

hists <- hists %>% left_join(fitted_dat %>% rename(t=age))

p1 <-  ggplot(data=line_dat) +
  geom_violin(data=hists, aes(x=t, y=value, group=t, fill=prob),scale="width",
              draw_quantiles = c(0.025,0.5,0.975)) +
  geom_ribbon(data=vl_quants, aes(x=t,ymin=lower,ymax=upper),alpha=0.25)+
  geom_line(data=tmp_vl_melted %>% filter(samp %in% sample(unique(tmp_trajs_melted$samp), 25)),
            aes(x=t,y=vl,group=samp),alpha=0.75,size=0.1)+
  scale_fill_gradient(name="Proportion detectable:", low="blue", high="red",
                      breaks=seq(0,1,by=.25), labels=seq(0,1,by=.25), limits=c(-.05,1.05),
                      guide=guide_colorbar(direction="horizontal",title.position="top",barwidth=5,barheight=1,
                                           ticks=FALSE)) +
  geom_line(aes(x=t, y=mean_load), col=cbbPalette[4],size=1.2) +
  #geom_point(aes(x=t, y=mean_load), col=cbbPalette[4],size=1) +
  #theme_light() +
  scale_x_continuous(breaks=seq(0,55,by=5)) +
  scale_y_continuous(breaks=seq(0,60,by=10), trans="reverse", expand=expansion(add=c(0,0)),
                     sec.axis=sec_axis(trans=~(.*(-1)+40)/log2(10)+pars1["LOD"],
                                       name=expression("Viral load ("*log["10"]*" RNA copies / mL)"),
                                       breaks=seq(-3,15,by=3))) +
  export_theme +
  theme(legend.position = c(0.7,0.9),panel.grid.major=element_line(size=0.05,color="grey40")) +
  coord_cartesian(ylim=c(pars1["intercept"],0)) +
  labs(x="Time since infection (days)", y="Cycle threshold (Ct) value") +
  geom_hline(yintercept=pars1["intercept"], lty="dashed")+
  labs(tag="A")

p2 <- ggplot(data=fitted_dat) +
  geom_ribbon(data=prob_quants, aes(x=t,ymin=lower,ymax=upper),alpha=0.1)+
  geom_line(data=tmp_trajs_melted %>% filter(samp %in% sample(unique(tmp_trajs_melted$samp), 25)),
            aes(x=t,y=prob,group=samp),alpha=0.5,size=0.1)+
  geom_line(aes(x=age, y=prob, col=prob),size=1.2) +
  geom_pointrange(data=to_fit_probs,aes(x=time_since_onset+assumed_incu_period, y=pos/100,
                                        ymin=lower_95/100,ymax=upper_95/100), col="black", size=0.1) +
  scale_color_gradient(name="Proportion detectable:", low="blue", high="red",
                       breaks=seq(0,1,by=.25), labels=seq(0,1,by=.25), limits=c(-.05,1.05)) +
  #theme_light() +
  export_theme +
  theme(legend.position="none",panel.grid.major=element_line(size=0.05,color="grey40")) +
  scale_x_continuous(breaks=seq(0,55,by=5)) +
  labs(x="Time since infection (days)", y="Proportion detectable")+
  labs(tag="B")

p_left <- p1/p2

colnames(all_pars) <- parTab$names
all_pars_melted <- reshape2::melt(all_pars)
colnames(all_pars_melted) <- c("samp","par","value")
parTab[parTab$names == "beta","fixed"] <- 1
all_pars_melted <- all_pars_melted %>% filter(par %in% parTab[parTab$fixed == 0, "names"])


par_key <- c("obs_sd"="sigma[obs]",
             "viral_peak"="VL[peak]",
             "t_switch"="t[switch]",
             "level_switch"="VL[switch]",
             "wane_rate2"="t[LOD]",
             "prob_detect"="p[addl]"
             )
all_pars_melted$par <- as.character(all_pars_melted$par)
all_pars_melted$label <- par_key[all_pars_melted$par]
all_pars_melted$label <- factor(all_pars_melted$label,
                                levels=c("sigma[obs]","VL[peak]","t[switch]",
                                         "VL[switch]","t[LOD]","p[addl]"))

p3 <- ggplot(all_pars_melted %>% filter(par %in% c("viral_peak","obs_sd",
                                                   "t_switch","prob_detect"))) +
  geom_density(aes(x=value), fill=cbbPalette[4],alpha=0.5) +
  ylab("Prior density") +
  xlab("Value") +
  facet_wrap(~label, scales="free", ncol=2,labeller=label_parsed) +
  export_theme+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6)) +
  labs(tag="C")

p_main <- (p_left | p3) + plot_layout(widths=c(1.5,1))
p_main
if(SAVE_PLOTS){
  if(run_ver == "bwh"){
    ggsave("figures/supplement/viral_kinetics_bwh.pdf",p_main,height=6,width=8)
    ggsave("figures/supplement/viral_kinetics_bwh.png",p_main,height=6,width=8)
    parTab$values <- pars1
    write_csv(parTab,"pars/partab_fitted_bwh.csv")
  } else {
    ggsave("figures/supplement/viral_kinetics_nh.pdf",p_main,height=6,width=8)
    ggsave("figures/supplement/viral_kinetics_nh.png",p_main,height=6,width=8)
    parTab$values <- pars1
    write_csv(parTab,"pars/partab_fitted_nh.csv")
  }
}
