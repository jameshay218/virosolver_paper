## Main script 1: fit kinetics model based on expected parameters
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

devtools::load_all("~/Documents/GitHub/virosolver/")

######################################################
## 1. FIT 2-STAGE HINGE MODEL SUBJECTIVELY
######################################################
parTab <- read.csv("pars/partab_for_optim.csv")

## Dealing with log10 per ul
parTab[parTab$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
to_fit_probs <- read_csv("data/borremands_urt_pos.csv")
ages_observed <- c(-assumed_incu_period, to_fit_probs$time_since_onset) + assumed_incu_period
desired_probs <- (c(0,to_fit_probs$pos)/100) - 0.15
desired_probs <- pmax(desired_probs, 0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_hinge2 <- parTab$values
names(pars_hinge2) <- parTab$names

ages <- ages_observed

######################################################
## COST FUNCTION
cost_function_hinge <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1, ct_upper=11.4){
  pars1 <- pars_hinge2
  ## Change true 0, gamma variance, observation variance and peak viral load
  #pars1["true_0"] <- explore_pars[1]
  pars1["t_switch"] <- explore_pars[1]
  pars1["obs_sd"] <- explore_pars[2]
  pars1["viral_peak"] <- explore_pars[3]

  pars1["level_switch"] <- explore_pars[4]
  pars1["wane_rate2"] <- explore_pars[5]
  pars1["prob_detect"] <- explore_pars[6]

  peak <- viral_load_func(pars1,pars1["desired_mode"]+pars1["tshift"],convert_ct=TRUE)
  cost1 <- 0

  ## Also include highest recorded Ct
  if(use_ct_cost){
    cost1 <- (ct_upper - qgumbel(0.01,peak,pars1["obs_sd"]))^2
  }
  vls <- viral_load_func(pars1, 1:100, convert_ct=TRUE)
  obs <- c(0,vapply(ages[ages > 0], function(a) prop_detectable(a,pars1, vls),FUN.VALUE=numeric(1)))
  costs2 <- (obs-desired_probs)^2
  sum(costs2 + cost1*ct_cost_weight)
}
upper_bounds <- c(30,10,5,5,25,1)
lhs_pars <- randomLHS(1000,6)
test_pars <- t(apply(lhs_pars,1, function(x) x*upper_bounds))

res <- apply(test_pars, 1, function(x){
  cost_function_hinge(x)
  })
start_par <- test_pars[which.min(res),]
######################################################
## A) RUN OPTIM
######################################################
fit_hinge2 <- optimx(c(12,5,9,3.5,30,0.05), cost_function_hinge,
                     lower=c(0,1,0,3,0,.000001),
                     upper=c(25,10,25,10,50,1),
                     method="L-BFGS-B",
                   use_ct_cost=TRUE, ct_cost_weight=1)

pars1 <- pars_hinge2

pars1["t_switch"] <- as.numeric(fit_hinge2[1,1])
pars1["obs_sd"] <- as.numeric(fit_hinge2[1,2])
pars1["viral_peak"] <- as.numeric(fit_hinge2[1,3])
pars1["level_switch"] <- as.numeric(fit_hinge2[1,4])
pars1["wane_rate2"] <- as.numeric(fit_hinge2[1,5])
pars1["prob_detect"] <- as.numeric(fit_hinge2[1,6])
pars1["t_unit"] <- 1

######################################################
## B) PLOT AGAINST PROPORTION DETECTABLE
######################################################
parTab$values <- pars1
write_csv(parTab,"pars/partab_fitted_lower_sens.csv")

test_ages <- seq(0,55,by=1)
vls <- viral_load_func(pars1, test_ages, TRUE)
probs <- c(0,vapply(test_ages[test_ages > 0], function(a) prop_detectable(a,pars1, vls),FUN.VALUE=numeric(1)))

fitted_dat <- tibble(age=test_ages,prob=probs)

## Generate draws from our priors and get quantiles
sample_pars_from_prior <- function(pars){
  #sds <- c(2, 3, 5, 3, 0.5)
  sds <- c(1, 1, 3, 3, 1, 0.03)
  tmp_pars <- pars
  #tmp_pars["desired_mode"] <- runif(1,2,5)
  tmp_pars["obs_sd"] <- max(rnorm(1, pars["obs_sd"], sds[1]), 1)
  tmp_pars["viral_peak"] <- rnorm(1, pars["viral_peak"], sds[2])
  tmp_pars["wane_rate2"] <- rnorm(1,pars["wane_rate2"],sds[3])
  tmp_pars["t_switch"] <- rnorm(1,pars["t_switch"],sds[4])
  tmp_pars["level_switch"] <- max(rnorm(1,pars["level_switch"],sds[5]),pars["LOD"])
  beta1_mean <- pars["prob_detect"]
  beta1_sd <- sds[6]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  tmp_pars["prob_detect"] <- rbeta(1,beta_alpha,beta_beta)
  tmp_pars
}
n_samp1 <- 10000
tmp_trajs <- matrix(0, nrow=n_samp1, ncol=length(test_ages))
tmp_trajs2 <- matrix(0, nrow=n_samp1, ncol=length(test_ages))
all_pars <- matrix(nrow=n_samp1,ncol=length(pars1))
for(i in 1:n_samp1){
  tmp_pars <- sample_pars_from_prior(pars1)
  tmp_vls <- viral_load_func(tmp_pars, test_ages, TRUE)
  tmp_vls2 <- viral_load_func(tmp_pars, test_ages, TRUE)
  tmp_vls2_mean <- tmp_vls2 - tmp_pars["obs_sd"]*(-digamma(1))
  probs <- c(0,vapply(test_ages[test_ages > 0], function(a) prop_detectable(a,tmp_pars, tmp_vls),FUN.VALUE=numeric(1)))
  if(any(is.na(probs))) print(tmp_pars)
  tmp_trajs[i,] <- probs
  tmp_trajs2[i,] <- tmp_vls2
  all_pars[i,] <- tmp_pars
}
prob_quants <- t(apply(tmp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
prob_quants <- as.data.frame(prob_quants)
colnames(prob_quants) <- c("lower","mid_low","median","mid_high","upper")
prob_quants$t <- test_ages
ggplot(prob_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=mid_low,ymax=mid_high),alpha=0.5) +
  geom_line(aes(x=t,y=median))

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


ggplot(vl_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper), alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=mid_low,ymax=mid_high),alpha=0.5) +
  geom_line(aes(x=t,y=median))
######################################################
## C) PLOT DISTRIBUTION OF VIRAL LOADS
######################################################
n_samps <- 100000
vls <- viral_load_func(pars1, test_ages, TRUE)
line_dat <- tibble(mean_load=vls,t=test_ages)
violin_times <- seq(0,50,by=5)
hists <- matrix(ncol=n_samps,nrow=length(violin_times))
for(i in 1:length(violin_times)){
  omg <- viral_load_func(pars1, violin_times[i],TRUE)
  #hists[i,] <- rnorm(10000,omg, pars1["obs_sd"])
  hists[i,] <- rgumbel(n_samps,omg, pars1["obs_sd"])

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
                     sec.axis=sec_axis(trans=~(.*(-1)+40)/log2(10)+3,
                                       name=expression("Viral load ("*log["10"]*" RNA copies / mL)"),
                                       breaks=seq(-3,15,by=3))) +
  export_theme +
  theme(legend.position = c(0.7,0.9),panel.grid.major=element_line(size=0.05,color="grey40")) +
  coord_cartesian(ylim=c(pars1["intercept"]+10,0)) +
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

p3 <- ggplot(all_pars_melted) +
  geom_density(aes(x=value), fill=cbbPalette[4],alpha=0.5) +
  ylab("Prior density") +
  xlab("Value") +
  facet_wrap(~label, scales="free", ncol=2,labeller=label_parsed) +
  export_theme+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6)) +
  labs(tag="C")

p_main <- (p_left | p3) + plot_layout(widths=c(1.5,1))

ggsave("figures/supplement/viral_kinetics_lower.pdf",p_main,height=6,width=8)
ggsave("figures/supplement/viral_kinetics_lower.png",p_main,height=6,width=8)

