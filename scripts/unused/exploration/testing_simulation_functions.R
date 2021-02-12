library(tidyverse)
library(ggplot2)
library(patchwork)
library(lazymcmc)
library(extraDistr)
library(ggthemes)

HOME_WD <- "~/Documents/GitHub/"
#HOME_WD <- "~"

set.seed(2)

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


model_pars <- read.csv(paste0(HOME_WD,"virosolver_paper/pars/massachusetts/partab_seir_model.csv"))
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulation parameters
population_n <- 1000000
times <- 0:200
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 

sim <- simulate_seir_process(pars, times, population_n)
plot(times,sim$raw_incidence/population_n)

sim <- solveSEIRModel_rlsoda_wrapper(pars, times)
lines(times,sim,col="blue")

sim <- simulate_seir_process(pars, times, population_n)
complete_linelist <- virosolver::simulate_observations_wrapper(sim$raw_incidence,times=times,symp_frac=0.35,population_n=population_n)

inc1 <- complete_linelist %>% group_by(infection_time) %>% tally()
set.seed(10)
lines(inc1$infection_time,inc1$n/population_n,col="red")
samp_times <- c(80,90,100)
tvr <- data.frame(prob=c(0.001,0.001/2,0.001/4),t=samp_times)
#tvr <- data.frame(prob=c(0.001,0.001,0.001),t=samp_times)

LL_R <- simulate_reporting(complete_linelist,
                           timevarying_prob=tvr,
                           solve_times=times,
                           symptomatic=FALSE)
#LL_R1<- LL_R$sampled_individuals %>% filter(sampled_time == 50) %>% sample_frac(1)
#LL_R2 <- LL_R$sampled_individuals %>% filter(sampled_time == 60) %>% sample_frac(0.5)
#LL_R3 <- LL_R$sampled_individuals %>% filter(sampled_time == 70) %>% sample_frac(0.25)
#LL_R_all <- bind_rows(LL_R1, LL_R2, LL_R3)
LL_R_all <- LL_R$sampled_individuals
SVL <- simulate_viral_loads_wrapper(LL_R_all,
                                    kinetics_pars=pars)

#SVL1 <- SVL %>% filter(sampled_time == 80) %>% sample_frac(1)
#SVL2 <- SVL %>% filter(sampled_time == 90) %>% sample_frac(0.5)
#SVL3 <- SVL %>% filter(sampled_time == 100) %>% sample_frac(0.25)
#SVL <- bind_rows(SVL1, SVL2, SVL3)

tmp <- LL_R$sampled_individuals %>% group_by(infection_time) %>% tally() 
lines(tmp$infection_time,tmp$n/population_n,col="green")


tmp <- SVL %>% group_by(infection_time) %>% tally() 
lines(tmp$infection_time,tmp$n/population_n,col="green")

vl <- viral_load_func(pars, 0:90)
prop <- tibble(age=0:90, p=prop_detectable(1:91,pars,vl))

beta1 <- ((6/pars["obs_sd"])/pi^2)


kinetics_pars <- pars
test_ages <- 0:1000
t_switch <-  kinetics_pars["t_switch"] + kinetics_pars["desired_mode"] + kinetics_pars["tshift"]
sd_mod <- rep(kinetics_pars["sd_mod"], max(test_ages))
unmod_vec <- 1:(min(t_switch,max(test_ages)))
sd_mod[unmod_vec] <- 1
decrease_vec <- (t_switch+1):(t_switch+kinetics_pars["sd_mod_wane"])
sd_mod[decrease_vec] <- 1 - ((1-kinetics_pars["sd_mod"])/kinetics_pars["sd_mod_wane"])*seq_len(kinetics_pars["sd_mod_wane"])

vl <- viral_load_func(pars, 0:90)
vl_dat <- tibble(t=0:90,vl=vl,mean1=vl + pars["obs_sd"]*(-digamma(1)),lower=qgumbel(0.95,vl,pars["obs_sd"]*sd_mod[1:91]),upper=qgumbel(0.05,vl,pars["obs_sd"]*sd_mod[1:91]))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mean_dat <- SVL %>% filter(last_detectable_day > sampled_time,days_since_infection >= 0,!is.na(infection_time)) %>% 
  group_by(days_since_infection) %>% summarize(y=mean(ct_obs_sim),lower=quantile(ct_obs_sim,0.05),upper=quantile(ct_obs_sim,0.95))

SVL %>% filter(last_detectable_day >= sampled_time,days_since_infection >= 0,days_since_infection < 60) %>% ggplot() + 
  geom_jitter(aes(y=ct_obs_sim,x=days_since_infection),alpha=0.1,size=0.1) + 
  geom_line(data=mean_dat,aes(x=days_since_infection,y=y),col="green") +
  geom_line(data=mean_dat,aes(x=days_since_infection,y=lower),col="green") +
  geom_line(data=mean_dat,aes(x=days_since_infection,y=upper),col="green") +
  #geom_line(data=vl_dat,aes(x=t,y=vl),col="red") +
  geom_line(data=vl_dat,aes(x=t,y=mean1),col="red",linetype="dashed") +
  geom_line(data=vl_dat,aes(x=t,y=lower),col="red",linetype="dashed") +
  geom_line(data=vl_dat,aes(x=t,y=upper),col="red",linetype="dashed") +
  scale_y_continuous(trans="reverse")


data_detect <- SVL %>% filter(days_since_infection >= 0) %>% mutate(is_detectable=ct_obs < 40) %>% group_by(days_since_infection) %>% 
  summarize(n=n(),pos=sum(is_detectable)) %>% mutate(prob=pos/n)



vl <- viral_load_func(pars, test_ages)
prop_detect <- tibble(age=test_ages,y=(sapply(test_ages,FUN = function(a) prop_detectable_cpp(a,vl[a+1],pars["obs_sd"]*sd_mod[a+1],pars["intercept"],
                                                                                              pars["t_switch"]+pars["tshift"]+pars["desired_mode"],pars["prob_detect"]))))
 

ggplot(data_detect) + geom_line(aes(x=days_since_infection,y=prob)) + coord_cartesian(xlim=c(0,100)) +
#  geom_line(data=prop,aes(x=age,y=p),col="red") +
  geom_line(data=prop_detect,aes(x=age,y=y),col="blue")


## Proportion positive over time

prob_infection_tmp <- solveSEIRModel_rlsoda_wrapper(pars,times)
preds <- pred_dist_wrapper(seq(0,40,by=1),obs_times=2:200,1:200,pars,prob_infection_tmp)

SVL %>% mutate(is_detectable=ct_obs < 40) %>% group_by(sampled_time) %>% 
  summarize(pos=sum(is_detectable),n=n(),prev=pos/n) %>% 
  ggplot() + geom_line(aes(x=sampled_time,y=prev)) +
  geom_line(data=preds %>% filter(ct == 40),aes(x=t,y=1-density),col="red")

cts_for_hist <- 
  SVL %>% filter(sampled_time %in% seq(0,200,by=10)) %>%
  filter(ct_obs < 40)

cts_n <- cts_for_hist %>% group_by(sampled_time) %>% tally() %>% rename(t=sampled_time)

ggplot(cts_for_hist) +
  geom_histogram(aes(x=ct_obs),binwidth=1,boundary=0) +
  geom_line(data=preds %>% filter(t %in% seq(0,200,by=10),ct < 40) %>% 
              group_by(t) %>% mutate(density=density/sum(density))%>%
              left_join(cts_n) %>% mutate(pred=density*n) %>% 
              rename(sampled_time=t),aes(x=ct,y=density*n),col="red") +
  facet_wrap(~sampled_time,scales="free_y")


## Test fitting
test_data <- SVL %>% filter(sampled_time %in%  samp_times) %>%
  dplyr::select(sampled_time, ct_obs) %>%
  rename(t=sampled_time,ct=ct_obs)

test_data %>% filter(ct < 40) %>% ggplot() + geom_histogram(aes(x=ct)) + facet_wrap(~t)

## Parameter table for this run
parTab <- read.csv(paste0(HOME_WD,"virosolver_paper/pars/massachusetts/partab_seir_model.csv"),stringsAsFactors=FALSE)
means <- parTab$values
names(means) <- parTab$names

f <- create_posterior_func(parTab,test_data,prior_func_hinge_seir,solveSEIRModel_rlsoda_wrapper,solve_ver="likelihood",use_pos=FALSE)
f(parTab$values)

parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(test_data$t)



n_temperatures <- 5
mcmcPars_ct <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                    "thin"=10,"adaptive_period"=20000,"save_block"=1000,"temperature" = seq(1,101,length.out=n_temperatures),
                    "parallel_tempering_iter" = 5,"max_adaptive_period" = 20000, 
                    "adaptiveLeeway" = 0.2, "max_total_iterations" = 30000)

covMat <- diag(nrow(parTab))
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
parTab_list <- rep(list(parTab),n_temperatures)
mvrPars <- rep(list(mvrPars), n_temperatures)

output <- run_MCMC(parTab=parTab_list,
                   data=test_data,
                   INCIDENCE_FUNC=solveSEIRModel_rlsoda_wrapper,
                   PRIOR_FUNC = prior_func_hinge_seir,
                   solve_likelihood=TRUE,
                   mcmcPars=mcmcPars_ct,
                   filename="~/Documents/GitHub/virosolver_paper/mcmc_chains/test",
                   CREATE_POSTERIOR_FUNC=create_posterior_func,
                   mvrPars=mvrPars,
                   OPT_TUNING=0.2,
                   use_pos=FALSE,
                   t_dist=NULL)

chain1 <- read.csv(output$file)
chain1 <- chain1[chain1$sampno > 10000,]
chain1$chain <- 1

p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
  mutate(chain = as.factor(chain)) %>%
  pivot_longer(-c(sampno,chain)) %>%
  ggplot() +
  geom_line(aes(x=sampno,y=value,col=chain)) +
  facet_wrap(~name,scales="free_y")+
  scale_x_continuous(breaks=seq(min(chain1$sampno),max(chain1$sampno),by=2000)) +
  export_theme

## Posterior density plots
p_densities <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
  mutate(chain = as.factor(chain)) %>%
  pivot_longer(-c(sampno,chain)) %>%
  ggplot() +
  geom_density(aes(x=value,fill=chain),alpha=0.25) +
  geom_vline(data=parTab[which(parTab$fixed == 0),] %>% rename(name=names), aes(xintercept=values)) +
  facet_wrap(~name,scales="free") +
  export_theme

## Get predicted incidence trends
predictions <- plot_prob_infection(chain1, 100, solveSEIRModel_rlsoda_wrapper,
                                   times,
                                   obs_dat=test_data)
p1 <- predictions$plot + geom_line(data=inc1,aes(x=infection_time,y=n/population_n),col="red")

## Get fits to observed Ct distribution
model_func <- create_posterior_func(parTab,test_data,NULL,solveSEIRModel_rlsoda_wrapper,"model")
p2 <- plot_distribution_fits(chain1, test_data, model_func,100,pos_only=FALSE)

samps <- sample(unique(chain1$sampno),100)
all_prop_detects <- NULL
all_rt_ests <- NULL

for(i in 1:100){
  pars <- get_index_pars(chain1, samps[i])
  vl <- viral_load_func(pars, test_ages)
  prop_detect_tmp <- tibble(samp=i,age=test_ages,y=(sapply(test_ages,FUN = function(a) prop_detectable_cpp(a,vl[a+1],pars["obs_sd"]*sd_mod[a+1],pars["intercept"],
                                                                                                pars["t_switch"]+pars["tshift"]+pars["desired_mode"],pars["prob_detect"]))))
  all_prop_detects[[i]] <- prop_detect_tmp
  
  
  init <- c(1-pars["I0"],0,pars["I0"],0,0)
  seir_pars <- c("beta"=pars["R0"]*(1/pars["infectious"]),"sigma"=1/pars["incubation"],"gamma"=1/pars["infectious"])
  names(seir_pars) <- c("beta","sigma","gamma")
  Rt_est <- data.frame(t=times,Rt=solveSEIRModel_rlsoda(times,init,seir_pars,compatible = TRUE)[,2]*pars["R0"],samp=i)
  all_rt_ests[[i]] <- Rt_est
}

all_prop_detects <- do.call("bind_rows",all_prop_detects)
all_rt_ests <- do.call("bind_rows",all_rt_ests)

ggplot(all_rt_ests) + geom_line(aes(x=t,y=Rt,group=samp),size=0.1) + 
  geom_line(data=sim$seir_outputs %>% filter(variable=="Rt"),aes(x=time,y=value),col="red")

ggplot(data_detect) + 
  geom_line(aes(x=days_since_infection,y=prob),col="red") + 
  coord_cartesian(xlim=c(0,100)) +
  #  geom_line(data=prop,aes(x=age,y=p),col="red") +
  geom_line(data=all_prop_detects,size=0.1,aes(x=age,y=y,group=samp)) +
  geom_line(data=prop_detect,aes(x=age,y=y),col="blue")
