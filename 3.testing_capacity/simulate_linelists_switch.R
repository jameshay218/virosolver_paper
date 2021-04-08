########################################
## SCRIPT TO SIMULATE LINE LIST DATA
## 2nd November 2020
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
library(lazymcmc)
library(extraDistr)
devtools::load_all("~/Documents/GitHub/virosolver")

## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/"
setwd(paste0(HOME_WD,"/virosolver_paper/"))

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

set.seed(1)


main_theme <- theme_classic() +
  theme(axis.ticks.length.x = unit(0.2,"cm"),
        text=element_text(family="sans"),
        plot.margin = margin(-1,-1,-1,-1,unit="cm"),
        plot.tag = element_text(size=10,face="bold"),
        panel.grid.minor.x=element_blank())

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
model_pars <- read.csv("pars/massachusetts/partab_seir_model.csv")
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulation parameters
population_n <- 1000000
times <- 0:365
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 

########################################
## 3. Full simulated line list
########################################
## Simulate SEIR dynamics, incidence, growth rates, Rt etc
## Use "ode" for deterministic
## Use "odin" for stochastic
seir_pars <- read.csv("pars/massachusetts/partab_seir_switch_model.csv")
seir_pars1 <- seir_pars$values
names(seir_pars1) <- seir_pars$names
seir_pars1["R0_1"] <- 2
seir_pars1["R0_2"] <- 0.7
seir_pars1["R0_3"] <- 1.5
seir_pars1["t_switch1"] <- 80
seir_pars1["t_switch2"] <- 150
seir_pars1["t0"] <- 0
seir_pars1["I0"] <- 0.00001

seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=seir_pars1, ver="odin",switch_model=TRUE)
## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                                  population_n=population_n)

########################################
## 4. Simulate observation process
########################################
## Choose which process to use
## You could combine these as well, just make sure that each individual is only sampled once.
## ie. use frac_report_increase on observed_individuals; subset observed_individuals by individuals 
## who weren't sampled with frac_report_increase; use prob_report_increase on the remaining individuals

## Simulate time-varying reporting **fraction**
## This represents a fraction of the population being tested each day
frac_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.00001,end_prob=0.001,
                                                                   growth_rate=0.1,switch_point=120), ver="increase")
frac_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.0001,end_prob=0.00001,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Simulate time-varying reporting **probability**
## This gives the probability of an individual getting reported if they have symptoms
prob_report_increase <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.01,end_prob=0.1,
                                                                   growth_rate=0.05,switch_point=250), ver="increase")
prob_report_decrease <- tibble(t=times_extended,prob=logistic_func(times_extended,start_prob=0.1,end_prob=0.01,
                                                                   growth_rate=0.05,switch_point=250), ver="decrease")

## Choose which version you want
reporting_prob <- 1
#observed_linelist <- simulate_reporting(complete_linelist, 
#                                        frac_report=reporting_prob,
#                                        timevarying_prob=NULL,
#                                        solve_times=times, 
#                                        symptomatic=FALSE)

## How to use other versions
# ## Do simulation recovery
# observed_indivs_increasing <- simulate_reporting(complete_linelist, timevarying_prob=frac_report_increase, 
#                                                  solve_times=times, symptomatic=FALSE)
# observed_indivs_decreasing <- simulate_reporting(complete_linelist, timevarying_prob=frac_report_decrease, 
#                                                  solve_times=times, symptomatic=FALSE)
# 
# ## Symptomatic surveillance
 observed_indivs_symptom <- simulate_reporting(complete_linelist %>% filter(is_infected == 1), frac_report=reporting_prob, 
                                               solve_times=times, symptomatic=TRUE)
# observed_indivs_symptom_increasing <- simulate_reporting(complete_linelist, timevarying_prob=prob_report_increase, 
#                                                          solve_times=times, symptomatic=TRUE)
# observed_indivs_symptom_decreasing <- simulate_reporting(complete_linelist, timevarying_prob=prob_report_decrease, 
#                                                          solve_times=times, symptomatic=TRUE)


########################################
## 5. Simulate Ct values
########################################
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_symptom$sampled_individuals,
                                                      kinetics_pars=pars)
simulated_viral_loads <- simulated_viral_loads %>% mutate(period = "Second decline",
                                 period= ifelse(sampled_time < 250 & sampled_time > 150,"Second growth",period),
                                 period= ifelse(sampled_time < 150 & sampled_time > 75,"First decline",period),
                                 period= ifelse(sampled_time < 75,"First growth",period))

simulated_viral_loads <- simulated_viral_loads %>%
  mutate(days_since_onset = sampled_time - onset_time)

obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)

R0_times <- tibble(x=c(0,seir_pars1[c("t_switch1","t_switch2")]),label=c("R[0] == 2","R[0] == 0.9","R[0] == 1.5"))

p1 <- ggplot(data=seir_dynamics$seir_outputs) + 
  geom_bar(aes(x=step, y=inc, fill="Incidence"),stat="identity",col="#808180FF",size=0.25) +
  geom_line(aes(x=step,y=Rt*2000,color="R[t]"), size=1) + 
  theme_light() +
  scale_color_manual(values=c("R[t]" = "#3B4992FF"),labels=c(parse(text="R[t]"))) +
  scale_fill_manual(values=c("Incidence"="#808180FF")) +
  geom_hline(yintercept=1*2000,linetype="dashed",col="black",size=1) +
  geom_segment(data=R0_times,aes(x=x,xend=x,y=0,yend=4000),col="#008280FF",size=1) +
  geom_text(data=R0_times,aes(x=x,y=4250,label=label),size=3,parse=TRUE) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,4350),
                     name="Incidence",
                     sec.axis=sec_axis(trans=~.*1/2000, name=expression(R[t]))) + 
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  main_theme +
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.1,colour="grey40")) + 
  labs(tag="A")


p2 <- simulated_viral_loads %>% 
  ggplot() + 
  geom_smooth(aes(x=sampled_time,y=days_since_infection),fill="#008B45FF",col="#008B45FF") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=1) +
  ylab("Smoothed days from\ninfection to testing") +
  main_theme + 
    theme(legend.title=element_blank(), 
          legend.position=c(0.8,0.8),
          axis.line.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major=element_line(size=0.1,colour="grey40")) +
  labs(tag="B")
p2_alt <- simulated_viral_loads %>% 
  ggplot() + 
  geom_smooth(aes(x=sampled_time,y=days_since_onset),fill="#008B45FF",col="#008B45FF") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=1) +
  ylab("Smoothed days from\ninfection to testing") +
  main_theme + 
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.1,colour="grey40")) +
  labs(tag="B")

p3 <- obs_dat %>% 
  filter(ct < 40) %>% 
  ggplot() +
  geom_smooth(aes(x=t,y=ct),col="#631879FF",fill="#631879FF") +
  scale_y_continuous(trans="reverse",expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=1) +
  ylab("Smoothed Ct values") +
  xlab("Days since start of outbreak") +
  main_theme + 
  main_theme + 
  theme(legend.title=element_blank(), 
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.1,colour="grey40")) +
  labs(tag="C")

p_main <- p1 / p2 / p3
ggsave("figures/supp_symptom_surveillance.pdf",p_main,height=8,width=8,units="in",device=cairo_pdf)
ggsave("figures/supp_symptom_surveillance.png",p_main,height=8,width=8,units="in",dpi=300)
