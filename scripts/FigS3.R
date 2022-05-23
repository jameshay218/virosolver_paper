########################################
## SCRIPT TO GENERATE FIG S3
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

AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="#1B1919FF")

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
        plot.margin = margin(0,0,0,0,unit="cm"),
        plot.tag = element_text(size=10,face="bold"),
        panel.grid.minor.x=element_blank())

main_theme2 <- theme(axis.text.x=element_text(size=7),
                     axis.text.y=element_text(size=7),
                     axis.title.x=element_text(size=8),
                     axis.title.y=element_text(size=8),
                     legend.text=element_text(size=7),
                     plot.tag=element_text(size=10,face="bold"))

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
complete_linelist <- virosolver::simulate_observations_wrapper(as.integer(seir_dynamics$incidence),times=times,
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

R0_times <- tibble(x=c(0,seir_pars1[c("t_switch1","t_switch2")]),label=c("R[0] == 2","R[0] == 0.7","R[0] == 1.5"))



## Look at all the various distributions over time
p_incu_overall <- simulated_viral_loads %>% ggplot() + 
  geom_histogram(aes(x=incu_period,y=..density..),fill="#008280FF",col="black") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(expand=c(0,0)) +
  xlab("Overall incubation period (days)") +
  ylab("Probability density") +
  main_theme +
  main_theme2 +
  theme(plot.background = element_blank(),
        panel.grid.major=element_line(size=0.01,color="grey70")) +
  labs(tag="A")
p_conf_overall <- simulated_viral_loads %>% ggplot() + 
  geom_histogram(aes(x=confirmation_delay,y=..density..),fill="#BB0021FF",col="black",binwidth=1) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,10,by=1)) +
  xlab("Overall testing delay (days)") +
  ylab("Probability density") +
  main_theme+
  main_theme2 +
  theme(plot.background = element_blank(),
        panel.grid.major=element_line(size=0.01,color="grey70")) +
  labs(tag="B")

p_simulation <- ggplot(data=seir_dynamics$seir_outputs) + 
  geom_bar(aes(x=step, y=inc, fill="Incidence"),stat="identity",col="#808180FF",size=0.25) +
  geom_line(aes(x=step,y=Rt*2000,color="R[t]"), size=1) + 
  geom_segment(data=R0_times,aes(x=x,xend=x,y=0,yend=4000),col="#008280FF",size=0.5,linetype='dashed') +
  theme_light() +
  scale_color_manual(values=c("R[t]" = "#3B4992FF"),labels=c(parse(text="R[t]"))) +
  scale_fill_manual(values=c("Incidence"="#808180FF")) +
  geom_hline(yintercept=1*2000,linetype="dashed",col="black",size=0.5) +
  geom_text(data=R0_times,aes(x=x,y=4300,label=label),size=3,parse=TRUE) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,4500),
                     name="Incidence",
                     sec.axis=sec_axis(trans=~.*1/2000, name=expression(R[t]))) + 
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  main_theme +
  main_theme2 + 
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.05,color="grey70")) + 
  labs(tag="C")


p_incu_time <- simulated_viral_loads %>% 
  ggplot() + 
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=0.5,linetype='dashed') +
  geom_smooth(aes(x=sampled_time,y=incu_period,fill="Grouped by sample date",col="Grouped by sample date"),size=0.5,alpha=0.25) +
  geom_smooth(aes(x=infection_time,y=incu_period,fill="Grouped by infection date",col="Grouped by infection date"),size=0.5,alpha=0.25) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  scale_fill_manual(name="group",values=c("Grouped by sample date"="#3B4992FF","Grouped by infection date"="#008280FF")) +
  scale_color_manual(name="group",values=c("Grouped by sample date"="#3B4992FF","Grouped by infection date"="#008280FF")) +
  xlab("Time") +
  ylab("Smoothed delay from infection\n to symptom onset (days)") +
  main_theme+
  main_theme2 +
  theme(legend.title=element_blank())+
  theme(legend.position=c(0.75,0.2),
        legend.direction = "horizontal",
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.05,color="grey70"))+
  labs(tag="D")
p_conf_time <- simulated_viral_loads %>% ggplot() + 
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=0.5,linetype='dashed') +
  geom_smooth(aes(x=sampled_time,y=confirmation_delay,fill="Grouped by sample date",color="Grouped by sample date"),size=0.5,alpha=0.25) +
  geom_smooth(aes(x=infection_time,y=confirmation_delay,fill="Grouped by infection date",color="Grouped by infection date"),size=0.5,alpha=0.25) +
  scale_fill_manual(name="group",values=c("Grouped by sample date"="#EE0000FF","Grouped by infection date"="#BB0021FF")) +
  scale_color_manual(name="group",values=c("Grouped by sample date"="#EE0000FF","Grouped by infection date"="#BB0021FF")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  xlab("Time") +
  ylab("Smoothed delay from symptom\n onset to test date (days)") +
  main_theme +
  main_theme2 +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(0.75,0.2),
        legend.direction = "horizontal",
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.05,color="grey70"))+
  labs(tag="E")

p_cts <- simulated_viral_loads %>% ggplot() + 
  geom_vline(data=R0_times,aes(xintercept=x),col="#008280FF",size=0.5,linetype='dashed') +
  geom_smooth(aes(x=sampled_time,y=ct_obs,fill="Grouped by sample date",color="Grouped by sample date"),size=0.5,alpha=0.25) +
  geom_smooth(aes(x=infection_time,y=ct_obs,fill="Grouped by infection date",color="Grouped by infection date"),size=0.5,alpha=0.25) +
  scale_fill_manual(name="group",values=c("Grouped by sample date"="#631879FF","Grouped by infection date"="#5F559BFF")) +
  scale_color_manual(name="group",values=c("Grouped by sample date"="#631879FF","Grouped by infection date"="#5F559BFF")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(breaks=seq(0,365,by=50),limits=c(0,365)) +
  xlab("Time") +
  ylab("Smoothed Ct values") +
  main_theme +
  main_theme2 +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(0.75,0.2),
        legend.direction = "horizontal",
        panel.grid.minor.x=element_blank(),
        panel.grid.major=element_line(size=0.05,color="grey70"))+
  labs(tag="F")


p_top <- p_incu_overall + p_conf_overall + plot_layout(nrow=1)

p_main <- ((p_top) / p_simulation / p_incu_time / p_conf_time / p_cts)

ggsave("~/Documents/GitHub/virosolver_paper/figures/supplement/FigS3_new.pdf",width=6.5,height=9)
