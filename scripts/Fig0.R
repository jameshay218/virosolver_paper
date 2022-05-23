########################################
## SCRIPT TO GENERATE COMPONENTS OF FIG 0
## 17th May 2021

########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lazymcmc)
library(extraDistr)
library(ggpubr)
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
        axis.ticks = element_line(color="black"),
        text=element_text(family="sans"),
        plot.margin = margin(0,0,0,0,unit="cm"),
        plot.tag = element_text(size=10,face="bold"),
        panel.grid.minor.x=element_blank())

main_theme2 <- theme(axis.text.x=element_text(size=7.5),
                     axis.text.y=element_text(size=7.5),
                     axis.title.x=element_text(size=8),
                     axis.title.y=element_text(size=8),
                     legend.text=element_text(size=7.5),
                     plot.tag=element_text(size=10,face="bold"))

save_plots <- TRUE


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
seir_pars <- read.csv("pars/massachusetts/partab_seir_model.csv")
seir_pars1 <- seir_pars$values
names(seir_pars1) <- seir_pars$names
seir_pars1["R0"] <- 1.8
seir_pars1["t0"] <- 0
seir_pars1["I0"] <- 0.0001

set.seed(2)
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=seir_pars1, ver="odin",switch_model=FALSE)

seir_dynamics$plot

y <- (seir_dynamics$incidence/population_n)
fit1 <- smooth.spline(y)
fit <- smooth.spline(seir_dynamics$incidence/population_n,spar=0.6)


midpoint <- which.max(y) - 1

p1 <- ggplot(data=tibble(x=times,y=y)) + 
  geom_ribbon(aes(x=x,ymax=y,ymin=0),fill="grey70",size=0.5,col="black") +
  geom_segment(aes(x=midpoint-25,xend=midpoint-25,y=0,yend=y[midpoint-25+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_segment(aes(x=midpoint+25,xend=midpoint+25,y=0,yend=y[midpoint+25+1]),size=0.25,linetype="dotted",col="grey40") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.018),breaks=seq(0,0.018,by=0.005),labels=seq(0,0.018,by=0.005)*100000) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,350,by=50),limits=c(0,170)) +
  ylab("Incidence \nper 100,000") +
  xlab("Time") +
  theme_pubr() +
  theme(text=element_text(family="sans"))+
  theme(legend.position=c(0.7,0.45),
        axis.ticks = element_line(color="black"),
        panel.grid=element_blank(),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        axis.text=element_text(size=7,color="black"),
        axis.title=element_text(size=7,color="black"),
        legend.box.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.2,0.2,unit="cm"),
        legend.background =  element_rect(color=NA,fill=NA)) 
p1
p1_alt <- ggplot(data=tibble(x=times,y=y)) + 
  #geom_segment(aes(x=midpoint-25,xend=midpoint-25,y=0,yend=y[midpoint-25+1]),size=0.25,linetype="dotted",col="grey40") +
  #geom_segment(aes(x=midpoint+25,xend=midpoint+25,y=0,yend=y[midpoint+25+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_bar(aes(x=x,y=y),size=0.25,stat="identity",col="black") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.018),breaks=seq(0,0.018,by=0.005),labels=seq(0,0.018,by=0.005)*100000) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,350,by=50),limits=c(0,170)) +
  ylab("Incidence \nper 100,000") +
  xlab("Time") +
  theme_pubr() +
  theme(text=element_text(family="sans"))+
  theme(legend.position=c(0.7,0.45),
        panel.grid=element_blank(),
        axis.ticks = element_line(color="black"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        axis.text=element_text(size=7,color="black"),
        axis.title=element_text(size=7,color="black"),
        legend.box.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.2,0.2,unit="cm"),
        legend.background =  element_rect(color=NA,fill=NA)) 
p1_alt



scale_amount <- (1:(midpoint+50+1) - (midpoint/2))^14
scale_amount <- scale_amount/sum(scale_amount)
scale_amount <- (15*scale_amount + 0.1) + 5*exp(-0.1*(1:(midpoint+50+1)))
plot(scale_amount)
scale_amount_bot <- pmin(scale_amount,1)

plot_dat <- tibble(x=times,y=pmax(fit$y,0),ymax=pmax(fit$y*1.1+0.0005,0),ymin=pmax(fit$y*0.9-0.0005,0))
y1 <- fit$y
plot_dat1 <- plot_dat %>% filter(x <= midpoint+50) %>%
  mutate(ymax = pmax((y1[1:(midpoint+50+1)]+0.0005)*(1+scale_amount),0),
         ymin = pmax((y1[1:(midpoint+50+1)]-0.0005)*(1-scale_amount_bot),0))

p2 <- ggplot(data=plot_dat) + 
  geom_ribbon(data=plot_dat1,aes(x=x,ymin=ymin,ymax=ymax),alpha=0.25,fill="#56B4E9") + 
  geom_line(aes(x=x,y=y),size=0.5) +
  geom_segment(aes(x=midpoint,xend=midpoint,y=0,yend=y[midpoint+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_segment(aes(x=midpoint+25,xend=midpoint+25,y=0,yend=y[midpoint+25+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_segment(aes(x=midpoint+50,xend=midpoint+50,y=0,yend=y[midpoint+50+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_segment(aes(x=midpoint-25,xend=midpoint-25,y=0,yend=y[midpoint-25+1]),size=0.25,linetype="dotted",col="grey40") +
  geom_segment(aes(x=midpoint-50,xend=midpoint-50,y=0,yend=y[midpoint-50+1]),size=0.25,linetype="dotted",col="grey40") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.018),breaks=seq(0,0.018,by=0.005),labels=seq(0,0.018,by=0.005)*100000) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,350,by=50),limits=c(0,170)) +
  ylab("Estimated incidence") +
  xlab("Time") +
  theme_pubr() +
  theme(text=element_text(family="sans"))+
  theme(legend.position=c(0.7,0.45),
        axis.ticks = element_line(color="black"),
        panel.grid=element_blank(),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        axis.text=element_text(size=7,color="black"),
        axis.title=element_text(size=7,color="black"),
        legend.box.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.2,0.2,unit="cm"),
        legend.background =  element_rect(color=NA,fill=NA)) 
p2

if(save_plots){
  ggsave("figures/Figure0/fig0_baseA.pdf",p1,height=1.3,width=3.5,colormode="cmyk")
  #ggsave("figures/Figure0/fig0_baseA_alt.pdf",p1_alt,height=1.3,width=3.7,colormode="cmyk")
  ggsave("figures/Figure0/fig0_baseB.pdf",p2,height=1.3,width=3.5,colormode="cmyk")
}

complete_linelist <- virosolver::simulate_observations_wrapper(incidence=seir_dynamics$incidence,times=times,
                                                               population_n=population_n)
reporting_prob <- 1
observed_linelist <- simulate_reporting(complete_linelist %>% filter(is_infected==1), 
                                        frac_report=reporting_prob,
                                        timevarying_prob=NULL,
                                        solve_times=midpoint + seq(-50,50,by=25), 
                                        symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,kinetics_pars=pars)

plot_dat_cts1 <- simulated_viral_loads %>% filter(ct_obs < 40) %>% 
  filter(sampled_time %in% c(midpoint-25,midpoint+25)) %>%
  mutate(sampled_time=ifelse(sampled_time==midpoint-25,"Sample 1","Sample 2")) %>%
  #mutate(sampled_time=ifelse(sampled_time==75,"Sample 1","Sample 2")) %>%
  mutate(ct_range="Low",ct_range=ifelse(ct_obs < 35,"Medium",ct_range),ct_range=ifelse(ct_obs < 28,"High",ct_range)) %>%
  mutate(ct_range="Recent",ct_range=ifelse(days_since_infection > 7,"Intermediate",ct_range),ct_range=ifelse(days_since_infection > 14,"Old",ct_range))
plot_dat_cts1 %>% group_by(sampled_time) %>%
  summarize(med_ct=median(ct_obs),skew_ct=moments::skewness(ct_obs))

set.seed(4)
p_option1 <- plot_dat_cts1 %>%  
  group_by(sampled_time) %>%
  sample_n(50) %>%
  ungroup() %>%
  ggplot() + 
  geom_jitter(aes(x=sampled_time,y=ct_obs,col=ct_range),width=0.1,height=0,size=0.75) +
  geom_hline(yintercept=40,linetype="dotted") +
  #scale_color_manual(values=c("Low"="darkgreen","Medium"="orange","High"="red")) +
  scale_color_manual(values=c("Old"="#1E88E5","Intermediate"="#FFC107","Recent"="#D81B60")) +
  scale_y_continuous(expand=c(0,0),trans="reverse",breaks=seq(15,40,by=5),limits=c(42,12),
                     sec.axis=sec_axis(trans=~.*-1/log2(10) + 40/log2(10) + 3, name="Viral load",breaks=seq(3,11,by=1))) +
  ylab("Ct value") +
  theme_pubr() +
  theme(text=element_text(family="sans"))+
  theme(legend.position="none",
        axis.ticks=element_line(color="black"),
        panel.grid=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        axis.text=element_text(size=6,color="black"),
        axis.title.y=element_text(size=6,color="black"),
        axis.title.x=element_blank(),
        legend.box.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.2,0.2,unit="cm"),
        legend.background =  element_rect(color=NA,fill=NA)) 
p_option1


set.seed(1)
plot_dat_cts1 %>%  
  group_by(sampled_time) %>%
  sample_n(60) %>%
  summarize(med_ct=median(ct_obs),skew_ct=moments::skewness(ct_obs))
set.seed(1)
p_option1_alt <- plot_dat_cts1 %>%  
  group_by(sampled_time) %>%
  sample_n(60) %>%
  ungroup() %>%
  mutate(ct_obs = floor(ct_obs)) %>%
  group_by(ct_obs,sampled_time) %>% 
  arrange(ct_range) %>%
  mutate(top_range=n(),dist=(top_range-1)/2) %>%
  mutate(i=1:n()-1,i=i-dist)%>%
  ggplot() + 
  geom_point(aes(x=as.numeric(as.factor(sampled_time))+i/20,y=ct_obs,col=ct_range),size=0.3)+
  scale_x_continuous(limits=c(0.75,2.25),breaks=seq(1,2,by=1)) +
  geom_hline(yintercept=40,linetype="dotted") +
  scale_color_manual(values=c("Old"="#1E88E5","Intermediate"="#FFC107","Recent"="#D81B60")) +
  scale_y_continuous(expand=c(0,0),trans="reverse",breaks=seq(15,40,by=5),limits=c(42,12),
                     sec.axis=sec_axis(trans=~.*-1/log2(10) + 40/log2(10) + 3, name="Viral load",breaks=seq(3,11,by=2))) +
  ylab("Ct value") +
  theme_pubr() +
  theme(text=element_text(family="sans"))+
  theme(legend.position="none",
        axis.ticks=element_line(color="black"),
        panel.grid=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        axis.text=element_text(size=6,color="black"),
        axis.title.y=element_text(size=6,color="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.box.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.2,0.2,unit="cm"),
        legend.background =  element_rect(color=NA,fill=NA)) 
p_option1_alt

p_option1_hist <- plot_dat_cts1 %>%
  group_by(sampled_time) %>%
  ungroup() %>%
  ggplot() + 
  geom_histogram(aes(x=ct_obs,y=..density..),binwidth=2,boundary=0,closed="left",
                 fill="grey70",col="black",size=0.25) +
  facet_wrap(~sampled_time) +
  scale_color_manual(values=c("Old"="#26d800","Intermediate"="#ffa626","Recent"="#ff0000")) +
  scale_x_continuous(trans="reverse",breaks=seq(10,40,by=5),limits=c(40,10)) +
  ylab("Density") +
  xlab("Ct value") +
  theme_void() +
  theme(strip.text=element_blank())
p_option1_hist


p_option1_hist_all <- simulated_viral_loads %>% filter(ct_obs < 40) %>%
  ggplot() + 
  geom_histogram(aes(x=ct_obs,y=..density..),binwidth=2,boundary=0,closed="left",
                 fill="grey70",col="black",size=0.25) +
  facet_wrap(~sampled_time,nrow=1) +
  scale_color_manual(values=c("Old"="#26d800","Intermediate"="#ffa626","Recent"="#ff0000")) +
  scale_x_continuous(trans="reverse",breaks=seq(5,40,by=5),limits=c(40,5)) +
  xlab("Ct value") +
  theme_void() +
  theme(strip.text=element_blank())
p_option1_hist_all

if(save_plots){
  ggsave("figures/Figure0/ct_dotplots.pdf",p_option1,height=1.2,width=3,colormode="cmyk")
  ggsave("figures/Figure0/ct_dotplots_alt.pdf",p_option1_alt,height=1.2,width=3,colormode="cmyk")
  ggsave("figures/Figure0/subset_hists.pdf",p_option1_hist,height=1.5,width=3.5,colormode="cmyk")
  ggsave("figures/Figure0/all_hists.pdf",p_option1_hist_all,height=1.25,width=8,colormode="cmyk")
}

