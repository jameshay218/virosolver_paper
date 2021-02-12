########################################
## SCRIPT TO COMBINE ALL NURSING HOME SIMULATION ANALYSES
## James Hay 1 December 2020
##  - Combines and the diagnostic table from all simulated runs
##  - Reads in MCMC chains from each run
##  - Population size analyses: Plots the posterior median growth rates, 95% CrI widths and posterior prob gr > 0
##  - Prior strength analyses: Plots the posterior median growth rates, 95% CrI widths and posterior prob gr > 0
##  - Main sim-recovery: Plots the posterior median growth rate estimates for the exponential and SEIR models, 
##                       with or without negative Ct values, against the simulated stochastic SEIR trajectory


########################################
## 1. Headers
########################################
library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
library(coda)

HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

source(paste0(HOME_WD,"/virosolver_paper/code/priors.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/plot_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/linelist_sim_funcs.R"))
source(paste0(HOME_WD,"/virosolver_paper/code/odin_funcs.R"))


mcmcPars_ct <- c("iterations"=100000,"popt"=0.234,"opt_freq"=1000,
                 "thin"=1,"adaptive_period"=100000,"save_block"=1000)


## Arguments for this run. Note that this will automatically find when Rt is 1 (for the peak) and fit to 20 days either side
control_table <- read_csv(paste0(HOME_WD,"/virosolver_paper/pars/nursing_homes/sim_control_nh.csv"))

## Only need first run index for each analysis
control_table <- control_table %>% filter(run_index == 1)

########################################
## 2. Compile population size analysis
########################################
## First, look at the population size analyses
control_table_pop <- control_table %>% filter(run_name=="sim_NH_pop")
pop_sizes <- unique(control_table_pop$population_n)

diagnostics_all <- NULL
gr_samples_all <- NULL
chains_all <- NULL

## For each row in control table
for(simno in 1:nrow(control_table_pop)){
  setwd(HOME_WD)
  ## Names and parameters for this run
  runname <- control_table_pop$run_name[simno] 
  use_pos <- control_table_pop$use_pos[simno] 
  model_version <- control_table_pop$model_version[simno] 
  top_chainwd <- control_table_pop$chainwd[simno] 
  top_plotwd <- control_table_pop$plotwd[simno] 
  prior_func_use <- get(control_table_pop$prior_func[simno])
  inc_func_use <- get(control_table_pop$inc_func[simno])
  age_max <- control_table_pop$age_max[simno]
  population_n <- control_table_pop$population_n[simno]
  prior_sd_mod <- control_table_pop$prior_mod[simno]
  
  ## Parameter table
  parTab <- read.csv(control_table$parTab_file[simno], stringsAsFactors=FALSE)
  
  ## Get full filename and working directories
  full_runname <- paste0(runname,"_",model_version,"_pos",use_pos,"_n", population_n,"_prior",prior_sd_mod,"_run")
  chainwd <- paste0(top_chainwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  plot_wd <- paste0(top_plotwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  
  ## Read in diagnostics
  setwd(plot_wd)
  all_runs <- list.files()
  diagnostics <- NULL
  gr_samples <- NULL
  for(index in seq_along(all_runs)){
    setwd(plot_wd)
    setwd(all_runs[index])
    if(file.exists("diagnostics.csv")){
      diagnostics[[index]] <- read_csv("diagnostics.csv")
    }
    if(file.exists("gr_samples.csv")){
      gr_samples[[index]] <- read_csv("gr_samples.csv")
    }
  }
  diagnostics_tmp <- do.call("bind_rows", diagnostics)
  diagnostics_tmp$pop_n <- population_n
  diagnostics_all[[simno]] <- diagnostics_tmp
  
  gr_samples_tmp <- do.call("bind_rows", gr_samples)
  gr_samples_tmp$pop_n <- population_n
  gr_samples_all[[simno]] <- gr_samples_tmp
  
  ## Read in chains
  tmp_chains_list <- NULL
  
  setwd(chainwd)
  all_runs <- list.files()
  for(index in seq_along(all_runs)){
    setwd(chainwd)
    setwd(all_runs[index])
    
    tmp_chains_list2 <- NULL
    phases <- list.files()
    for(phase in phases){
      ## Growth chains
      if(length(list.files(phase) > 0)){
        tmp_chains <- load_mcmc_chains(phase, parTab,FALSE,10,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=FALSE)
        tmp_chain <- as.data.frame(tmp_chains$chain)
        tmp_chain$sampno <- 1:nrow(tmp_chain)
        tmp_chain <- tmp_chain %>% 
          pivot_longer(-sampno) %>% 
          rename(par=name) %>% 
          mutate(phase=phase, pop_n=population_n, index=index)
        tmp_chains_list2[[phase]] <- tmp_chain
      } else {
        tmp_chains_list2[[phase]] <- NULL
      }
    }
    if(is.list(tmp_chains_list2)){
      tmp_chains_list[[index]] <- do.call("bind_rows",tmp_chains_list2)
    }
  }
  chains_all[[simno]] <- do.call("bind_rows", tmp_chains_list)
}

diagnostics_all_comb_pop <- do.call("bind_rows", diagnostics_all)
all_chains_comb_pop <- do.call("bind_rows",chains_all)
########################
## NOTE HERE - A COUPLE SIMULATIONS WENT WHACKY. TRUNCATING ESTIMATES FOR PLOT
########################
gr_samples_all_comb_pop <- do.call("bind_rows", gr_samples_all) %>% filter(daily_gr < 1 & daily_gr > -1)
use_pars <- c("viral_peak","level_switch","obs_sd")

########################

########################
## MANUSCRIPT PLOTS
########################
pop_dat <- gr_samples_all_comb_pop %>% group_by(label, index, pop_n) %>%
  summarize(post_med=median(daily_gr),
            lower=quantile(daily_gr,0.025),
            upper=quantile(daily_gr,0.975),
            width=upper-lower) 

label_key <- c("growth"="Growth phase","decline"="Decline phase","peak"="Peak")
pop_dat$label <- label_key[as.character(pop_dat$label)]
pop_dat$label <- factor(pop_dat$label, levels=c("Growth phase","Peak","Decline phase"))

pop_growth_dat <- gr_samples_all_comb_pop %>% group_by(label, index, pop_n) %>%
  mutate(is_growth=daily_gr > 0) %>%
  summarize(prob_growth=sum(is_growth)/n()) 
pop_growth_dat$label <- label_key[as.character(pop_growth_dat$label)]
pop_growth_dat$label <- factor(pop_growth_dat$label, levels=c("Growth phase","Peak","Decline phase"))

p_pop_medians <-ggplot(pop_dat) + 
  geom_violin(aes(x=as.factor(pop_n),y=post_med),draw_quantiles = c(0.025,0.5,0.975),fill="grey70",scale="width") +
  geom_jitter(aes(x=as.factor(pop_n),y=post_med),height=0,size=0.1,width=0.15) +
  geom_hline(yintercept=0,linetype="dashed") +
  facet_wrap(~label)+
  xlab("") +
  ylab("Latest growth rate (posterior median)") +
  export_theme +
  labs(tag="A")

p_pop_widths <- ggplot(pop_dat) + 
  geom_violin(aes(x=as.factor(pop_n),y=width),draw_quantiles = c(0.025,0.5,0.975),fill="grey70",scale="width") +
  geom_jitter(aes(x=as.factor(pop_n),y=width),height=0,size=0.1,width=0.15) +
  facet_wrap(~label) +
  ylab("Latest growth rate 95% CrI width") +
  xlab("") +
  export_theme +
  labs(tag="B")


p_pop_growth <-ggplot(pop_growth_dat) + 
  geom_violin(aes(x=as.factor(pop_n),y=prob_growth),draw_quantiles = c(0.025,0.5,0.975),scale="width",fill="grey70") +
  geom_jitter(aes(x=as.factor(pop_n),y=prob_growth),height=0,size=0.1,width=0.15) +
  facet_wrap(~label)+
  geom_hline(yintercept=0.5,linetype="dashed") +
  ylab("Posterior probability in growth phase") +
  xlab("Simulated population size") +
  export_theme +
  labs(tag="C")

p_pop <- p_pop_medians/p_pop_widths/p_pop_growth
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_pop.pdf"),p_pop,width=8,height=8)
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_pop.png"),p_pop,width=8,height=8,units="in",dpi=300)
########################################
## 3. Compile prior power analysis
########################################
## First, look at the population size analyses
control_table_prior <- control_table %>% filter(run_name=="sim_NH_prior")

diagnostics_all <- NULL
gr_samples_all <- NULL
chains_all <- NULL

## For each row in control table
for(simno in 1:nrow(control_table_prior)){
  setwd(HOME_WD)
  ## Names and parameters for this run
  runname <- control_table_prior$run_name[simno] 
  use_pos <- control_table_prior$use_pos[simno] 
  model_version <- control_table_prior$model_version[simno] 
  top_chainwd <- control_table_prior$chainwd[simno] 
  top_plotwd <- control_table_prior$plotwd[simno] 
  prior_func_use <- get(control_table_prior$prior_func[simno])
  inc_func_use <- get(control_table_prior$inc_func[simno])
  age_max <- control_table_prior$age_max[simno]
  population_n <- control_table_prior$population_n[simno]
  prior_sd_mod <- control_table_prior$prior_mod[simno]
  
  ## Parameter table
  parTab <- read.csv(control_table$parTab_file[simno], stringsAsFactors=FALSE)
  
  ## Get full filename and working directories
  full_runname <- paste0(runname,"_",model_version,"_pos",use_pos,"_n", population_n,"_prior",prior_sd_mod,"_run")
  chainwd <- paste0(top_chainwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  plot_wd <- paste0(top_plotwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  
  ## Read in diagnostics
  setwd(plot_wd)
  all_runs <- list.files()
  diagnostics <- NULL
  gr_samples <- NULL
  for(index in seq_along(all_runs)){
    setwd(plot_wd)
    setwd(all_runs[index])
    if(file.exists("diagnostics.csv")){
      diagnostics[[index]] <- read_csv("diagnostics.csv")
    }
    if(file.exists("gr_samples.csv")){
      gr_samples[[index]] <- read_csv("gr_samples.csv")
    }
  }
  diagnostics_tmp <- do.call("bind_rows", diagnostics)
  diagnostics_tmp$prior_strength <- prior_sd_mod
  diagnostics_all[[simno]] <- diagnostics_tmp
  
  gr_samples_tmp <- do.call("bind_rows", gr_samples)
  gr_samples_tmp$prior_strength <- prior_sd_mod
  gr_samples_all[[simno]] <- gr_samples_tmp
  
  ## Read in chains
  tmp_chains_list <- NULL
  
  setwd(chainwd)
  all_runs <- list.files()
  for(index in seq_along(all_runs)){
    setwd(chainwd)
    setwd(all_runs[index])
    
    tmp_chains_list2 <- NULL
    phases <- list.files()
    for(phase in phases){
      ## Growth chains
      if(length(list.files(phase) > 0)){
      tmp_chains <- load_mcmc_chains(phase, parTab,FALSE,10,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=FALSE)
      tmp_chain <- as.data.frame(tmp_chains$chain)
      tmp_chain$sampno <- 1:nrow(tmp_chain)
      tmp_chain <- tmp_chain %>% 
        pivot_longer(-sampno) %>% 
        rename(par=name) %>% 
        mutate(phase=phase, prior_strength=prior_sd_mod, index=index)
      tmp_chains_list2[[phase]] <- tmp_chain
      } else {
        tmp_chains_list2[[phase]] <- NULL
      }
    }
    if(is.list(tmp_chains_list2)){
      tmp_chains_list[[index]] <- do.call("bind_rows",tmp_chains_list2)
    }
  }
  chains_all[[simno]] <- do.call("bind_rows", tmp_chains_list)
}

diagnostics_all_comb_prior <- do.call("bind_rows", diagnostics_all)
all_chains_comb_prior <- do.call("bind_rows",chains_all)
gr_samples_all_comb_prior <- do.call("bind_rows", gr_samples_all) %>% filter(daily_gr < 1 & daily_gr > -1)
use_pars <- c("viral_peak","level_switch","obs_sd")

########################
## MANUSCRIPT PLOTS
########################
pop_dat_prior <- gr_samples_all_comb_prior %>% group_by(label, index, prior_strength) %>%
  summarize(post_med=median(daily_gr),
            lower=quantile(daily_gr,0.025),
            upper=quantile(daily_gr,0.975),
            width=upper-lower) 

label_key <- c("growth"="Growth phase","decline"="Decline phase","peak"="Peak")
pop_dat_prior$label <- label_key[as.character(pop_dat_prior$label)]
pop_dat_prior$label <- factor(pop_dat_prior$label, levels=c("Growth phase","Peak","Decline phase"))

pop_growth_dat_prior <- gr_samples_all_comb_prior %>% group_by(label, index, prior_strength) %>%
  mutate(is_growth=daily_gr > 0) %>%
  summarize(prob_growth=sum(is_growth)/n()) 
pop_growth_dat_prior$label <- label_key[as.character(pop_growth_dat_prior$label)]
pop_growth_dat_prior$label <- factor(pop_growth_dat_prior$label, levels=c("Growth phase","Peak","Decline phase"))

p_prior_medians <-ggplot(pop_dat_prior) + 
  geom_violin(aes(x=as.factor(prior_strength),y=post_med),draw_quantiles = c(0.025,0.5,0.975),fill="grey70",scale="width") +
  geom_jitter(aes(x=as.factor(prior_strength),y=post_med),height=0,size=0.1,width=0.15) +
  geom_hline(yintercept=0,linetype="dashed") +
  facet_wrap(~label)+
  xlab("") +
  ylab("Latest growth rate (posterior median)") +
  export_theme +
  labs(tag="A")

p_prior_widths <- ggplot(pop_dat_prior) + 
  geom_violin(aes(x=as.factor(prior_strength),y=width),draw_quantiles = c(0.025,0.5,0.975),fill="grey70",scale="width") +
  geom_jitter(aes(x=as.factor(prior_strength),y=width),height=0,size=0.1,width=0.15) +
  facet_wrap(~label) +
  ylab("Latest growth rate 95% CrI width") +
  xlab("") +
  export_theme +
  labs(tag="B")

p_prior_growth <-ggplot(pop_growth_dat_prior) + 
  geom_violin(aes(x=as.factor(prior_strength),y=prob_growth),draw_quantiles = c(0.025,0.5,0.975),scale="width",fill="grey70") +
  geom_jitter(aes(x=as.factor(prior_strength),y=prob_growth),height=0,size=0.1,width=0.15) +
  facet_wrap(~label)+
  geom_hline(yintercept=0.5,linetype="dashed") +
  ylab("Posterior probability in growth phase") +
  xlab("Relative prior strength") +
  export_theme +
  labs(tag="C")

p_pop_prior <- p_prior_medians/p_prior_widths/p_prior_growth
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_prior.pdf"),p_pop_prior,width=8,height=8)
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_prior.png"),p_pop_prior,width=8,height=8,units="in",dpi=300)


########################################
## 3. Compile main analyses
########################################
## First, look at the population size analyses
control_table_main <- control_table %>% filter(run_name=="sim_NH_normal")
pop_sizes <- unique(control_table_pop$population_n)

diagnostics_all <- NULL
gr_samples_all <- NULL
chains_all <- NULL
linelists_all <- NULL
for(simno in 1:nrow(control_table_main)){
  setwd(HOME_WD)
  ## Names and parameters for this run
  runname <- control_table_main$run_name[simno] 
  use_pos <- control_table_main$use_pos[simno] 
  model_version <- control_table_main$model_version[simno] 
  top_chainwd <- control_table_main$chainwd[simno] 
  top_plotwd <- control_table_main$plotwd[simno] 
  prior_func_use <- get(control_table_main$prior_func[simno])
  inc_func_use <- get(control_table_main$inc_func[simno])
  age_max <- control_table_main$age_max[simno]
  population_n <- control_table_main$population_n[simno]
  prior_sd_mod <- control_table_main$prior_mod[simno]
  
  ## Parameter table
  parTab <- read.csv(control_table_main$parTab_file[simno], stringsAsFactors=FALSE)
  
  ## Get full filename and working directories
  full_runname <- paste0(runname,"_",model_version,"_pos",use_pos,"_n", population_n,"_prior",prior_sd_mod,"_run")
  chainwd <- paste0(top_chainwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  plot_wd <- paste0(top_plotwd, "/",runname,"/",model_version,"_pos",use_pos,"_n", population_n,"/prior_strength",prior_sd_mod)
  
  ## Read in diagnostics
  setwd(plot_wd)
  all_runs <- list.files()
  diagnostics <- NULL
  gr_samples <- NULL
  linelists <- NULL
  for(index in seq_along(all_runs)){
    setwd(plot_wd)
    setwd(all_runs[index])
    if(file.exists("diagnostics.csv")){
      diagnostics[[index]] <- read_csv("diagnostics.csv")
    }
    if(file.exists("gr_samples.csv")){
      gr_samples[[index]] <- read_csv("gr_samples.csv")
    }
    if(file.exists("full_linelist.csv")){
      linelists[[index]] <- read_csv("full_linelist.csv") %>% mutate(index = index)
    }
    
  }
  diagnostics_tmp <- do.call("bind_rows", diagnostics)
  diagnostics_tmp$use_pos <- use_pos
  diagnostics_tmp$model_version <- model_version
  diagnostics_all[[simno]] <- diagnostics_tmp
  
  gr_samples_tmp <- do.call("bind_rows", gr_samples)
  gr_samples_tmp$use_pos <- use_pos
  gr_samples_tmp$model_version <- model_version
  gr_samples_all[[simno]] <- gr_samples_tmp
  
  linelists_tmp <- do.call("bind_rows", linelists)
  linelists_tmp$use_pos <- use_pos
  linelists_tmp$model_version <- model_version
  linelists_all[[simno]] <- linelists_tmp
  
  ## Read in chains
  tmp_chains_list <- NULL
  
  setwd(chainwd)
  all_runs <- list.files()
  for(index in seq_along(all_runs)){
    setwd(chainwd)
    setwd(all_runs[index])
    
    tmp_chains_list2 <- NULL
    phases <- list.files()
    for(phase in phases){
      ## Growth chains
      if(length(list.files(phase) > 0)){
        tmp_chains <- load_mcmc_chains(phase, parTab,FALSE,10,mcmcPars_ct["adaptive_period"],multi=TRUE,chainNo=FALSE)
        tmp_chain <- as.data.frame(tmp_chains$chain)
        tmp_chain$sampno <- 1:nrow(tmp_chain)
        tmp_chain <- tmp_chain %>% 
          pivot_longer(-sampno) %>% 
          rename(par=name) %>% 
          mutate(phase=phase, use_pos=use_pos, model_version=model_version,index=index)
        tmp_chains_list2[[phase]] <- tmp_chain
      } else {
        tmp_chains_list2[[phase]] <- NULL
      }
    }
    if(is.list(tmp_chains_list2)){
      tmp_chains_list[[index]] <- do.call("bind_rows",tmp_chains_list2)
    }
  }
  chains_all[[simno]] <- do.call("bind_rows", tmp_chains_list)
}

diagnostics_all_comb_normal <- do.call("bind_rows", diagnostics_all)
all_chains_comb_normal <- do.call("bind_rows",chains_all)
gr_samples_all_comb_normal <- do.call("bind_rows", gr_samples_all) %>% filter(daily_gr < 1 & daily_gr > -1)
linelists_all_comb <- do.call("bind_rows", linelists_all)


########################
## MANUSCRIPT PLOTS
########################
use_pos_key <- c("TRUE"="Detectable Cts only","FALSE"="All Ct values")
model_version_key <- c("exp"="Exponential growth","seir"="SEIR model")
gr_samples_all_comb_normal$use_pos <- use_pos_key[as.character(gr_samples_all_comb_normal$use_pos)]
gr_samples_all_comb_normal$model_version <- model_version_key[as.character(gr_samples_all_comb_normal$model_version)]
linelists_all_comb$use_pos <- use_pos_key[as.character(linelists_all_comb$use_pos)]
linelists_all_comb$model_version <- model_version_key[as.character(linelists_all_comb$model_version)]

## Get peak times
peak_time <- times[which.max(seir_dynamics$incidence)] #seir_dynamics$seir_outputs$step[which.min(abs(1-seir_dynamics$seir_outputs$Rt))]
growth_time <- peak_time - 21
decline_time <- peak_time + 21

ll_counts <- linelists_all_comb %>% group_by(index, use_pos, model_version, infection_time) %>% tally() %>% drop_na()
ll_counts$label <- paste0(ll_counts$model_version, "; ",ll_counts$use_pos)
ll_counts <- ll_counts %>% group_by(index,use_pos,model_version, label) %>% complete(infection_time=1:125,fill=list(n=0))
ll_summaries <- ll_counts %>% group_by(use_pos, model_version, infection_time) %>% summarize(med=median(n),
                                                                                             lower=quantile(n,0.025),
                                                                                             upper=quantile(n,0.975))
ll_summaries$label <- paste0(ll_summaries$model_version, "; ",ll_summaries$use_pos)

peak_times <- ll_counts %>% group_by(index,use_pos,model_version,label) %>% filter(n==max(n)) %>%
  rename(peak_time=infection_time) %>%
  mutate(growth_phase=peak_time-21,
         decline_phase=peak_time+21)
peak_time_summaries <- peak_times %>% group_by(use_pos,model_version,label) %>%
  summarize(median_peak=median(peak_time),
            median_growth=median(growth_phase),
            median_decline=median(decline_phase),
            lower_peak=quantile(peak_time,c(0.025)),
            lower_growth=quantile(growth_phase,c(0.025)),
            lower_decline=quantile(decline_phase,c(0.025)),
            upper_peak=quantile(peak_time,c(0.975)),
            upper_growth=quantile(growth_phase,c(0.975)),
            upper_decline=quantile(decline_phase,c(0.975)))

p_dat <- ggplot(ll_counts %>% filter(index <= 25)) + 
  geom_line(aes(x=infection_time,y=n,group=index),size=0.1) + 
  geom_line(data=ll_summaries,aes(x=infection_time,y=med,col="Median")) +
  geom_pointrange(data=peak_time_summaries,aes(xmin=lower_growth,xmax=upper_growth,x=median_growth,col="Growth phase"),y=28,size=0.25)+
  geom_pointrange(data=peak_time_summaries,aes(xmin=lower_peak,xmax=upper_peak,x=median_peak,col="Peak"),y=30,size=0.25)+
  geom_pointrange(data=peak_time_summaries,aes(xmin=lower_decline,xmax=upper_decline,x=median_decline,col="Decline phase"),y=32,size=0.25)+
  scale_color_manual(values=c("Growth phase"="#EE0000FF","Peak"="#3B4992FF","Decline phase"="#008B45FF","Median"="#631879FF"),
                     breaks=c("Growth phase","Peak","Decline phase","Median"))+
  guides(color=guide_legend()) +
  facet_wrap(~label,ncol=1) +
  export_theme +
  coord_cartesian(xlim=c(0,90)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,35)) +
  theme(legend.position="bottom",legend.title=element_blank()) +
  ylab("New infections") +
  xlab("Day") +
  labs(tag="A")


pop_dat_normal <- gr_samples_all_comb_normal %>% group_by(label, index, use_pos,model_version) %>%
  summarize(post_med=median(daily_gr),
            lower=quantile(daily_gr,0.025),
            upper=quantile(daily_gr,0.975),
            width=upper-lower) 

label_key <- c("growth"="Growth phase","decline"="Decline phase","peak"="Peak")
pop_dat_normal$label <- label_key[as.character(pop_dat_normal$label)]
pop_dat_normal$label <- factor(pop_dat_normal$label, levels=c("Growth phase","Peak","Decline phase"))

pop_dat_normal <- pop_dat_normal %>% rename(Data=use_pos)

pop_growth_dat_normal <- gr_samples_all_comb_normal %>% group_by(label, index, use_pos,model_version) %>%
  mutate(is_growth=daily_gr > 0) %>%
  summarize(prob_growth=sum(is_growth)/n()) 
pop_growth_dat_normal$label <- label_key[as.character(pop_growth_dat_normal$label)]
pop_growth_dat_normal$label <- factor(pop_growth_dat_normal$label, levels=c("Growth phase","Peak","Decline phase"))

pop_growth_dat_normal <- pop_growth_dat_normal %>% rename(Data=use_pos)

p_normal_medians <-ggplot(pop_dat_normal) + 
  geom_violin(aes(x=as.factor(model_version),y=post_med,fill=Data),draw_quantiles = c(0.025,0.5,0.975),scale="width",alpha=0.25) +
  geom_jitter(aes(x=as.factor(model_version),y=post_med,group=Data,col=Data),size=0.1,position=position_jitterdodge()) +
  geom_hline(yintercept=0,linetype="dashed") +
  ggsci::scale_color_aaas() +
  ggsci::scale_fill_aaas() +
  facet_wrap(~label)+
  xlab("") +
  ylab("Latest growth rate (posterior median)") +
  export_theme +
  theme(legend.position="none") +
  labs(tag="B")

p_normal_widths <- ggplot(pop_dat_normal) + 
  geom_violin(aes(x=as.factor(model_version),y=width,fill=Data),draw_quantiles = c(0.025,0.5,0.975),scale="width",alpha=0.25) +
  geom_jitter(aes(x=as.factor(model_version),y=width,group=Data,col=Data),size=0.1,position=position_jitterdodge()) +
  facet_wrap(~label) +
  ggsci::scale_color_aaas() +
  ggsci::scale_fill_aaas() +
  ylab("Latest growth rate 95% CrI width") +
  xlab("") +
  export_theme +
  theme(legend.position="none") +
  labs(tag="C")

p_normal_growth <-ggplot(pop_growth_dat_normal) + 
  geom_violin(aes(x=as.factor(model_version),y=prob_growth,fill=Data),draw_quantiles = c(0.025,0.5,0.975),scale="width",alpha=0.25) +
  geom_jitter(aes(x=as.factor(model_version),y=prob_growth,group=Data,col=Data),size=0.1,position=position_jitterdodge()) +
  ggsci::scale_color_aaas() +
  geom_hline(yintercept=0.5,linetype="dashed") +
  ggsci::scale_fill_aaas() +
  facet_wrap(~label)+
  ylab("Posterior probability in growth phase") +
  xlab("Model version") +
  export_theme +
  theme(legend.position="bottom") +
  labs(tag="D")

p_pop_normal <- p_dat + (p_normal_medians/p_normal_widths/p_normal_growth) + plot_layout(widths=c(1,2))
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_normal.pdf"),p_pop_normal,width=10,height=8)
ggsave(paste0(HOME_WD,"virosolver_paper/figures/supplement/nh_sim_normal.png"),p_pop_normal,width=10,height=8,units="in",dpi=300)
