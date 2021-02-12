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

setwd("~/Documents/GitHub/virosolver_paper")

cumu_data <- FALSE
use_pos <- TRUE
i <- 7
min_date <- as.Date("2020-03-01")
max_age <- NA
## Testing SEIR model
parTab <- read.csv("pars/nursing_homes/partab_seir_model.csv")
pars <- parTab$values
names(pars) <- parTab$names
means <- pars

## Test with nursing home data
dat_wide <- read_csv(file = "data/nursing_home_data.csv")
dat_wide$collection_date <- as.Date(dat_wide$collection_date,format="%m/%d/%y")
dat_wide <- dat_wide %>% left_join(
  dat_wide %>% filter(is_nursinghome==1) %>% 
    group_by(location, week) %>% 
    count(collection_date) %>% 
    filter(n == max(n)) %>% 
    rename(mean_week_date=collection_date) %>% 
    ungroup() %>% 
    dplyr::select(-n))
## 3rd sampling time for nursing home 3 is not random collection
dat_wide <- dat_wide %>% filter(!(location == "NH 4" & mean_week_date == "2020-04-29"))
dat_use <- dat_wide %>% filter(is_nursinghome==1) %>%
  dplyr::select(location, mean_week_date, N2) %>%
  rename(ct=N2)


unique_data_combs <- dat_use %>% dplyr::select(location, mean_week_date) %>% distinct()

## If using cumulative cross sections or independent
if(cumu_data) {
  dat_tmp <- dat_use %>%
    filter(mean_week_date <= unique_data_combs$mean_week_date[i]) %>%
    filter(location == unique_data_combs$location[i]) %>%
    mutate(t=as.numeric(mean_week_date - min_date)) %>%
    mutate(ct=ifelse(is.na(ct),40,ct))
} else {
  dat_tmp <- dat_use %>%
    filter(mean_week_date == unique_data_combs$mean_week_date[i]) %>%
    filter(location == unique_data_combs$location[i]) %>%
    mutate(t=as.numeric(mean_week_date - min_date)) %>%
    mutate(ct=ifelse(is.na(ct),40,ct))
}

## If only using positive Cts
if(use_pos) {
  dat_tmp <- dat_tmp %>% filter(ct < 40)
}

dat_tmp_used_in_run <- dat_tmp

## Observation times
if(!is.na(max_age)){
  dat_tmp_used_in_run <- dat_tmp_used_in_run %>% mutate(t = t - max(t),
                                                        t = t + max_age)
}
ages <- 1:max(dat_tmp_used_in_run$t)
times <- 0:max(dat_tmp_used_in_run$t)

## Epidemic cannot start after first observation time
parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(dat_tmp_used_in_run$t) -7 

inc_func_use <- solveSEIRModel_rlsoda_wrapper

## Test that model can be run without failure
f <- create_posterior_func(parTab, dat_tmp_used_in_run, CREATE_POSTERIOR_FUNC=create_posterior_func,
                           INCIDENCE_FUNC=inc_func_use,
                           PRIOR_FUNC=NULL,
                           use_pos=use_pos,
                           t_dist=NULL,solve_ver = "model")



chain <- read.csv("~/Documents/GitHub/virosolver_paper/mcmc_chains/1.nursing_home_ct/seir_TRUE_FALSE/ct_NH 3_2020-04-09/ct_NH 3_2020-04-09_chainno_1_chain.csv")
pars <- get_best_pars(chain)
all_solved <- NULL
inc_all <- NULL
R0s <- seq(1,25,by=0.5)
pars["desired_mode"] <- 2
pars["intercept"] <- 40
pars["t_switch"] <- 5
pars["prob_detect"] <- 0.8
pars["viral_peak"] <- 15
pars["wane_rate2"] <- 5

plot(viral_load_func(pars, seq(0,35,by=1),convert_vl = FALSE))

for(j in seq_along(R0s)){
  R0 <- R0s[j]
  pars1 <- pars
  pars1["t0"] <- 0
  pars1["R0"] <- R0
  obs <- f(pars1)
  obs$R0 <- R0
  all_solved[[j]] <- obs
  inc_all[[j]] <- tibble(t=times,y=inc_func_use(pars1,times))
}

all_solved <- do.call("bind_rows",all_solved)

obs_detect <- all_solved %>% filter(ct < pars["intercept"]) %>% group_by(R0) %>% mutate(density1=density/sum(density))

p1 <- ggplot(obs_detect) + geom_histogram(data=dat_tmp_used_in_run, aes(x=ct,y=..density..)) + geom_line(aes(x=ct,y=density1,col=R0,group=R0)) + scale_color_viridis_c()
p2 <- ggplot(all_solved %>% filter(ct == pars["intercept"])) + geom_bar(aes(x=R0,y=1-density),stat="identity") + scale_y_continuous(limits=c(0,1))
p1/p2
