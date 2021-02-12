library(lazymcmc)
library(tidyverse)
library(extraDistr)
library(ggthemes)
devtools::load_all("~/Documents/GitHub/virosolver/")
source("~/Documents/GitHub/virosolver_paper/code/priors.R")
source("~/Documents/GitHub/virosolver_paper/code/plot_funcs.R")
#setwd("~/Documents/GitHub/virosolver_paper/mcmc_chains/4.real_ma_ct/real_ma_free_gp_excl_last_week/")
#setwd("~/Downloads/1/")
setwd("~/Documents/GitHub/virosolver_paper/mcmc_chains/5.real_ma_single_timepoint/real_ma_seir_longer/")


parTab<- read.csv("~/Documents/GitHub/virosolver_paper/pars/partab_seir_model.csv")

means <- parTab$values
names(means) <- parTab$names

PTchain <- TRUE
multiChain <- TRUE
top_wd <- getwd()

n_samp <- 1000

timepoints <- list.files()

chains_all_list <- NULL
vl_traj_all <- NULL
prop_detect_all <- NULL
prior_values_all <- NULL
likelihood_values_all <- NULL

for(timepoint in timepoints){
  setwd(top_wd)
  setwd(timepoint)
  chains <- load_mcmc_chains(getwd(),parTab,FALSE,1,30000,multiChain,FALSE,PTchain)
  
  chain <- as.data.frame(chains$chain)
  chain$sampno <- 1:nrow(chain)
  samps <- sample(unique(chain$sampno),n_samp)
  ages <- 1:100
  vl_trajs <- matrix(0, nrow=length(samps),ncol=length(ages))
  trajs <- matrix(0, nrow=length(samps),ncol=length(ages))
  prior_vals <- numeric(length(samps))
  likelihood_vals <- numeric(length(samps))
  for(i in seq_along(samps)) {
    tmp_pars <- get_index_pars(chain, samps[i])
    prior_vals[i] <- prior_func_hinge_seir(tmp_pars)
    likelihood_vals[i] <- chain[chain$sampno == samps[i],"lnlike"] - prior_vals[i]
    tmp_pars["t_unit"] <- 1
    vl <- viral_load_func(tmp_pars,ages,TRUE)
    vl_trajs[i,]  <- rgumbel(length(vl),vl,tmp_pars["obs_sd"])
    trajs[i,] <- sapply(ages, function(a) prop_detectable(a, tmp_pars,vl))
  }
  
  
  trajs1 <- t(apply(trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
  trajs1 <- as.data.frame(trajs1)
  colnames(trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
  trajs1$t <- ages
  
  vl_trajs1 <- t(apply(vl_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
  vl_trajs1 <- as.data.frame(vl_trajs1)
  colnames(vl_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
  vl_trajs1$t <- ages
  
  trajs1$timepoint <- timepoint
  vl_trajs1$timepoint <- timepoint
  prior_vals <- tibble(prior=prior_vals,timepoint=timepoint)
  likelihood_vals <- tibble(lnlike=likelihood_vals,timepoint=timepoint)
  
  vl_traj_all[[timepoint]] <- vl_trajs1
  prop_detect_all[[timepoint]] <- trajs1
  prior_values_all[[timepoint]] <- prior_vals
  likelihood_values_all[[timepoint]] <- likelihood_vals
  
  chains <- load_mcmc_chains(getwd(),parTab,FALSE,1,30000,multiChain,FALSE,PTchain)
  chain <- as.data.frame(chains$chain)
  chain$sampno <- 1:nrow(chain)
  chain <- chain[,colnames(chain) != "prob"]
  chain_melted <- reshape2::melt(chain,id.vars="sampno")
  chain_melted$timepoint <- timepoint
  chains_all_list[[timepoint]] <- chain_melted
  
}

chains_all <- do.call("bind_rows", chains_all_list)
trajs1_all <- do.call("bind_rows",prop_detect_all)
prior_values_all_comb <- do.call("bind_rows",prior_values_all)
likelihood_values_all_comb <- do.call("bind_rows",likelihood_values_all)
vl_traj_comb <- do.call("bind_rows",vl_traj_all)

a <- prior_values_all_comb %>% mutate(ver="Prior probability") %>%
  rename(value=prior)
b <- likelihood_values_all_comb %>% mutate(ver="Log likelihood") %>%
  rename(value=lnlike)
c <- bind_rows(a, b)
ggplot(c) + 
  geom_violin(aes(x=as.numeric(timepoint),group=as.numeric(timepoint),
                  y=value,fill=ver),scale="width") +
  facet_wrap(~ver,nrow=2,scales="free")


time_vec <- seq(as.Date("2020-04-15")-35,as.Date("2020-12-01"),by="1 day")
chains_all$t <- time_vec[as.numeric(chains_all$timepoint)]
trajs1_all$date <- time_vec[as.numeric(trajs1_all$timepoint)]
vl_traj_comb$date <- time_vec[as.numeric(vl_traj_comb$timepoint)]

## Entire fit results
setwd("~/Documents/GitHub/virosolver_paper/mcmc_chains/4.real_ma_ct/real_ma_free_gp_excl_last_week/")

parTab<- read.csv("~/Documents/GitHub/virosolver_paper/pars/partab_gp_model.csv")
means <- parTab$values
names(means) <- parTab$names

chains <- load_mcmc_chains(getwd(),parTab,FALSE,1,30000,FALSE,FALSE,FALSE)

chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
samps <- sample(unique(chain$sampno),1000)
ages <- 1:100
vl_trajs <- matrix(0, nrow=length(samps),ncol=length(ages))
trajs <- matrix(0, nrow=length(samps),ncol=length(ages))
for(i in seq_along(samps)) {
  tmp_pars <- get_index_pars(chain, samps[i])
  tmp_pars["t_unit"] <- 1
  vl <- viral_load_func(tmp_pars,ages,TRUE)
  vl_trajs[i,]  <- rgumbel(length(vl),vl,tmp_pars["obs_sd"])
  trajs[i,] <- sapply(ages, function(a) prop_detectable(a, tmp_pars,vl))
}


trajs1 <- t(apply(trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
trajs1 <- as.data.frame(trajs1)
colnames(trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
trajs1$t <- ages
trajs1$timepoint <- as.character(min(as.numeric(timepoints))-7)
trajs1$ver <- "GP"

vl_trajs1 <- t(apply(vl_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
vl_trajs1 <- as.data.frame(vl_trajs1)
colnames(vl_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
vl_trajs1$t <- ages
vl_trajs1$timepoint <- as.character(min(as.numeric(timepoints))-7)
vl_trajs1$ver <- "GP"

chains <- load_mcmc_chains(getwd(),parTab,FALSE,1,30000,FALSE,FALSE,FALSE)
chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
chain <- chain[,colnames(chain) != "prob"]
chain_melted <- reshape2::melt(chain,id.vars="sampno")
chain_melted$timepoint <- as.character(min(as.numeric(timepoints))-7)
chain_melted$ver <- "GP"

chain_melted$t <- time_vec[as.numeric(chain_melted$timepoint)]
trajs1$date <- time_vec[as.numeric(trajs1$timepoint)]
vl_trajs1$date <- time_vec[as.numeric(vl_trajs1$timepoint)]


## Get prior draws
parTab <- read.csv("~/Documents/GitHub/virosolver_paper/pars/partab_fitted.csv")
pars <- parTab$values
names(pars) <- parTab$names
free_pars <- parTab[parTab$fixed == 0 & parTab$names != "prob", "names"]

sds <- c("beta"=0.25,"obs_sd"=1,"viral_peak"=1,
         "wane_rate2"=3,"t_switch"=3,"level_switch"=1,
         "prob_detect"=0.03)

ages <- 1:100

trajs_prior <- matrix(0, n_samp, length(ages))
vl_trajs_prior <- matrix(pars["true_0"], n_samp,length(ages))
sampd_pars <- matrix(0, n_samp, length(pars))

for(j in 1:n_samp){
  tmp_pars <- pars
  for(i in 1:length(pars)) {
    par_name <- names(pars)[i]
    if(par_name %in% free_pars) {
      if(par_name == "prob_detect"){
        beta1_mean <- pars["prob_detect"]
        beta1_sd <- sds["prob_detect"]
        beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
        beta_beta <- beta_alpha*(1/beta1_mean - 1)
        tmp_pars[i] <- -1
        while(tmp_pars[i] < 0){
          tmp_pars[i] <- rbeta(1,beta_alpha,beta_beta)
        }
        
      } else {
        tmp_pars[par_name] <- -1
        while(tmp_pars[par_name] < 0){
          tmp_pars[par_name] <- rnorm(1,pars[par_name],sds[par_name])
        }
      }
    }
  }
  names(tmp_pars) <- names(pars)
  sampd_pars[j,] <- tmp_pars
  vl <- viral_load_func(tmp_pars, ages, TRUE)
  vl1 <- viral_load_func(tmp_pars, ages, FALSE)
  vl_trajs_prior[j,] <- extraDistr::rgumbel(length(vl), vl,tmp_pars["obs_sd"])
  trajs_prior[j,] <- sapply(ages, function(a) prop_detectable(a, tmp_pars,vl))
}


trajs_prior <- t(apply(trajs_prior, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
trajs_prior <- as.data.frame(trajs_prior)
colnames(trajs_prior) <- c("lower","lower_mid","median","upper_mid","upper")
trajs_prior$t <- ages
trajs_prior$timepoint <- as.character(min(as.numeric(timepoints))-14)
trajs_prior$ver <- "Prior"

vl_trajs_prior <- t(apply(vl_trajs_prior, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
vl_trajs_prior <- as.data.frame(vl_trajs_prior)
colnames(vl_trajs_prior) <- c("lower","lower_mid","median","upper_mid","upper")
vl_trajs_prior$t <- ages
vl_trajs_prior$timepoint <- as.character(min(as.numeric(timepoints))-14)
vl_trajs_prior$ver <- "Prior"

colnames(sampd_pars) <- names(tmp_pars)
sampd_pars <- as.data.frame(sampd_pars)
sampd_pars$sampno <- 1:nrow(sampd_pars)
sampd_pars_melted <- reshape2::melt(sampd_pars,id.vars="sampno")
sampd_pars_melted$timepoint <- as.character(min(as.numeric(timepoints))-14)
sampd_pars_melted$ver <- "Prior"


sampd_pars_melted$t <- time_vec[as.numeric(sampd_pars_melted$timepoint)]
vl_trajs_prior$date <- time_vec[as.numeric(vl_trajs_prior$timepoint)]
trajs_prior$date <- time_vec[as.numeric(trajs_prior$timepoint)]



chains_all$ver <- "SEIR"
chains_all <- bind_rows(chains_all, chain_melted,sampd_pars_melted)

p_posteriors <- ggplot(chains_all %>% filter(variable %in% c("level_switch","obs_sd","prob_detect","t_switch",
                                                             "viral_peak","wane_rate2"))) +
  geom_violin(aes(x=t,y=value,group=t,fill=ver),draw_quantiles=c(0.025,0.5,0.975),scale="width") +
  geom_hline(data=parTab[parTab$names %in% c("level_switch","obs_sd","prob_detect","t_switch",
                                             "viral_peak","wane_rate2"),] %>% rename(variable=names),
             aes(yintercept=values))+
  facet_wrap(~variable,scales="free_y",ncol=1) +
  scale_x_date(breaks="2 month") +
  export_theme


vl_traj_comb$ver <- "SEIR"
vl_traj_comb <- bind_rows(vl_traj_comb, vl_trajs1,vl_trajs_prior)

p_vl <- ggplot(vl_traj_comb) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=ver),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid,fill=ver),alpha=0.5) +
  geom_line(aes(x=t,y=median,col=ver)) +
  xlab("Days since infection") +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  coord_cartesian(ylim=c(40,0)) +
  ylab("Ct value") +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme+
  facet_wrap(~date,ncol=5)

trajs1_all$ver <- "SEIR"
trajs1_all <- bind_rows(trajs1_all, trajs1,trajs_prior)

p_detect <- ggplot(trajs1_all) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=ver),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid,fill=ver),alpha=0.5) +
  geom_line(aes(x=t,y=median,col=ver)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme+
  facet_wrap(~date,ncol=5)

ggsave("~/Documents/GitHub/virosolver_paper/figures/posterior_checks/ma_posteriors.png",
       p_posteriors,height=10,width=8,units="in",dpi=300)
ggsave("~/Documents/GitHub/virosolver_paper/figures/posterior_checks/ma_prop_detect.png",
       p_detect,height=10,width=8,units="in",dpi=300)
ggsave("~/Documents/GitHub/virosolver_paper/figures/posterior_checks/ma_vl.png",
       p_vl,height=10,width=8,units="in",dpi=300)
