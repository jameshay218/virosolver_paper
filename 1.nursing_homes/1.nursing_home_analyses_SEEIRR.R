##########################################
# 2. Fit SEEIRR model
##########################################
parTab <- read.csv("pars/generic/partab_seeirr_combined.csv",stringsAsFactors=FALSE)

test_institutions <- c("NH 1","NH 2", "NH 3","NH 4")

## Make min_date 0
nh_prev1 <- dat_wide %>%
  filter(result %in% c("POS","NEG")) %>%
  filter(location %in% test_institutions) %>%
  group_by(location, result, week, mean_week_date) %>%
  tally() %>%
  group_by(location, result, mean_week_date, week) %>%
  summarize(n=sum(n)) %>%
  pivot_wider(values_from=n, names_from=result,values_fill=list(n=0)) %>%
  group_by(location, mean_week_date) %>%
  mutate(prev=POS/(POS+NEG),
         lower_confint=prop.test(POS, POS+NEG)$conf.int[1],
         upper_confint=prop.test(POS,POS+NEG)$conf.int[2]) %>%
  ungroup()  %>%
  rename(date=mean_week_date) %>%
  mutate(date = as.numeric(date - min_date))

## Test that function works
posterior <- create_post_func_seeirr_combined(parTab, data=nh_prev1,
                                              ts=times,PRIOR_FUNC=prior_func_seeirr,ver="likelihood")
dat <- posterior(parTab$values)

## Optional flag to rerun MCMC or not
if(rerun_mcmc_seeirr){
  res <- foreach(i=1:nchains,.packages = c("lazymcmc","rethinking","extraDistr","ggthemes","tidyverse","deSolve")) %dopar% {
    posterior <- create_post_func_seeirr_combined(parTab, data=nh_prev1,
                                                  ts=times,PRIOR_FUNC=prior_func_seeirr,ver="likelihood")
    dat <- posterior(parTab$values)
    
    devtools::load_all("~/Documents/GitHub/virosolver")
    ## Generate starting parameters
    startTab <- generate_start_tab(parTab)
    start_pars <- startTab$values

    while(!is.finite(posterior(start_pars))){
      startTab <- generate_start_tab(par_tab=startTab)
      start_pars <- startTab$values
    }
  
    ## Run MCMC part 1
    output <- run_MCMC(parTab=startTab, data=nh_prev1, mcmcPars=mcmcPars1,
                       filename=paste0(chainwd,"/chain_",i),
                       CREATE_POSTERIOR_FUNC = create_post_func_seeirr_combined, mvrPars = NULL,
                       PRIOR_FUNC=prior_func_seeirr, ts=times)
    chain <- read.csv(output$file)
    best_pars <- get_best_pars(chain)
    chain <- chain[chain$sampno >= mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]
    covMat <- cov(chain)
    mvrPars <- list(covMat,0.5,w=0.8)
    
    ## Start from best location of previous chain
    parTab$values <- best_pars
    
    ## Run MCMC part 2
    output <- run_MCMC(parTab=parTab, data=nh_prev1, mcmcPars=mcmcPars2,
                       filename=paste0(chainwd,"/chain_",i),
                       CREATE_POSTERIOR_FUNC = create_post_func_seeirr_combined, mvrPars = mvrPars,
                       PRIOR_FUNC=prior_func_seeirr, ts=times)
    chain <- read.csv(output$file)
    chain <- chain[chain$sampno > mcmcPars2["adaptive_period"],]
    chain$sampno <- chain$sampno + max(chain$sampno)*(i-1)
    chain
  }
  chain <- do.call("bind_rows",res)
} else {
  ## If not rerunning, use stored chains
  chain <- load_mcmc_chains(chainwd, parTab,FALSE,1,mcmcPars2["adaptive_period"],multi=TRUE)$chain
  chain <- as.data.frame(chain)
  chain$sampno <- 1:nrow(chain)
}

pdf("figures/supplement/seeirr_chains.pdf")
plot(coda::as.mcmc(chain))
dev.off()

## Get credible intervals etc
quants_all <- NULL
random_traj_all <- NULL
t0_all <- NULL
best_traj_all <- NULL
n_samp <- 1000 ## Number of posterior samples
samps <- sample(unique(chain$sampno), n_samp)
times1 <- 1:(max(nh_prev1$date) + 5) ## Final time to solve model up until

## Function to solve prevalence
model_func <- create_post_func_seeirr_combined(parTab, data=nh_prev1,
                                               ts=times1,INCIDENCE_FUNC=detectable_SEEIRRModel,
                                               PRIOR_FUNC=prior_func_seeirr,ver="model")
## Function to solve incidence
model_func_inc <- create_post_func_seeirr_combined(parTab, data=nh_prev1,
                                                   ts=times1,INCIDENCE_FUNC=solveSEEIRRModel_rlsoda_wrapper,
                                                   PRIOR_FUNC=prior_func_seeirr,ver="model")

## Get MAP trajectories
best_pars <- get_best_pars(chain)
best_traj <- model_func(best_pars) %>%
  mutate(time=as.Date(t + min_date)) %>%
  rename(location = loc)
best_inc_traj <- model_func_inc(best_pars) %>%
  mutate(time=as.Date(t + min_date)) %>%
  rename(location = loc,
         inc=prev)


dat <- NULL
dat_inc <- NULL 
## Get posterior draw trajectories
for(j in seq_along(samps)){
  pars <- get_index_pars(chain, samps[j])
  dat[[j]] <- model_func(pars)
  dat[[j]]$samp <- j
  
  tmp_inc <- model_func_inc(pars)
  dat_inc[[j]] <- tmp_inc
  dat_inc[[j]]$samp <- j
  
}
## Prevalence - quantiles
dat <- do.call("bind_rows",dat)
quants <- dat %>% group_by(loc, t) %>%
  summarize(lower=quantile(prev, 0.025),
            median=median(prev),
            upper=quantile(prev, 0.975))
quants <- quants %>%
  mutate(time=as.Date(t + min_date)) %>%
  rename(location=loc)

## Incidence - quantiles
dat_inc <- do.call("bind_rows",dat_inc)
quants_inc <- dat_inc %>% group_by(loc, t) %>%
  summarize(lower=quantile(prev, 0.025),
            median=median(prev),
            upper=quantile(prev, 0.975))
quants_inc <- quants_inc %>%
  mutate(time=as.Date(t + min_date)) %>%
  rename(location=loc)

## Get random trajectories
## Prevalence
random_traj <- dat %>%
  filter(samp %in% sample(unique(dat$samp), 25)) %>%
  mutate(time = as.Date(t + min_date)) %>%
  rename(location=loc)
## Incidence
random_traj_inc <- dat_inc %>%
  filter(samp %in% sample(unique(dat$samp), 25)) %>%
  mutate(time = as.Date(t + min_date)) %>%
  rename(location=loc,
         inc=prev)

## Get point-range for t0 posterior
t0_1 <- as.Date(quantile(chain$t0_1,c(0.025,0.5,0.975)) + min_date)
t0_2 <- as.Date(quantile(chain$t0_2,c(0.025,0.5,0.975)) + min_date)
t0_3 <- as.Date(quantile(chain$t0_3,c(0.025,0.5,0.975)) + min_date)
t0_4 <- as.Date(quantile(chain$t0_4,c(0.025,0.5,0.975)) + min_date)
t0_pointrange <- bind_rows(t0_1,t0_2,t0_3,t0_4)
colnames(t0_pointrange) <- c("lower","median","upper")
t0_pointrange$location <- test_institutions

## Get growth rates
dat_inc <- dat_inc %>% group_by(loc, samp) %>% mutate(gr = log(lead(prev,1)/prev)) %>%
  mutate(gr = ifelse(!is.finite(gr),NA,gr))
dat_inc <- dat_inc %>% mutate(avg_gr=zoo::rollapply(gr,max(35),mean,align='right',fill=NA))

## Quantiles on growth rates
gr_quants <- dat_inc %>% group_by(loc, t) %>%
  summarize(lower_gr=quantile(gr, 0.025, na.rm=TRUE),
            median_gr = median(gr, na.rm=TRUE),
            upper_gr = quantile(gr, 0.975, na.rm=TRUE),
            
            lower_gr_avg=quantile(avg_gr, 0.025, na.rm=TRUE),
            median_gr_avg = median(avg_gr, na.rm=TRUE),
            upper_gr_avg = quantile(avg_gr, 0.975, na.rm=TRUE)) %>%
  filter(t > min(nh_prev1$date) -15) %>%
  mutate(time = min_date + t)

