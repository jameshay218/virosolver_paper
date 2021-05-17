##########################################
# 3. Fit Ct model
##########################################
## For each analysis being run
betas_all_models <- NULL ## Growth rates for all analyses
chains_all_models <- NULL ## MCMC chains
prop_detect_models <- NULL ## Proportion detectable over time
vl_all_models <- NULL ## Predicted viral kinetics
prior_draws_all <- NULL ## Prior draws
preds_all_models <- NULL ## Predicted trajectories (prob infection)
best_prob_models <- NULL ## MAP predicted trajectories
obs_dat_models <- NULL ## Data used, after cleaning
preds_fits_all_comb_models <- NULL ## Predicted Ct distributions

###
for(run in 1:nrow(analysis_control)){
  print(analysis_control[run,])
  ## Extract analysis control parameters
  run_version <- analysis_control$version[run]
  use_pos <- analysis_control$use_pos[run]
  cumu_data <- analysis_control$cumu_data[run]
  par_tab_file <- analysis_control$par_table_file[run]
  prior_func_use <- analysis_control$prior_func[[run]]
  inc_func_use <- analysis_control$inc_func[[run]]
  max_age <- analysis_control$age_max[run]
  
  ## Read in parameter table
  parTab <- read.csv(par_tab_file,stringsAsFactors = FALSE)
  pars <- parTab$values
  names(pars) <- parTab$names
  
  ## Means for priors
  means <- parTab$values
  names(means) <- parTab$names
  
  dir_name_tmp <- paste0(run_version,"_",use_pos,"_",cumu_data)
  
  ## If flag set to rerun, set up data and time periods, then fit the model
  if(rerun_mcmc_ct){
    ## For each location and each week, fit 3 chains
    res <- foreach(i=1:nrow(unique_data_combs),.packages = c("rethinking","extraDistr","ggthemes","tidyverse","deSolve","patchwork")) %dopar% {
      runname <- paste0("ct_", unique_data_combs$location[i], "_",unique_data_combs$mean_week_date[i])
      dir.create(paste0(chainwd2,"/",dir_name_tmp,"/",runname),recursive = TRUE)
      ## The level above this should be the folder with all of your Git repos
      ## NOTE that foreach can sometimes be funny with working directories with dopar
      devtools::load_all(paste0(GIT_WD,"/lazymcmc")) ## Parallel tempering branch
      mcmc_pars_use <- mcmcPars_ct_pt
      devtools::load_all(paste0(GIT_WD,"/virosolver"))
      
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
      
      ## Test that model can be run without failure
      f <- create_posterior_func(parTab, dat_tmp_used_in_run, CREATE_POSTERIOR_FUNC=create_posterior_func,
                                 INCIDENCE_FUNC=inc_func_use,
                                 PRIOR_FUNC=prior_func_use,
                                 use_pos=use_pos,
                                 t_dist=NULL,solve_ver = "model")
      f(means)
      
      ## Run each chain sequentially here
      chains <- NULL
      for(j in 1:nchains){
        ## Get random starting values
        startTab <- rep(list(parTab),n_temperatures)
        for(k in 1:length(startTab)){
          startTab[[k]] <- generate_viable_start_pars(parTab=parTab,
                                                      obs_dat=dat_tmp_used_in_run,
                                                      CREATE_POSTERIOR_FUNC=create_posterior_func,
                                                      INCIDENCE_FUNC=inc_func_use,
                                                      PRIOR_FUNC=prior_func_use,
                                                      use_pos=use_pos,
                                                      t_dist=NULL)
        }
        covMat <- diag(nrow(startTab[[1]]))
        mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
        mvrPars <- rep(list(mvrPars), n_temperatures)
        
        chains <- NULL
        output <- run_MCMC(parTab=startTab,
                           data=dat_tmp_used_in_run,
                           INCIDENCE_FUNC=inc_func_use,
                           PRIOR_FUNC = prior_func_use,
                           solve_likelihood=TRUE,
                           mcmcPars=mcmc_pars_use,
                           filename=paste0(chainwd2,"/", dir_name_tmp,"/",runname,"/",runname,"_chainno_",j),
                           CREATE_POSTERIOR_FUNC=create_posterior_func,
                           mvrPars=mvrPars,
                           OPT_TUNING=0.2,
                           use_pos=use_pos,
                           t_dist=t_dist)
        chain <- read.csv(output$file)
        chain <- chain[chain$sampno > mcmc_pars_use["adaptive_period"],]
        chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
        chains[[j]] <- chain
      }
      chain <- do.call("bind_rows",chains)
    }
  } 
  
  ## For each model fit, go through and read in the MCMC chains
  res <- NULL
  diagnostics_all <- NULL
  prob_infection_preds <- NULL
  prob_infection_preds_map <- NULL
  pred_fits_all <- NULL
  chains_all <- NULL
  obs_dat_all <- NULL
  model_func_all <- NULL
  
  for(i in 1:nrow(unique_data_combs)){
    PTchain <- TRUE
    mcmc_pars_use <- mcmcPars_ct_pt
    runname <- paste0("ct_", unique_data_combs$location[i], "_",unique_data_combs$mean_week_date[i])
    chainwd_tmp <- paste0(chainwd2,"/", dir_name_tmp,"/",runname)
    
    ## Convergence diagnostics
    chains_diag <- load_mcmc_chains(chainwd_tmp, parTab,TRUE,1,
                                    mcmc_pars_use["adaptive_period"],multi=TRUE,
                                    chainNo=FALSE,PTchain=PTchain)
    if(length(chains_diag$list) > 1){
      gelman_tmp <- lazymcmc::gelman_diagnostics(chains_diag$list)
      ess_tmp <- coda::effectiveSize(chains_diag$list)
      diagnostics_all[[i]] <- tibble(psrf_point=gelman_tmp$GelmanDiag$psrf[,1], psrf_upper=gelman_tmp$GelmanDiag$psrf[,2],
                                     ess=ess_tmp, runname=paste0(dir_name_tmp,"/",runname),parname=names(ess_tmp))
    }
    
    ## Read in chain for plotting and summaries
    chain <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmc_pars_use["adaptive_period"],
                              multi=TRUE,chainNo=TRUE,PTchain=PTchain)$chain
    chain <- as.data.frame(chain)
    
    chain_save <- chain
    chain_save$loc <- unique_data_combs$location[i]
    chain_save$t <- unique_data_combs$mean_week_date[i]
    chain_save$ver <- run_version
    chain_save$use_pos <- use_pos
    res[[i]] <- chain
  }
  ## Combine convergence diagnostics and save combined results
  if(nchains > 1){
    diagnostics_all <- do.call("bind_rows", diagnostics_all)
    dir.create(paste0(plot_wd, "/",dir_name_tmp),recursive = TRUE)
    write_csv(diagnostics_all, paste0(paste0(plot_wd, "/",dir_name_tmp), "/mcmc_diagnostics.csv"))
  }
  
  ## For each time point and nursing home fit, generate trajectories, predictions etc
  betas <- NULL
  chains_all_comb <- NULL
  vl_trajs_all <- NULL
  prop_detect_traj_all <- NULL
  for(i in 1:nrow(unique_data_combs)) {
    print(i)
    runname <- paste0("ct_", unique_data_combs$location[i], "_",unique_data_combs$mean_week_date[i])
    plot_wd_tmp <- paste0(plot_wd, "/",dir_name_tmp)
    
    #######
    ## Set up data as above. This is copied and pasted from the section in the rerun_mcmc part
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
    
    ## Store data used here for later plotting
    obs_dat_all[[i]] <- dat_tmp_used_in_run
    
    ## This is for the GP version
    mat <- matrix(rep(times, each=length(times)),ncol=length(times))
    t_dist <- abs(apply(mat, 2, function(x) x-times))
    if(run_version == "gp"){
      parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
      pars <- parTab$values
      names(pars) <- parTab$names
      
      ## Means for priors
      means <- parTab$values
      names(means) <- parTab$names
    }
    chain <- res[[i]]
    chain_comb <- chain
    chain_comb$sampno <- 1:nrow(chain_comb)
    
    chain1 <- chain
    if("prob" %in% colnames(chain)){
      use_cols <- c(which(colnames(chain) != "prob"), which(colnames(chain) == "prob")[1])
      chain1 <- chain[,use_cols]
    }
    chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
    
    ## Plot MCMC traces and densities
    p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
      mutate(chain = as.factor(chain)) %>%
      pivot_longer(-c(sampno,chain)) %>%
      ggplot() +
      geom_line(aes(x=sampno,y=value,col=chain)) +
      facet_wrap(~name,scales="free_y")+
      scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),by=2000)) +
      export_theme
    
    p_densities <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
      mutate(chain = as.factor(chain)) %>%
      pivot_longer(-c(sampno,chain)) %>%
      ggplot() +
      geom_density(aes(x=value,fill=chain),alpha=0.25) +
      facet_wrap(~name,scales="free") +
      export_theme
    
    ## Get model predicted Ct distribution 
    predictions <- plot_prob_infection(chain_comb, n_samp, inc_func_use,times,obs_dat=dat_tmp_used_in_run)
    
    ## Save projected trajectories
    predictions_traj <- predictions$predictions
    predictions_traj$loc <- unique_data_combs$location[i]
    predictions_traj$samp_t <- unique_data_combs$mean_week_date[i]
    predictions_traj$ver <- run_version
    predictions_traj$use_pos <- use_pos 
    
    ## Save projected MAP trajectories
    predictions_traj_map <- predictions$map_prediction
    predictions_traj_map$loc <- unique_data_combs$location[i]
    predictions_traj_map$samp_t <- unique_data_combs$mean_week_date[i]
    predictions_traj_map$ver <- run_version
    predictions_traj_map$use_pos <- use_pos 
    
    prob_infection_preds[[i]] <- predictions_traj
    prob_infection_preds_map[[i]] <- predictions_traj_map
    
    ## Get model predicted Ct distributions and plots
    p1 <- predictions$plot
    model_func <- create_posterior_func(parTab,dat_tmp_used_in_run,NULL,inc_func_use,"model")
    model_func_all[[i]] <- model_func
    
    predicted_fits_dat <- predicted_distribution_fits(chain_comb, dat_tmp_used_in_run, model_func,n_samp)
    predicted_fits_dat$loc <- unique_data_combs$location[i]
    predicted_fits_dat$samp_t <- unique_data_combs$mean_week_date[i]
    predicted_fits_dat$ver <- run_version
    predicted_fits_dat$use_pos <- use_pos 
    pred_fits_all[[i]] <- predicted_fits_dat
    
    p2 <- plot_distribution_fits(chain_comb, dat_tmp_used_in_run, model_func,n_samp)
    
    samps <- sample(unique(chain_comb$sampno),n_samp)
    trajs <- matrix(0, nrow=n_samp,ncol=length(times))
    vl_trajs <- matrix(0, nrow=n_samp, ncol=100)
    prop_detect_tmp <- matrix(0, nrow=n_samp, ncol=100)
    
    ## Get incidence curve and viral load samples
    for(ii in seq_along(samps)){
      tmp_pars <- get_index_pars(chain_comb, samps[ii])
      trajs[ii,] <- pmax(inc_func_use(tmp_pars,times),0.0000001)
      vl <- viral_load_func(tmp_pars,1:100,FALSE)
      vl_noise <- extraDistr::rgumbel(length(vl),vl,tmp_pars["obs_sd"])
      vl_trajs[ii,] <- vl_noise
      
      vl1 <- viral_load_func(tmp_pars,1:100,FALSE)
      prop_detect_tmp[ii,] <- sapply(1:100, 
                                     function(a) prop_detectable_cpp(a,vl[a],tmp_pars["obs_sd"],
                                                                     tmp_pars["intercept"],tmp_pars["t_switch"],
                                                                     tmp_pars["prob_detect"]))
      
    }
    
    ## Store viral kinetics trajectory
    vl_trajs_melted <- reshape2::melt(vl_trajs)
    colnames(vl_trajs_melted) <- c("sampno","age","ct")
    vl_trajs_melted$t <-unique_data_combs$mean_week_date[i]
    vl_trajs_melted$loc <- unique_data_combs$location[i]
    vl_trajs_melted$ver <- dir_name_tmp
    vl_trajs_all[[i]] <- vl_trajs_melted
    
    ## Proportion detectable
    prop_detect_melted <- reshape2::melt(prop_detect_tmp)
    colnames(prop_detect_melted) <- c("sampno","age","p")
    prop_detect_melted$t <-unique_data_combs$mean_week_date[i]
    prop_detect_melted$loc <- unique_data_combs$location[i]
    prop_detect_melted$ver <- dir_name_tmp
    prop_detect_traj_all[[i]] <- prop_detect_melted
    
    trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
    trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
    trajs1_quants <- as.data.frame(trajs1_quants)
    trajs1_quants$t <- 1:nrow(trajs1_quants)
    colnames(trajs1_quants) <- c("lower","median","upper","t")
    p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
      geom_line(aes(x=t,y=median)) + 
      coord_cartesian(ylim=c(-0.5,0.5))
    
    if(run_version == "exp"){
      betas[[i]] <- tibble(beta=chain_comb$beta,loc=unique_data_combs$location[i],
                           t=unique_data_combs$mean_week_date[i])
    } else {
      betas[[i]] <- tibble(beta=trajs1[,ncol(trajs1)],loc=unique_data_combs$location[i],
                           t=unique_data_combs$mean_week_date[i])
    }
    chain_melted <- reshape2::melt(chain_comb[,c("sampno","viral_peak","obs_sd","t_switch","level_switch",
                                                 "wane_rate2","prob_detect")], id.vars="sampno")
    chain_melted$t <- unique_data_combs$mean_week_date[i]
    chain_melted$loc <- unique_data_combs$location[i]
    chain_melted$ver <- dir_name_tmp
    chains_all_comb[[i]] <- chain_melted
    
    ## Save plots
    dir.create(paste0(plot_wd_tmp,"/traces/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/predictions/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/distributions/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/posteriors/"),recursive = TRUE)
    dir.create(paste0(plot_wd_tmp,"/grs/"),recursive = TRUE)
    
    ggsave(paste0(plot_wd_tmp,"/traces/",runname,"_trace.png"),p_trace,width=7,height=4)
    ggsave(paste0(plot_wd_tmp,"/predictions/",runname,"_predictions.png"),p1,width=7,height=4)
    ggsave(paste0(plot_wd_tmp,"/distributions/",runname,"_distributions.png"),p2,
           width=(7/5) * length(unique(dat_tmp$t)),height=6)
    ggsave(paste0(plot_wd_tmp,"/posteriors/",runname,"_densities.png"),p_densities,width=7,height=4)
    ggsave(paste0(plot_wd_tmp,"/grs/",runname,"_grs.png"),p_gr,width=7,height=4)
  }
  betas_all <- do.call("bind_rows",betas)
  betas_all$ver <- dir_name_tmp
  betas_all_models[[run]] <- betas_all
  
  ## Combine viral load trajectory estimates
  vl_all <- do.call("bind_rows",vl_trajs_all)
  vl_all$ver <- dir_name_tmp
  vl_all_models[[run]] <- vl_all
  
  ## Combine proportion detectable estimates
  prop_detect_all <- do.call("bind_rows",prop_detect_traj_all)
  prop_detect_all$ver <- dir_name_tmp
  prop_detect_models[[run]] <- prop_detect_all
  
  ## Combine MCMC chains
  chains_all <- do.call("bind_rows",chains_all_comb)
  chains_all$ver <- dir_name_tmp
  chains_all_models[[run]] <- chains_all
  
  ## Generate prior draws for later
  if(run_version == "exp"){
    sds_use <- sds_exp
  } else {
    sds_use <- sds_seir
  }
  means <- parTab$values
  names(means) <- parTab$names
  means["level_switch"] <- 33
  
  prior_draws <- NULL
  for(par in parTab[parTab$fixed == 0,"names"]){
    if(par %in% names(sds_use)) {
      if(par == "prob_detect"){
        beta_alpha <- ((1-means[par])/sds_use[par]^2 - 1/means[par])*means[par]^2
        beta_beta <- beta_alpha*(1/means[par] - 1)
        tmp <- rbeta(10000,beta_alpha,beta_beta)
      } else {
        tmp <- rnorm(10000,mean=means[par],sd=sds_use[par])
      }
    } else {
      tmp <- runif(10000,parTab[parTab$names == par, "lower_bound"], parTab[parTab$names == par, "upper_bound"])
    }
    prior_draws[[par]] <- tibble(value=tmp,ver=dir_name_tmp,variable=par)
  }
  prior_draws_all[[run]] <- do.call("bind_rows",prior_draws)
  
  ## Plot supplementary figure of model fits
  preds_all <- do.call("bind_rows", prob_infection_preds)
  best_prob_dat <- do.call("bind_rows", prob_infection_preds_map)
  preds_all$t <- as.Date(preds_all$t, origin=min_date)
  best_prob_dat$t <- as.Date(best_prob_dat$t, origin=min_date)
  
  preds_all_models[[run]] <- preds_all
  best_prob_models[[run]] <- best_prob_dat
  obs_dat_models[[run]] <- do.call("bind_rows",obs_dat_all) %>% rename(samp_t = t) %>% mutate(samp_t = as.Date(samp_t, origin=min_date))
  preds_fits_all_comb_models[[run]] <- do.call("bind_rows",pred_fits_all)
  
  #if(run_version == "seir"){
  #  source("scripts/nursing_home_fit_supp.R")
  #}
  
}

## Growth rate estimates
betas_all_models_comb <- do.call("bind_rows",betas_all_models)
betas_all_quants <- betas_all_models_comb %>% group_by(loc, t, ver) %>%
  summarize(lower=quantile(beta,0.025),
            median=median(beta),
            upper=quantile(beta,0.975))

## Viral kinetics estimates
vl_trajs_all_comb <- do.call("bind_rows",vl_all_models)
vl_trajs_all_quants <- vl_trajs_all_comb %>% group_by(loc, age,t, ver) %>%
  summarize(lower=quantile(ct,0.025),
            lower_mid=quantile(ct,0.25),
            median=median(ct),
            upper_mid=quantile(ct,0.75),
            upper=quantile(ct,0.975))

## Prop detectable estimates
prob_detect_all_comb <- do.call("bind_rows",prop_detect_models)
prob_detect_all_quants <- prob_detect_all_comb %>% group_by(loc, age,t, ver) %>%
  summarize(lower=quantile(p,0.025),
            lower_mid=quantile(p,0.25),
            median=median(p),
            upper_mid=quantile(p,0.75),
            upper=quantile(p,0.975))

## Viral kinetics plots
p_vls <- ggplot(vl_trajs_all_quants %>% filter(ver=="seir_FALSE_FALSE")) +
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=age,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=age,y=median)) +
  xlab("Days since infection") +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  coord_cartesian(ylim=c(40,0)) +
  ylab("Viral load (log10 RNA copies/ml)") +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme +
  facet_grid(loc~t)

## Proportion detectable
p_detect <- ggplot(prob_detect_all_quants%>% filter(ver=="seir_FALSE_FALSE")) + 
  geom_ribbon(aes(x=age,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=age,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=age,y=median)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme +
  facet_grid(loc~t)



## MCMC chains
mcmc_all <- do.call("bind_rows", chains_all_models)
priors_all_comb <- do.call("bind_rows",prior_draws_all)

mcmc_all_thinned <- mcmc_all %>% group_by(loc, t, ver, variable) %>% sample_n(1000) %>%
  filter(variable %in% c("viral_peak","obs_sd","t_switch","level_switch","prob_detect"))

name_key <- c("viral_peak"="Ct[peak]","obs_sd"="sigma","t_switch"="t[switch]","level_switch"="Ct[switch]","prob_detect"="p[addl]",
              "beta"="beta","R0"="R[0]","t0"="t[0]")
loc_key <- c("NH 1"="NH1", "NH 2"="NH2", "NH 3"="NH3", "NH 4"="NH4")
mcmc_all_thinned$variable <- name_key[as.character(mcmc_all_thinned$variable)]
mcmc_all_thinned$loc <- loc_key[mcmc_all_thinned$loc]
parTab1 <- parTab %>% filter(names %in% c("viral_peak","obs_sd","t_switch","level_switch","prob_detect"))
parTab1$names <- name_key[parTab1$names]
mcmc_all_thinned <- mcmc_all_thinned %>% rename(Location=loc)

priors_all_comb1 <- priors_all_comb %>% filter(variable %in% c("viral_peak","obs_sd","t_switch","level_switch","prob_detect"))
priors_all_comb1$variable <- name_key[as.character(priors_all_comb1$variable)]
prior_ranges <- priors_all_comb1 %>% filter(ver=="seir_FALSE_FALSE") %>% group_by(ver, variable) %>%
  mutate(lower_bound=min(value),upper_bound=max(value))

p_par_priors <- ggplot(priors_all_comb1 %>% filter(ver=="seir_FALSE_FALSE")) + 
  geom_violin(aes(y=value,x=1),scale="width",draw_quantiles=c(0.025,0.5,0.975),fill="grey70") +
  geom_hline(data=parTab1 %>% rename(variable=names),
             aes(yintercept=values),linetype="dashed",col=AAAS_palette["grey1"]) +
  #geom_blank(data=parTab1 %>% rename(variable=names),aes(ymin=lower_bound,ymax=upper_bound)) +
  xlab("") +
  ylab("Prior estimate") +
  ggsci::scale_fill_aaas() +
  export_theme +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        strip.text=element_text(size=8,face="bold")) +
  facet_wrap(~variable,scales="free_y",ncol=1,labeller = label_parsed,strip.position="right") +
  ggtitle("Assumed priors") +
  labs(tag="C")

p_par_estimates_seir <- ggplot(mcmc_all_thinned %>% filter(ver=="seir_FALSE_FALSE")) + 
  geom_violin(aes(x=as.factor(t),y=value,group=as.factor(t),fill=Location),draw_quantiles=c(0.025,0.5,0.975)) +
  geom_hline(data=parTab1 %>% rename(variable=names),
             aes(yintercept=values),linetype="dashed",col=AAAS_palette["grey1"]) +
  geom_blank(data=prior_ranges,aes(ymin=lower_bound,ymax=upper_bound)) +
  xlab("Date") +
  ylab("Posterior estimate") +
  ggsci::scale_fill_aaas() +
  export_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),legend.position="bottom",
        strip.text=element_text(size=8,face="bold")) +
  facet_grid(variable~Location,scales="free_y",labeller = label_parsed) +
  ggtitle("SEIR fits") +
  labs(tag="A")

p_par_estimates_exp <- ggplot(mcmc_all_thinned %>% filter(ver=="exp_FALSE_FALSE")) + 
  geom_violin(aes(x=as.factor(t),y=value,group=as.factor(t),fill=Location),draw_quantiles=c(0.025,0.5,0.975)) +
  geom_hline(data=parTab1 %>% rename(variable=names),
             aes(yintercept=values),linetype="dashed",col=AAAS_palette["grey1"]) +
  geom_blank(data=prior_ranges,aes(ymin=lower_bound,ymax=upper_bound)) +
  xlab("Date") +
  ylab("Value") +
  ggsci::scale_fill_aaas() +
  export_theme +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),legend.position="bottom",
        strip.text=element_text(size=8,face="bold")) +
  facet_grid(variable~Location,scales="free_y",labeller = label_parsed)+
  ggtitle("Exponential model fits")+
  ylab("") +
  labs(tag="B")

lhs_p <- (p_par_estimates_seir| p_par_estimates_exp)
fig_S_posteriors <- (lhs_p | p_par_priors) + plot_layout(widths=c(3,3,1))

ggsave("figures/supplement/nh_posteriors.pdf",height=8,width=9)
ggsave("figures/supplement/nh_posteriors.png",height=8,width=9,units="in",dpi=300)