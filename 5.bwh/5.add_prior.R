## Generate prior samples
######################################################
test_ages <- seq(0,50,by=1)

## Generate draws from our priors and get quantiles
sample_pars_from_prior <- function(pars){
  sds <- sds_gp
  ## Draw parameters from normal priors
  tmp_pars <- pars
  tmp_pars["obs_sd"] <- max(rnorm(1, pars["obs_sd"], sds["obs_sd"]), 1)
  tmp_pars["viral_peak"] <- rnorm(1, pars["viral_peak"], sds["viral_peak"])
  tmp_pars["t_switch"] <- rnorm(1,pars["t_switch"],sds["t_switch"])
  
  #tmp_pars["wane_rate2"] <- rnorm(1,pars["wane_rate2"],sds["wane_rate2"])
  #tmp_pars["level_switch"] <- min(rnorm(1,pars["level_switch"],sds["level_switch"]), pars["intercept"])
  
  beta1_mean <- pars["prob_detect"]
  beta1_sd <- sds["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  tmp_pars["prob_detect"] <- rbeta(1,beta_alpha,beta_beta)
  tmp_pars
}

## Draw 10000 samples from prior
n_samp1 <- 10000
tmp_trajs <- matrix(0, nrow=n_samp1, ncol=length(test_ages))
tmp_trajs2 <- matrix(0, nrow=n_samp1, ncol=length(test_ages))

pars <- parTab$values
names(pars) <- parTab$names
all_pars <- matrix(nrow=n_samp1,ncol=length(pars))

## For 10000 samples
for(i in 1:n_samp1){
  ## Get parameters
  tmp_pars <- sample_pars_from_prior(pars)
  ## Get trajectory Ct scale
  tmp_vls <- viral_load_func(tmp_pars, test_ages, FALSE)
  ## Get trajectory viral load scale
  tmp_vls2 <- viral_load_func(tmp_pars, test_ages, TRUE)
  ## Get proportion detectable, but 0 on day 0
  probs <- c(0,prop_detectable(test_ages[test_ages > 0], tmp_pars, tmp_vls))
  if(any(is.na(probs))) print(tmp_pars)
  ## Save
  tmp_trajs[i,] <- probs
  tmp_trajs2[i,] <- tmp_vls
  all_pars[i,] <- tmp_pars
}

colnames(all_pars) <- names(pars)
all_pars <- all_pars[,c("viral_peak","obs_sd","t_switch","prob_detect")]
all_pars <- as.data.frame(all_pars)
all_pars$sampno <- 1:nrow(all_pars)
all_pars_melted <- all_pars %>% pivot_longer(-sampno)
all_pars_melted$chain <- "Prior"


## Get median, 95 and 50% quantiles
prior_prob_quants <- t(apply(tmp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
prior_prob_quants <- as.data.frame(prior_prob_quants)
colnames(prior_prob_quants) <- c("lower","mid_low","median","mid_high","upper")
prior_prob_quants$t <- test_ages


prior_vl_quants <- t(apply(tmp_trajs2, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
prior_vl_quants <- as.data.frame(prior_vl_quants)
colnames(prior_vl_quants) <- c("lower","mid_low","median","mid_high","upper")
prior_vl_quants$t <- test_ages

