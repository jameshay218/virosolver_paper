## Standard deviations of viral kinetics priors
sds_exp <- sds_seir <- sds_gp <- c("beta"=0.25,"R0"=0.6,
                                   "obs_sd"=0.5,"viral_peak"=2,
                                   "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
                                   "prob_detect"=0.03,
                                   "incubation"=0.25, "infectious"=0.5,
                                   "rho"=2,"nu"=0.5)


## Prior for exponential growth model
prior_func_hinge_exp <- function(pars,...){
  prior1 <- dnorm(pars["beta"], 0, sds_exp["beta"],log=TRUE)
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_exp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_exp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_exp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_exp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_exp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_exp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  prior1 + obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior + level_prior + beta_prior
}
## Prior for SEIR model
prior_func_hinge_seir <- function(pars,...){
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
  #r0_prior <- dlnorm(pars["R0"],log(2),sds_seir["R0"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_seir["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
  infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
  
  obs_sd_prior + viral_peak_prior +
    wane_2_prior + tswitch_prior + level_prior + beta_prior +
    incu_prior + infectious_prior# + r0_prior
}

## Prior for GP version
prior_func_hinge_gp <- function(pars, t_dist){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ## Gaussian process prior, un-centered version
  k <- pars[which(par_names=="prob")]
  ## Leave this - correct for uncentered version as per Chapter 14 Statistical Rethinking
  prob_priors <- sum(dnorm(k, 0, 1, log=TRUE))
  
  nu_prior <- dexp(pars["nu"], 1/means[which(names(means) == "nu")],log=TRUE)
  rho_prior <- dexp(pars["rho"], 1/means[which(names(means) == "rho")],log=TRUE)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + prob_priors +
    nu_prior + rho_prior
}

## Prior for SEEIRR transmission model
## Note that the R0s and t0s are independent, the infection course parameters are shared
prior_func_seeirr <- function(pars){
  sds <- c("latent"=0.3,"incubation"=0.3,"infectious"=0.3,"recovery"=0.3)
  names(pars) <- parTab$names
  ## Uniform priors in parTab
  p1 <- dlnorm(pars["R0_1"], log(2), 0.6,log=TRUE)
  p1b <- dlnorm(pars["R0_2"], log(2), 0.6,log=TRUE)
  p1c <- dlnorm(pars["R0_3"], log(2), 0.6,log=TRUE)
  p1d <- dlnorm(pars["R0_4"], log(2), 0.6,log=TRUE)
  p2 <- dlnorm(pars["latent"],log(2),sds["latent"],log=TRUE)
  p3 <- dlnorm(pars["incubation"],log(2),sds["incubation"],log=TRUE)
  p4 <- dlnorm(pars["infectious"],log(4),sds["infectious"],log=TRUE)
  p5 <- dlnorm(pars["recovery"],log(11), sds["recovery"], log=TRUE)
  return(sum(p1,p1b,p1c,p1d, p2,p3,p4,p5))
}
## Prior for SEEIRR transmission model
## Note that the R0s and t0s are independent, the infection course parameters are shared
prior_func_seeirr_single <- function(pars){
  sds <- c("latent"=0.25,"incubation"=0.25,"infectious"=0.5,"recovery"=0.25)
  names(pars) <- parTab$names
  ## Uniform priors in parTab
  #p1 <- dlnorm(pars["R0"], log(2), 0.6,log=TRUE)
  p2 <- dlnorm(pars["latent"],log(2),sds["latent"],log=TRUE)
  p3 <- dlnorm(pars["incubation"],log(2),sds["incubation"],log=TRUE)
  p4 <- dlnorm(pars["infectious"],log(4),sds["infectious"],log=TRUE)
  p5 <- dlnorm(pars["recovery"],log(11), sds["recovery"], log=TRUE)
  return(sum(p2,p3,p4,p5))
}

## Prior for SEIR model
prior_func_seir_single <- function(pars,...){
  incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
  infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
  
  incu_prior + infectious_prior# + r0_prior
}
