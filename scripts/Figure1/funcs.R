
#### Functions for Viral Load Distributions: ####
## Describes average viral load as a function of time since infection
viral_load_func <- function(pars, a, convert_ct=TRUE){
  viral_load <- viral_load_func_asymp(pars, obs_t=a, convert_ct)
  viral_load
}

## Function to give probability of observing x given age a and the viral kinetics curve
p_a_use <- function(x, a, pars, obs_model="normal", log_prob=FALSE, ct=TRUE, additional_detect) {
  if (!ct) {
    x <- pars["intercept"]-log2(10)*(x-pars["LOD"])
  }
  p <- p_a(x, a, pars, viral_load_func(pars, a, convert_ct=TRUE), obs_model, log_prob, additional_detect)
  p
}

## Proportion of individuals with detectable viral loads for a given age
prop_det_use <- function(a, pars, obs_model="normal", additional_detect) {
  phi <- prop_detectable(a, pars, viral_load_func(pars, a, convert_ct=TRUE), obs_model, additional_detect)
  phi
}

## Function to generate random viral loads given a single age a and the viral kinetics curve
rx_a <- function(n, a, pars, obs_model="normal", additional_detect=FALSE) {
  ct_mean <- viral_load_func(pars, a, convert_ct=TRUE)
  if(obs_model=="normal") {
    xs <- rnorm(n, mean=ct_mean, sd=pars["obs_sd"])
  } else {
    xs <- rgumbel(n, mu=ct_mean, sigma=pars["obs_sd"])
  }
  xs <- xs[xs < pars["intercept"]]
  if(additional_detect){
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    if (a > t_switch) {
      ts <- a - t_switch
      additional_prob <- (1-pars["prob_detect"]*pars["t_unit"])^ts
      xs <- xs[rbinom(n, size=1, prob=additional_prob)==1]
    }
  }
  xs
}

## Function to generate random viral load trajectories:
rx <- function(pars, obs_model="normal", additional_detect=FALSE, a.range=1:35, corr=.8) {
  ct_mean <- viral_load_func(pars, a.range, convert_ct=TRUE)
  if (obs_model=="normal") {
    Sigma <- matrix(pars["obs_sd"]^2, nrow=length(a.range), ncol=length(a.range))
    for (i in 1:length(a.range)) {
      for (j in 1:length(a.range)) {
        Sigma[i,j] <- Sigma[i,j]*corr^(abs(i-j))
      }
    }
    xs <- rmvnorm(1, mean=ct_mean, sigma=Sigma)[1,]
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  } else {
    ind <- extraDistr::rgumbel(1, mu=0, sigma=pars["obs_sd"])
    xs <- ct_mean + ind
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  }
  if (additional_detect) {
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    t_switch_int <- ceiling(t_switch)
    t_switch_inc <- t_switch_int-t_switch
    ts <- a.range[a.range >= t_switch_int]
    addl_det <- c(rbinom(1, 1, (1-pars["prob_detect"]*pars["t_unit"])^t_switch_inc),
                  rbinom(length(ts)-1, 1, (1-pars["prob_detect"]*pars["t_unit"])))
    last.day <- min(c(max(a.range)+1,ts[addl_det==0]))-1
    xs[a.range > last.day] <- NA
  }
  xs
}