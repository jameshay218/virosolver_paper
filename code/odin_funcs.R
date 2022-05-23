## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  beta <- user(0.2)
  sigma <- user(0.15)
  gamma <- user(0.1)
  
}, verbose = FALSE)



## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_switch <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  beta <- beta1*(step <= t_switch1) + beta2*(step > t_switch1 && step < t_switch2) + beta3*(step >= t_switch2)
  
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  beta1 <- user(0.2)
  beta2 <- user(0.2)
  beta3 <- user(0.2)
  sigma <- user(0.15)
  gamma <- user(0.1)
  t_switch1 <- user(25)
  t_switch2 <- user(50)
  
  
}, verbose = FALSE)


## Adapted from http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_interpolate <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(inc) <- n_SE
  
  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * I / N) # S to E
  p_EI <- 1 - exp(-sigma) # E to I
  p_IR <- 1 - exp(-gamma) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)
  
  ## Total population size
  N <- S + E + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- 0
  initial(inc) <- 0
  
  ## User defined parameters - default in parentheses:
  beta = interpolate(betat,betay,"constant")
  output(beta) = beta
  betat[] = user()
  betay[] = user()
  dim(betat) <- user()
  dim(betay) <- length(betat)
  
  S_ini <- user(1000)
  E_ini <- user(0)
  I_ini <- user(1)
  sigma <- user(0.15)
  gamma <- user(0.1)
  
  
}, verbose = FALSE)


## Wrapper for SEIR model simulation, deterministic or stochastic
## INPUTS: 
##      1. population_n: size of population to simulate
##      2. solve_times: vector of times to solve model over
##      3. pars: named vector of SEIR model parameters
##      4. switch_model: if TRUE, uses the SEIR model with 2 switch points in transmission intensity
##      5. beta_smooth: spar parameter fed to smooth.spline, used to smooth the time-varying beta
## OUTPUTS: 
##      1. Data frame of the SEIR solution
##      2. Vector of absolute incidence per time point
##      3. Overall probability of infection
##      4. Plot of incidence over time
##      5. Average growth rates over different time periods
##      6. Plot of average growth rates over different time periods
simulate_seir_wrapper <- function(population_n, solve_times, pars, switch_model=FALSE,beta_smooth=0.8){
    ####################################################
    ## Stochastic model
    ####################################################
    if(switch_model){
      gamma1 <- 1/pars["infectious"]
      sigma1 <- 1/pars["incubation"]
      beta1 <- pars["R0_1"]*gamma1
      beta2 <- pars["R0_2"]*gamma1
      beta3 <- pars["R0_3"]*gamma1
      I0 <- ceiling(pars["I0"]*population_n)
      ## Odin stochastic SEIR model generator
      #seir <- seir_generator_switch(beta1=beta1,beta2=beta2,beta3=beta3,
      #                              sigma=sigma1,gamma=gamma1,
      #                              S_ini=population_n-I0,I_ini=I0,
      #                              t_switch1=pars["t_switch1"],t_switch2=pars["t_switch2"])
      betas <- rep(beta3, length(solve_times))
      betas[which(solve_times < pars["t_switch2"])] <- beta2
      betas[which(solve_times < pars["t_switch1"])] <- beta1
      betas <- smooth.spline(betas,spar=beta_smooth)$y
      seir <- seir_generator_interpolate$new(betat=solve_times,betay=betas,sigma=sigma1,gamma=gamma1,S_ini=population_n-I0,I_ini=I0)
      
    } else {
      gamma1 <- 1/pars["infectious"]
      sigma1 <- 1/pars["incubation"]
      beta1 <- pars["R0"]*gamma1
      I0 <- ceiling(pars["I0"]*population_n)
      ## Odin stochastic SEIR model generator
      seir <- seir_generator$new(beta=beta1,sigma=sigma1,gamma=gamma1,S_ini=population_n-I0,I_ini=I0)
    }
    ## Solve model
    res <- seir$run(solve_times)
    ## Make sure we get a simulation with an outbreak - keep trying until it takes off
    while(max(res[,"I"]) <= I0) res <- seir$run(solve_times)
    res <- as.data.frame(res)
    
    ## Shift for start
    res$step <- res$step + floor(pars["t0"])
    
    ## Dummy rows from pre-seeding
    if(pars["t0"] > 0){
      dummy_row <- data.frame("step"=0:(floor(unname(pars["t0"]))-1),"S"=population_n,"E"=0,"I"=0,"R"=0,"inc"=0)
      res <- bind_rows(dummy_row, res)
    }
    res <- res[res$step %in% times,]
    
    ## Get raw incidence and overall probability of infection
    incidence <- res$inc/population_n
    overall_prob <- max(res$R)/population_n
    
    if(!switch_model){
      res$beta <- pars["R0"]/pars["infectious"]
    }
    
    res$Rt <- (res$S/population_n) * res$beta * pars["infectious"]
  } 
  ## Get absolute incidence
  incidence <- incidence * population_n
  
  ## Reshape solution
  res_melted <- res %>% pivot_longer(-step)
  
  ## Compartment plot
  p_compartments <- res_melted %>% 
    filter(name %in% c("S","E","I","R","cumulative_incidence")) %>%
    ggplot() + geom_line(aes(x=step,y=value,col=name)) +
    ylab("Per capita") +
    xlab("Date") +
    theme_bw() +
    theme(legend.position="top")
  
  ## Incidence plot
  p_inc <- ggplot(data.frame(x=solve_times,y=incidence)) +
    geom_line(aes(x=x,y=y),col="red") +
    ylab("True incidence") +
    xlab("Date") +
    theme_bw()
  
  ## Rt plot
  p_rt <- res_melted %>% filter(name == "Rt") %>%
    ggplot() +
    geom_line(aes(x=step,y=value),col="blue") +
    scale_y_continuous(limits=c(0,pars["R0"]+1)) +
    geom_hline(yintercept=1,linetype="dashed") +
    ylab("Rt") +
    xlab("Date") +
    theme_bw()
  
  ## Combine plots
  inc_plot <- p_compartments / p_inc / p_rt
  
  ## Get growth rates
  GR_daily <- log(incidence[2:length(incidence)]/incidence[1:(length(incidence)-1)])
  GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
  GR_daily_dat <- data.frame(t=solve_times[2:length(solve_times)],GR=GR_daily,ver="daily")
  
  ## Get average growth rate over different size windows
  lastdays <- seq(10,50,by=10)
  
  GR_all <- GR_daily_dat
  for(t in seq_along(lastdays)){
    GR_full <- NULL
    lastday <- lastdays[t]
    for (i in (lastday+1):length(solve_times)) {
      end_index <- i-1
      start_index <- max(1, (i-lastday))
      GR_full <- c(GR_full,mean(GR_daily[start_index:end_index], na.rm=TRUE))
    }
    GR_full_dat <- data.frame(t=(lastday+1):length(solve_times), GR=GR_full,ver=as.character(lastday))
    GR_all <- bind_rows(GR_all, GR_full_dat)
  }
  
  ## Get daily growth rate around the peak
  gr_crossover <- GR_all %>% filter(ver == "daily") %>%
    filter(t < 250 & t > 100) %>%
    mutate(abs_gr = abs(GR)) %>%
    filter(abs_gr == min(abs_gr, na.rm=TRUE)) %>% pull(t)
  
  ## Average growth rates
  p_gr <- ggplot(GR_all %>% filter(ver != "daily")) +
    geom_line(aes(x=t,y=GR,col=ver)) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    xlab("Date") +    
    ggtitle("Average growth rate over different windows") +
    theme_bw()
  
  ## Daily growth rate
  p_gr1 <- ggplot(GR_all %>% filter(ver == "daily")) +
    geom_line(aes(x=t,y=GR),col="black") +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_vline(xintercept=gr_crossover,linetype="dotted")+
    coord_cartesian(ylim=c(-0.2,0.2)) +
    ylab("Growth rate") +
    ggtitle("Daily growth rate") +
    xlab("Date") +    
    theme_bw()
  
  
  list(seir_outputs=res, 
       incidence=incidence,
       overall_prob=overall_prob,
       plot=inc_plot, 
       growth_rates=GR_all,
       growth_rate_p =p_gr1/p_gr)
}
