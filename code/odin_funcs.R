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
