#### Running SEIR Model: ####
set.seed(0)
## Entire population size
population_n <- 6900000
## Sample size
sample_n <- 500
## Duration of epidemic in days
times <- 0:365

seir_pars <- c("R0"=R0,"infectious"=4,"incubation"=4,"I0"=0.0001,"t0"=0)
names(seir_pars) <- c("R0","infectious","incubation","I0","t0")
epidemic_process <- simulate_seir_process(seir_pars,times,population_n)

x <- seq(0,40,by=1)
lastday <- 35
incidence <- epidemic_process$raw_incidence[1:200]/population_n

#### Calculating Growth Rates: ####
GR_daily <- log(incidence[2:200]/incidence[1:199])
GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
GR_full <- NULL
GR_25 <- NULL
for (i in (lastday+1):length(GR_daily)) {
  GR_full <- c(GR_full,mean(GR_daily[(i-lastday):(i-1)], na.rm=TRUE))
  GR_25 <- c(GR_25,mean(GR_daily[(i-25):(i-1)], na.rm=TRUE))
}
GR <- data.frame(Days=1:length(GR_daily), GR=c(rep(NA,lastday),GR_full), GR25=c(rep(NA,lastday),GR_25), DailyGR=GR_daily)

#### Calculating viral load and age distributions: ####
Ct_res <- matrix(0,nrow=length(incidence), ncol=length(x))
age_res <- matrix(0,nrow=length(incidence),ncol=lastday)

# Detectable probabilities over [lastday] days prior to test
viral_loads <- viral_load_func(pars, 1:lastday)
detectable_props <- prop_detectable(1:lastday, pars, viral_loads)
ages <- 1:lastday
## Time at which standard deviation is reduced
t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
sd_mod <- rep(pars["sd_mod"], max(ages))

## Prior to t_switch, full standard deviation
## Max sure we don't go past maximum of ages vector
unmod_vec <- 1:min(t_switch,max(ages)+1)
sd_mod[unmod_vec] <- 1

## For the next sd_mod_wane days, decrease linearly
decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])

# Create matrices of relative frequencies of observing age of infection and viral load by day of testing:
dat <- pred_dist_wrapper(x, 1:200,1:35,pars, incidence)

for (i in 2:(dim(Ct_res)[1])) {
  print(i)
  past_inc <- incidence[(i-1):(max(i-lastday,1))]
  days <- 1:length(past_inc)
  age_res[i,days] <- past_inc*detectable_props[days]
}


# Summaries of age and VL distributions:
age_res_std <- age_res/apply(age_res, 1, sum, na.rm=TRUE)
age_mean <- apply(age_res_std, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))
age_res_std_csum <- t(apply(age_res_std, 1, FUN=function(res) cumsum(res)))
age_median <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.5]))
plot(x=0:199, y=age_mean, type="l", lwd=2, col="blue")
lines(x=0:199, y=age_median, lwd=2, col="green", lty=2)

## Plot mean of detectable Cts
dat %>% filter(ct < 40) %>% group_by(t) %>% 
  mutate(density_scaled=density/sum(density)) %>% 
  summarize(y=sum(ct*density_scaled)) %>% 
  ggplot() + geom_line(aes(x=t,y=y))


dat %>% filter(ct < 40) %>% group_by(t) %>% 
  mutate(density_scaled=density/sum(density)) %>% 
  mutate(cumu_density=cumsum(density_scaled)) %>%
  filter(cumu_density >= 0.5) %>%
  filter(ct == min(ct)) %>%
  ggplot() + geom_line(aes(x=t,y=ct))

