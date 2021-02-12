#seir_pars <- c("R0"=R0,"infectious"=9,"incubation"=6,"I0"=0.0001,"t0"=0)
#seir_pars <- c("R0"=R0,"infectious"=7,"incubation"=4,"I0"=0.001,"t0"=0)
#names(seir_pars) <- c("R0","infectious","incubation","I0","t0")
#epidemic_process <- simulate_seir_process(seir_pars,times,population_n)

## Simulate onset times, confirmation delays etc
observed_individuals <- simulate_observations_wrapper(incidence=epidemic_process$raw_incidence,times=0:365,
                                                      population_n=population_n)

## Random surveillance
observed_indivs_flat <- simulate_reporting(observed_individuals, frac_report=0.2,timevarying_prob=NULL,
                                           solve_times=times, symptomatic=FALSE)

simulated_viral_loads <- simulate_viral_loads_wrapper(observed_indivs_flat$sampled_individuals,
                                                      kinetics_pars=pars)


seir_weekly <- epidemic_process$seir_outputs %>% 
  pivot_wider(names_from=variable,values_from=value) %>%
  as_tibble() %>% 
  rename(sampled_time=time) %>%
  mutate(week=floor(sampled_time/7))%>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  group_by(week, first_day) %>%
  summarize(mean_rt=mean(Rt))

combined_dat <- simulated_viral_loads %>% 
  mutate(week=floor(sampled_time/7))%>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < pars["intercept"]) %>%
  #filter(sampled_time >= min(use_days) & sampled_time <= max(use_days)) %>%
  #group_by(sampled_time) %>% 
  #filter(ct_obs < kinetics_pars["intercept"]) %>% 
  group_by(week) %>%
  summarize(median_ct=median(ct_obs),
            skew_ct=moments::skewness(ct_obs),
            n=n()) %>%
  #rename(date=first_day) %>%
  full_join(seir_weekly) %>%
  filter(n >= 25)


p1 <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=skew_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=skew_ct),size=2) +
  geom_vline(xintercept=1,linetype="dashed") +
  #scale_y_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  main_theme +
  ylab("Skewness of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="A")

p2 <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=median_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=median_ct),size=2)+
  geom_vline(xintercept=1,linetype="dashed") +
  main_theme +
  #scale_y_continuous(trans="reverse",limits=c(35,28),breaks=seq(28,35,by=1)) +
  scale_y_continuous(trans="reverse") +
  ylab("Median of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="B")

p3 <- ggplot(combined_dat) +
  geom_point(aes(x=skew_ct,y=median_ct,col=mean_rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1,
                        limits=c(0,2),
                        guide=guide_colorbar(title="Rt",
                                             barwidth=1,ticks=FALSE,barheight=4))+
  #scale_x_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  scale_y_continuous(trans="reverse")+#,limits=c(35,28),breaks=seq(28,35,by=1)) +
  main_theme +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  theme(legend.position=c(0.2,0.8),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6),
        panel.grid.major=element_line(colour="grey40",size=0.1)) +
  labs(tag="H")


p1 | p2 | p3
