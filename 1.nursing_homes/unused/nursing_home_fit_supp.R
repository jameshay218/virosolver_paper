use_loc <- "NH 1"
use_pars <- c("viral_peak","obs_sd","level_switch")

## Plot predicted trajectories
p_predictions <- ggplot(preds_all_models[[2]] %>% filter(loc == use_loc) %>% filter(sampno <= 100)) + 
  geom_line(aes(x=t,y=40-prob_infection/0.004,group=paste0(sampno,samp_t),col=as.factor(samp_t)),size=0.1) +
  #geom_line(data=best_prob_models[[2]]%>% filter(loc == use_loc),aes(x=t,y=40-prob_infection/0.004,col="MAP"),size=0.5) +
  geom_violin(data=obs_dat_models[[2]]%>% filter(location == use_loc),aes(x=mean_week_date, y=ct,group=mean_week_date),
              width=8,scale="width",fill=AAAS_palette[1],alpha=0.25,col="grey10",
              draw_quantiles=c(0.025,0.5,0.975)) +
  geom_jitter(data=obs_dat_models[[2]]%>% filter(location == use_loc),aes(x=mean_week_date, y=ct),size=0.25,height=0,width=1.25,col=AAAS_palette[1]) +
  
  scale_y_continuous(trans="reverse",
                     sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                       breaks=seq(0,0.12,by=0.03))) +
  xlab("Date") +
  ylab("Ct value") +
  export_theme + 
  #scale_color_manual(values=c("Posterior draw"="#EE0000FF","MAP"="#008B45FF")) +
  scale_linetype_manual(values=c("Sample date"="dashed")) +
  guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
  theme(legend.position=c(0.1,0.8),
        plot.margin = margin(0,0,0,0, "cm")) +
  #facet_wrap(~samp_t, nrow=1)+
  labs(tag="A")

summary_posterior_dat <- pred_fits_all_comb %>% 
  filter(loc == use_loc) %>%
  left_join(obs_dat %>% filter(location == use_loc) %>% group_by(samp_t) %>% tally()) %>%
  mutate(sim=rbinom(n(),n,density)) %>%
  mutate(density = density*n) %>%
  group_by(ct, samp_t, loc, ver) %>%
  summarize(lower=quantile(density,0.025),
            median=quantile(density,0.5),
            upper=quantile(density,0.975),
            lower_count=quantile(sim,0.025),
            median_count=quantile(sim,0.5),
            upper_count=quantile(sim,0.975))

## Plot fits to distributions
p_fits <- ggplot(obs_dat_models[[2]] %>% filter(location == use_loc)) +
  geom_histogram(aes(x=ct),binwidth=1,fill="grey70",col="grey20",boundary=0) +
  geom_ribbon(data=summary_posterior_dat %>% filter(loc == use_loc),aes(x=ct+0.5,ymin=lower_count,ymax=upper_count),fill="blue",alpha=0.1)+
  geom_ribbon(data=summary_posterior_dat %>% filter(loc == use_loc),aes(x=ct+0.5,ymin=lower,ymax=upper),fill="blue",alpha=0.25)+
  geom_line(data=summary_posterior_dat %>% filter(loc == use_loc),
            aes(x=ct+0.5,y=median)) +
  facet_wrap(~samp_t) +
  scale_x_continuous(trans="reverse",expand=c(0,0),limits=c(39,5),breaks=seq(0,40,by=5)) +
  coord_cartesian(xlim=c(39,5),ylim=c(0,8)) +
  scale_y_continuous(breaks=seq(0,8,by=2),expand=c(0,0)) +
  export_theme +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        panel.grid.minor = element_line(color="grey80",size=0.25),
        plot.margin = margin(0,0,0,0, "cm")) +
  xlab("Ct value") +
  ylab("Count") +
  labs(tag="B")

## Proportion detectable
summary_prop_detectable <- preds_fits_all_comb_models[[2]] %>% filter(ct == 40) %>%
  filter(loc == use_loc) %>%
  group_by(t) %>%
  mutate(density=1-density) %>%
  summarize(lower=quantile(density,0.025),
            median=quantile(density,0.5),
            upper=quantile(density,0.975)) %>%
  rename(samp_t = t) %>%
  mutate(samp_t = as.Date(samp_t,origin=min_date))

p_fit_detect <-  ggplot(obs_dat_models[[2]] %>% filter(location == "NH 1") %>%
                group_by(samp_t) %>%
                mutate(is_detectable=ct < 40) %>%
                summarize(prop_detectable=sum(is_detectable)/n())) +
  geom_point(aes(y=prop_detectable,x=0.4,col="Ground truth"),size=3,shape=18) +
  geom_point(data=summary_prop_detectable,aes(x=0.5,y=median,col="Posterior median & 95% CI"),size=1) +
  geom_errorbar(data=summary_prop_detectable,aes(x=0.5,ymin=lower,ymax=upper),
                width=0.025, col="blue") +
  scale_y_continuous() +
  scale_x_continuous(limits=c(0.3,0.6)) +
  scale_color_manual(values=c("Ground truth"="grey40",
                              "Posterior median & 95% CI"="blue")) +
  guides(color=guide_legend(title=NULL)) +
  facet_wrap(~samp_t,nrow=1) +
  ylab("Proportion detectable") +
  xlab("") +
  export_theme +
  theme(panel.grid.major = element_line(color="grey80",size=0.25),
        panel.grid.minor = element_line(color="grey80",size=0.25),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x=element_blank(),
        legend.background = element_rect(fill="white"),
        legend.margin = margin(0.1,0.1,0.1,0.1, "cm"),
        legend.position=c(0.15,0.3),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.ticks.x=element_blank()) +
  labs(tag="C")

par_key <- c("viral_peak" = "Ct[peak]",
             "obs_sd"="sigma",
             "level_switch"="Ct[switch]")

## Plot parameter posterior densities
chains_all_models[[2]]$variable <- as.character(chains_all_models[[2]]$variable)
chains_all_models[[2]]$t <- as.character(chains_all_models[[2]]$t)
p_densities <-  chains_all_models[[2]] %>% 
  filter(loc == use_loc) %>%
  filter(variable %in% use_pars) %>%
  mutate(variable = par_key[variable]) %>%
  ggplot() +
  geom_density(aes(x=value),alpha=0.25,fill="grey50") +
  facet_wrap(t~variable,scales="free",nrow=1,labeller = label_parsed) +
  ylab("Posterior density") +
  xlab("Value") +
  export_theme +
  theme(strip.text=element_text(size=6),
        plot.margin = margin(0,0,0,0, "cm"),
        axis.text.x=element_text(size=4)) +
  labs(tag="D")

fig <- (p_predictions/p_fits/p_fit_detect/p_densities) + plot_layout(heights=c(2,2,2,1.5))

ggsave("figures/supplement/nursing_home_seir_fits.pdf",height=7,width=8)
ggsave("figures/supplement/nursing_home_seir_fits.png",height=7,width=8,units="in",dpi=300)

