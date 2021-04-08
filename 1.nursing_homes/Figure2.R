##########################################
# GENERATE FIGURE 2 
## - Nursing home fits, show trajectories for NH 1
## - Inset Ct distribution fits
## - Plot empirical relationship for all nursing homes
## - SEEIRR fit for NH 1
## - Compare growth rates for the different models
##########################################
## Get posterior summaries and simulated counts for inset
## Second item in list is from SEIR fit
run_use <- 2
summary_posterior_dat <- preds_fits_all_comb_models[[run_use]] %>% 
  left_join(obs_dat_models[[run_use]] %>% group_by(samp_t,location) %>% tally() %>% rename(loc=location)) %>%
  filter(ct < 40) %>%
  mutate(sim=rbinom(n(),n,density)) %>%
  mutate(density = density*n) %>%
  group_by(ct, samp_t, loc, ver) %>%
  summarize(lower=quantile(density,0.025),
            median=quantile(density,0.5),
            upper=quantile(density,0.975),
            lower_count=quantile(sim,0.025),
            median_count=quantile(sim,0.5),
            upper_count=quantile(sim,0.975))


summary_prob_detect <- preds_fits_all_comb_models[[run_use]] %>% 
  left_join(obs_dat_models[[run_use]] %>% group_by(samp_t,location) %>% tally() %>% rename(loc=location)) %>%
  filter(ct == 40) %>%
  mutate(prop_detectable=1-density) %>%
  group_by(t,loc, samp_t, ver) %>%
  summarize(med=median(prop_detectable),
            lower = quantile(prop_detectable,0.025),
            upper=quantile(prop_detectable,0.95)) %>%
  left_join(obs_dat_models[[run_use]] %>% group_by(samp_t,location) %>% summarize(dat_detectable=sum(ct<40)/n()) %>% rename(loc=location))




## Figure 2
use_loc <- "NH 1"
obs_dat <- obs_dat_models[[run_use]]
timepoints <- unique(obs_dat[obs_dat$location == use_loc,"samp_t"]) %>% pull(samp_t)

## Insets
inset1 <- plot_inset_fig2("NH 1", timepoints[1],obs_dat,summary_posterior_dat,summary_prob_detect,8,1)
inset2 <- plot_inset_fig2("NH 1", timepoints[2],obs_dat,summary_posterior_dat,summary_prob_detect,8,1)
inset3 <- plot_inset_fig2("NH 1", timepoints[3],obs_dat,summary_posterior_dat,summary_prob_detect,8,1)


p_prev_fits <- plot_inset_fig2_prev_combined("NH 1", summary_prob_detect,8,1) + labs(tag="D")

inset1_comb <- inset1[[1]] #+ inset1[[2]] + plot_layout(widths=c(3,1))
inset2_comb <- inset2[[1]]# + inset2[[2]] + plot_layout(widths=c(3,1))
inset3_comb <- inset3[[1]] #+ inset3[[2]] + plot_layout(widths=c(3,1))


remove_x_axis_theme <- theme(axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())

p_LHS <- 
  plot_spacer() +
  (inset1[[1]]+remove_x_axis_theme + labs(tag="D")) + #(inset1[[2]] + remove_x_axis_theme) +
  (inset2[[1]] + remove_x_axis_theme) +  #(inset2[[2]]+ remove_x_axis_theme) +
  (inset3[[1]]) +  #(inset3[[2]]) +
  plot_layout(ncol=1)


## Trajectories
p_predictions <- plot_predictions_fig2(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50) + labs(tag="C")
set.seed(1)
p_predictions1 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[1]) + 
  remove_x_axis_theme + theme(axis.title.y=element_blank()) + labs(tag="C")
p_predictions2 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[2]) + 
  remove_x_axis_theme+ theme(legend.position = "none")
p_predictions3 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[3]) + 
  theme(axis.title.y=element_blank(),legend.position="none",plot.margin = margin(0,0,0.0,0, "cm"),axis.title.x=element_blank())



## Inset fitted distributions
#p_predictions <- p_predictions +  inset_element(inset1_comb, 0.01,0.8,0.28,0.96) + inset_element(inset2_comb, 0.01,0.45,0.28,0.61) + inset_element(inset3_comb, 0.01,0.11,0.28,0.27)
#p_predictions <- p_predictions +  inset_element(inset1[[2]], 0.29,0.8,0.42,0.96) + inset_element(inset2[[2]], 0.29,0.45,0.42,0.61) + inset_element(inset3[[2]], 0.29,0.1,0.42,0.27)

## SEEIRR model fit
## NOTE - random_traj, quants, nh_prev1, best_traj, best_inc_traj, quants_inc and t0_pointrange are all generated in the SEEIRR script first
p_seeirr <- ggplot() +
  #geom_rect(data=tibble(xmin=as.Date("2020-01-01"),xmax=as.Date(min(gr_quants$time)),ymin=0,ymax=1),
  #          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  #geom_rect(data=tibble(xmin=as.Date(max(gr_quants$time)),xmax=as.Date("2020-05-16"),ymin=0,ymax=1),
  #          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  geom_line(data=random_traj %>% filter(location == use_loc),aes(x=time,y=prev/(0.7/0.12),group=samp),size=0.05,col="orange") +
  geom_ribbon(data=quants%>% filter(location == use_loc),
              aes(ymin=lower/(0.7/0.12),ymax=upper/(0.7/0.12),x=time),alpha=0.1,fill="orange") +
  geom_point(data=nh_prev1 %>%
               mutate(date=as.Date(date + min_date))%>% filter(location == use_loc), aes(x=date, y=prev/(0.7/0.12)),shape=8,
             col="black") +
  geom_errorbar(data=nh_prev1 %>%
                  mutate(date=as.Date(date + min_date))%>% filter(location == use_loc),
                aes(x=date,ymin=lower_confint/(0.7/0.12), ymax=upper_confint/(0.7/0.12)),
                width=2, col="black") +
  geom_line(data=best_traj%>% filter(location == use_loc), aes(x=time,y=prev/(0.7/0.12)),col="orange") +
  #geom_pointrange(data=t0_pointrange%>% filter(location == use_loc),aes(xmin=lower,xmax=upper,x=median),y=0.1, 
  #                col=AAAS_palette["blue1"],size=0.25,shape=18) +
  geom_ribbon(data=quants_inc%>% filter(location == use_loc),
              aes(ymin=lower,ymax=upper,
              #aes(ymin=lower*(0.7/0.12),ymax=upper*(0.7/0.12),
                  x=time),alpha=0.1,fill=AAAS_palette["red1"]) +
  geom_line(data=best_inc_traj%>% filter(location == use_loc), aes(x=time,
                                                                   y=inc),
            #y=inc*(0.7/0.12)),
            col=AAAS_palette["red1"]) +
  coord_cartesian(ylim=c(0,0.12),
                  xlim=as.Date(c("2020-03-01","2020-05-09"))) +
  #scale_y_continuous(expand=c(0,0), sec.axis=sec_axis(~./(0.7/0.12),name="PCR detectable prevalence",
  #                                                    breaks=seq(0,0.12,by=0.03))) +
  scale_x_date(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.12,by=0.03),sec.axis=sec_axis(~.*(0.7/0.12),name="PCR detectable prevalence")) +
  ylab("Per capita incidence") +
  xlab("Date") +
  export_theme +
  theme(strip.background = element_blank(),
        plot.margin = margin(0,0,0.1,0, "cm"),
        panel.grid.major = element_line(color="grey80",size=0.25),
        strip.text=element_text(face="bold",size=10),
        axis.text.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x = element_blank(),
        plot.tag.position = "topleft") +
  labs(tag="A")



p_growth_rates <- dat_inc %>% filter(loc == use_loc) %>%
  ggplot() +
  geom_hline(yintercept=0,col="grey50") +
  geom_ribbon(data=gr_quants%>% filter(loc == use_loc),aes(x=time,ymin=lower_gr,ymax=upper_gr),fill=AAAS_palette["purple3"],alpha=0.25) +
  geom_ribbon(data=gr_quants%>% filter(loc == use_loc),aes(x=time,ymin=lower_gr_avg,ymax=upper_gr_avg),fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=gr_quants%>% filter(loc == use_loc),aes(x=time,y=median_gr),col=AAAS_palette["purple3"]) +
  geom_line(data=gr_quants%>% filter(loc == use_loc),aes(x=time,y=median_gr_avg),col=AAAS_palette["green1"]) +
  #scale_x_date(limits=range(gr_quants$time)) +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "seir_FALSE_FALSE")%>% filter(loc == use_loc),
              aes(x=t+1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["purple3"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "exp_FALSE_FALSE")%>% filter(loc == use_loc),
              aes(x=t-1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["green1"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  coord_cartesian(ylim=c(-0.55,0.55),
                  xlim=as.Date(c("2020-03-01","2020-05-09"))) +
  #scale_y_continuous(expand=c(0,0), sec.axis=sec_axis(~./(0.7/0.12),name="PCR detectable prevalence",
  #                                                    breaks=seq(0,0.12,by=0.03))) +
  scale_x_date(expand=c(0,0)) +
  export_theme +
  ylab("Growth rate") +
  xlab("Date") +
  theme(strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.major = element_line(color="grey80",size=0.25),
        strip.background = element_blank()) +
  labs(tag="E")



fig2 <- (plot_spacer() + theme(plot.margin = margin(0,0,0,0, "cm"))) +
  (inset1[[1]]+remove_x_axis_theme + labs(tag="B") + theme(axis.title.y=element_blank())) + 
  (inset2[[1]] + remove_x_axis_theme) +  
  (inset3[[1]] + theme(axis.title.y=element_blank())) + 
  (p_prev_fits + theme(plot.margin = margin(0,0,0,0, "cm"))) +
  p_seeirr + p_predictions1 + p_predictions2 + p_predictions3 + p_growth_rates +
  plot_layout(ncol=2,byrow=FALSE,widths=c(1,2.25),heights=c(1,1,1,1,1.25))


p_prev <- obs_dat %>% group_by(location, mean_week_date) %>%
  mutate(pos=ct < 40) %>%
  summarize(n=n(),pos_all=sum(pos)) %>%
  mutate(prev=pos_all/n) %>%
  group_by(location, mean_week_date) %>%
  mutate(lower_confint = prop.test(pos_all, n)$conf.int[1],
         upper_confint=prop.test(pos_all,n)$conf.int[2]) %>%
  rename(Location=location) %>%
  ggplot() +
  geom_point(aes(x=mean_week_date,y=prev,group=Location,col=Location),size=2) +
  geom_line(aes(x=mean_week_date,y=prev,group=Location,col=Location),size=1) +
  geom_errorbar(aes(x=mean_week_date,ymin=lower_confint,ymax=upper_confint,group=Location,col=Location)) +
  scale_color_viridis_d() +
  scale_x_date(limits=range(gr_quants$time)) +
  export_theme  +
  ylab("Proportion of samples positive") +
  xlab("Swab collection date") +
  labs(tag="B") +
  theme(legend.position=c(0.9,0.9),axis.text.x=element_blank(),axis.title.x=element_blank())

p_median <- dat_properties %>% 
  ggplot() +
  geom_point(aes(x=mean_week_date,y=ct_median,group=Location,col=Location),size=2) +
  geom_line(aes(x=mean_week_date,y=ct_median,group=Location,col=Location),size=1) +
  scale_color_viridis_d() +
  scale_x_date(limits=range(gr_quants$time)) +
  scale_y_continuous(trans="reverse") +
  export_theme  +
  ylab("Median of N2 Cts") +
  xlab("Swab collection date") +
  labs(tag="C") +
  theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())
p_skew <- dat_properties %>% 
  ggplot() +
  geom_point(aes(x=mean_week_date,y=ct_skew,group=Location,col=Location),size=2) +
  geom_line(aes(x=mean_week_date,y=ct_skew,group=Location,col=Location),size=1) +
  scale_x_date(limits=range(gr_quants$time)) +
  scale_color_viridis_d() +
  export_theme  +
  ylab("Skewness of N2 Cts") +
  xlab("Swab collection date") +
  labs(tag="D") +
  theme(legend.position="none")

p_nh_relationships <- p_prev/p_median/p_skew + plot_layout()

#fig2 <- p_seeirr + p_predictions + p_growth_rates + p_RHS + plot_layout(byrow=FALSE, widths=c(1.75,1),heights=c(1,3))


#p_lhs <-  p_predictions/p_seeirr + plot_layout(heights=c(2.5,1))
#p_rhs <- ((p_nh_relationships)+p_growth_rates) + plot_layout(heights=c(1,1,1,1))
#fig2 <- (p_lhs| p_rhs) + plot_layout(widths=c(2,1))



if(TRUE){
  ggsave("figures/Figure2.pdf",fig2,width=7,height=8)
  ggsave("figures/Figure2.png",fig2,width=7,height=8)
}

