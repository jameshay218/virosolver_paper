##########################################
# GENERATE supplementary figure for nursing home fits
## - Nursing home fits, show trajectories for NH 1
## - Inset Ct distribution fits
## - Plot empirical relationship for all nursing homes
## - SEEIRR fit for NH 1
## - Compare growth rates for the different models
##########################################
## Get posterior summaries and simulated counts for inset
## Second item in list is from SEIR fit

remove_x_axis_theme <- theme(axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())
remove_y_axis_theme <- theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y = element_blank())

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

## NH 2
use_loc <- "NH 2"
obs_dat <- obs_dat_models[[run_use]]
timepoints <- unique(obs_dat[obs_dat$location == use_loc,"samp_t"]) %>% pull(samp_t)
## Insets
inset1_nh2 <- plot_inset_fig2(use_loc, timepoints[1],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset2_nh2 <- plot_inset_fig2(use_loc, timepoints[2],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset3_nh2 <- plot_inset_fig2(use_loc, timepoints[3],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)

p_predictions1_NH2 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[1]) + 
  remove_x_axis_theme + theme(axis.title.y=element_blank(),legend.position=c(0.8,0.8)) + labs(tag="B")
p_predictions2_NH2 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[2]) + 
  remove_x_axis_theme+ theme(legend.position = "none")
p_predictions3_NH2 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[3]) + 
  theme(axis.title.y=element_blank(),legend.position="none",plot.margin = margin(unit(c(0.1,0.1,0.1,0.1,0.1), "cm")),axis.title.x=element_blank()) 

p_predictions1_NH2_alt <- p_predictions1_NH2 + inset_element(inset1_nh2[[1]] + theme(axis.title.x=element_blank()), 0.01,0.45,0.5,0.95) + inset_element(inset1_nh2[[2]], 0.7,0.4,1,0.85)
p_predictions2_NH2_alt <- p_predictions2_NH2 + inset_element(inset2_nh2[[1]]+ theme(axis.title.x=element_blank()), 0.01,0.45,0.5,0.95) + inset_element(inset2_nh2[[2]], 0.7,0.4,1,0.85)
p_predictions3_NH2_alt <- p_predictions3_NH2 + inset_element(inset3_nh2[[1]]+ theme(axis.title.x=element_blank()+ remove_x_axis_theme), 0.01,0.45,0.5,0.95) + inset_element(inset3_nh2[[2]]+ theme(plot.background = element_blank()), 0.7,0.4,1,0.85) 


## NH 3
use_loc <- "NH 3"
obs_dat <- obs_dat_models[[run_use]]
timepoints <- unique(obs_dat[obs_dat$location == use_loc,"samp_t"]) %>% pull(samp_t)
## Insets
inset1_nh3 <- plot_inset_fig2(use_loc, timepoints[1],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset2_nh3 <- plot_inset_fig2(use_loc, timepoints[2],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset3_nh3 <- plot_inset_fig2(use_loc, timepoints[3],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)

p_predictions1_NH3 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[1]) + 
  remove_x_axis_theme + theme(axis.title.y=element_blank(),legend.position="none") + remove_y_axis_theme
p_predictions2_NH3 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[2]) + 
  remove_x_axis_theme+ theme(axis.title.y=element_blank(),legend.position = "none") + remove_y_axis_theme
p_predictions3_NH3 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[3]) + 
  theme(axis.title.y=element_blank(),legend.position="none",plot.margin = margin(unit(c(0.1,0.1,0.1,0.1,0.1), "cm")),axis.title.x=element_blank()) + remove_y_axis_theme


p_predictions1_NH3_alt <- p_predictions1_NH3 + inset_element(inset1_nh3[[1]]+ theme(axis.title.x=element_blank()), 0.01,0.45,0.5,0.95) + inset_element(inset1_nh3[[2]], 0.7,0.4,1,0.85)
p_predictions2_NH3_alt <- p_predictions2_NH3 + inset_element(inset2_nh3[[1]]+ theme(axis.title.x=element_blank()), 0.01,0.45,0.5,0.95) + inset_element(inset2_nh3[[2]], 0.7,0.4,1,0.85)
p_predictions3_NH3_alt <- p_predictions3_NH3 + inset_element(inset3_nh3[[1]]+ theme(axis.title.x=element_blank()+ remove_x_axis_theme), 0.01,0.45,0.5,0.95) + inset_element(inset3_nh3[[2]], 0.7,0.4,1,0.85)


## NH 4
use_loc <- "NH 4"
obs_dat <- obs_dat_models[[run_use]]
timepoints <- unique(obs_dat[obs_dat$location == use_loc,"samp_t"]) %>% pull(samp_t)
## Insets
inset1_nh4 <- plot_inset_fig2(use_loc, timepoints[1],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset2_nh4 <- plot_inset_fig2(use_loc, timepoints[2],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)
inset3_nh4 <- plot_inset_fig2(use_loc, timepoints[3],obs_dat,summary_posterior_dat,summary_prob_detect,10,1)

p_predictions1_NH4 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[1]) + 
  theme(axis.title.y=element_blank(),legend.position="none") + remove_y_axis_theme+ remove_x_axis_theme
p_predictions2_NH4 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[2]) + 
  theme(axis.title.y=element_blank(),legend.position = "none") + remove_y_axis_theme+ remove_x_axis_theme
p_predictions3_NH4 <- plot_predictions_fig2_indiv(use_loc,preds_all_models[[run_use]],best_prob_models[[run_use]],obs_dat %>% filter(ct < 40),50,timepoint=timepoints[3]) + 
  theme(axis.title.y=element_blank(),legend.position="none",plot.margin = margin(unit(c(0.1,0.1,0.1,0.1,0.1), "cm")),axis.title.x=element_blank()) + 
  remove_y_axis_theme


p_predictions1_NH4_alt <- p_predictions1_NH4 + inset_element(inset1_nh4[[1]]+ theme(axis.title.x=element_blank() ), 0.01,0.45,0.5,0.95) + inset_element(inset1_nh4[[2]], 0.7,0.4,1,0.85)
p_predictions2_NH4_alt <- p_predictions2_NH4 + inset_element(inset2_nh4[[1]]+ theme(axis.title.x=element_blank() ), 0.01,0.45,0.5,0.95) + inset_element(inset2_nh4[[2]], 0.7,0.4,1,0.85)
p_predictions3_NH4_alt <- p_predictions3_NH4 + inset_element(inset3_nh4[[1]]+ theme(axis.title.x=element_blank() + remove_x_axis_theme), 0.01,0.45,0.5,0.95) + inset_element(inset3_nh4[[2]], 0.7,0.4,1,0.85)


main_p_alt <- p_predictions1_NH2_alt + p_predictions2_NH2_alt + p_predictions3_NH2_alt +
  p_predictions1_NH3_alt + p_predictions2_NH3_alt + p_predictions3_NH3_alt +
  p_predictions1_NH4_alt + p_predictions2_NH4_alt + p_predictions3_NH4_alt + plot_layout(byrow=FALSE)

no_y_theme <- theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.line.y=element_blank(),axis.ticks.y=element_blank()) 
main_p_top <- ((inset1_nh2[[1]]+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2)) + xlab("") + ggtitle("2020-04-09") + labs(tag="A")) +  
  (inset2_nh2[[1]] + no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ ggtitle("2020-04-17")) +
  (inset3_nh2[[1]] + no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ xlab("")+ ggtitle("2020-05-04")) +
  (inset1_nh3[[1]] +  scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ xlab("")+ ggtitle("2020-04-09")+ labs(tag="B")) +
  (inset2_nh3[[1]] +  no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ ggtitle("2020-04-16")) +
  (inset3_nh3[[1]] + no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ xlab("")+ ggtitle("2020-05-04")) +
  (inset1_nh4[[1]] +  scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ xlab("")+ ggtitle("2020-04-10")+ labs(tag="C")) +
  (inset2_nh4[[1]] +  no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ ggtitle("2020-04-16")) +
  (inset3_nh4[[1]] + no_y_theme+ scale_y_continuous(expand=c(0,0),limits=c(0,10.5),breaks=seq(0,10,by=2))+ xlab("")+ ggtitle("2020-05-05")) +
     plot_layout(nrow=1))

main_p <- ((p_predictions1_NH2 + labs(tag="D")) +
  p_predictions2_NH2 + 
  p_predictions3_NH2 +
  (p_predictions1_NH3  + labs(tag="E")) + 
  p_predictions2_NH3 + 
  (p_predictions3_NH3+ xlab("Date") + theme(axis.title.x=element_text(size=10,family="sans"))) +
  (p_predictions1_NH4  + labs(tag="F"))+
  p_predictions2_NH4 +
  p_predictions3_NH4 + 
  plot_layout(byrow=FALSE))

figSY_new <- (main_p_top / main_p) + plot_layout(heights=c(1,5))

## SEEIRR model fit
## NOTE - random_traj, quants, nh_prev1, best_traj, best_inc_traj, quants_inc and t0_pointrange are all generated in the SEEIRR script first
p_seeirr_supp <- ggplot() +
  #geom_rect(data=tibble(xmin=as.Date("2020-01-01"),xmax=as.Date(min(gr_quants$time)),ymin=0,ymax=1),
  #          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  #geom_rect(data=tibble(xmin=as.Date(max(gr_quants$time)),xmax=as.Date("2020-05-16"),ymin=0,ymax=1),
  #          aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  geom_ribbon(data=quants%>% filter(location != "NH 1"),
              aes(ymin=lower/(0.7/0.12),ymax=upper/(0.7/0.12),x=time),alpha=0.1,fill="orange") +
  geom_line(data=random_traj %>% filter(location != "NH 1"),aes(x=time,y=prev/(0.7/0.12),group=samp),size=0.05,col="orange") +
  geom_line(data=best_traj%>% filter(location != "NH 1"), aes(x=time,y=prev/(0.7/0.12)),col="orange") +
  geom_point(data=nh_prev1 %>%
               mutate(date=as.Date(date + min_date))%>% filter(location != "NH 1"), aes(x=date, y=prev/(0.7/0.12)),shape=8,
             col=AAAS_palette["black"]) +
  geom_errorbar(data=nh_prev1 %>%
                  mutate(date=as.Date(date + min_date))%>% filter(location != "NH 1"),
                aes(x=date,ymin=lower_confint/(0.7/0.12), ymax=upper_confint/(0.7/0.12)),
                width=2, col=AAAS_palette["black"]) +
  #geom_pointrange(data=t0_pointrange%>% filter(location != "NH 1"),aes(xmin=lower,xmax=upper,x=median),y=0.1, 
  #                col=AAAS_palette["blue1"],size=0.25,shape=18) +
  geom_ribbon(data=quants_inc%>% filter(location != "NH 1"),
              aes(ymin=lower,ymax=upper,x=time),alpha=0.1,fill=AAAS_palette["red1"]) +
  geom_line(data=best_inc_traj%>% filter(location != "NH 1"), aes(x=time,y=inc),col=AAAS_palette["red1"]) +
  
  coord_cartesian(xlim=as.Date(c("2020-03-01","2020-05-09")),ylim=c(0,0.12)) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.12,by=0.03),sec.axis=sec_axis(~.*(0.7/0.12),name="PCR detectable prevalence")) +
  scale_x_date(expand=c(0,0)) +
  ylab("Per capita incidence") +
  xlab("Date") +
  export_theme +
  theme(strip.background = element_blank(),
        strip.text=element_text(face="bold",size=10),
        panel.grid.major=element_line(size=0.1,color="grey40"),
        plot.tag.position = "topleft") +
  remove_x_axis_theme +
  labs(tag="A") +
  facet_wrap(~location)

p_growth_rates_supp <- dat_inc %>% filter(loc != "NH 1") %>%
  ggplot() +
  geom_hline(yintercept=0,col="grey50") +
  geom_rect(data=data.frame(ymin=-0.6,ymax=0.6,xmin=as.Date("2020-03-01"),xmax=min(gr_quants$time)),aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax),fill="grey50",alpha=0.25) +
  geom_ribbon(data=gr_quants%>% filter(loc != "NH 1"),aes(x=time,ymin=lower_gr,ymax=upper_gr),fill=AAAS_palette["purple3"],alpha=0.25) +
  geom_ribbon(data=gr_quants%>% filter(loc != "NH 1"),aes(x=time,ymin=lower_gr_avg,ymax=upper_gr_avg),fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=gr_quants%>% filter(loc != "NH 1"),aes(x=time,y=median_gr),col=AAAS_palette["purple3"]) +
  geom_line(data=gr_quants%>% filter(loc != "NH 1"),aes(x=time,y=median_gr_avg),col=AAAS_palette["green1"]) +
  scale_x_date(limits=as.Date(c("2020-03-01","2020-05-09")),expand=c(0,0)) +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "seir_FALSE_FALSE")%>% filter(loc != "NH 1"),
              aes(x=t+1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["purple3"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "exp_FALSE_FALSE")%>% filter(loc != "NH 1"),
              aes(x=t-1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["green1"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  coord_cartesian(ylim=c(-0.55,0.55)) +
  export_theme +
  ylab("Growth rate") +
  xlab("Date") +
  theme(strip.text = element_blank(),
        panel.grid.major=element_line(size=0.1,color="grey40"),
        strip.background = element_blank()) +
  labs(tag="B") + facet_wrap(~loc)

figSX_new <- p_seeirr_supp/p_growth_rates_supp

#figSX <- p_seeirr_supp/main_p/p_growth_rates_supp + plot_layout(heights=c(1,3,1))

if(TRUE){
  ggsave("figures/supplement/FigureS4_a.pdf",figSX_new,width=7,height=4)
  ggsave("figures/supplement/FigureS4_a.png",figSX_new,width=7,height=4,dpi=300,units="in")
  
  
  ggsave("figures/supplement/FigureS4_b.pdf",figSY_new,width=8,height=7)
  ggsave("figures/supplement/FigureS4_b.png",figSY_new,width=8,height=7,dpi=300,units="in")
}

