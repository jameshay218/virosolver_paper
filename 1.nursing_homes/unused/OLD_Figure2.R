
################################
## Additional plots
################################
## Prevalence over time
p_prev1 <- ggplot() +
  geom_rect(data=tibble(xmin=as.Date("2020-01-01"),xmax=as.Date(min(gr_quants$time)),ymin=0,ymax=1),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  geom_rect(data=tibble(xmin=as.Date(max(gr_quants$time)),xmax=as.Date("2020-05-16"),ymin=0,ymax=1),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="grey70",alpha=0.2) +
  geom_line(data=random_traj,aes(x=time,y=prev,group=samp),size=0.05,col=AAAS_palette["black"]) +
  geom_ribbon(data=quants,
              aes(ymin=lower,ymax=upper,x=time),alpha=0.1,fill=AAAS_palette["black"]) +
  geom_point(data=nh_prev1 %>%
               mutate(date=as.Date(date + min_date)), aes(x=date, y=prev),
             col=AAAS_palette["red2"]) +
  geom_errorbar(data=nh_prev1 %>%
                  mutate(date=as.Date(date + min_date)),
                aes(x=date,ymin=lower_confint, ymax=upper_confint),
                width=2, col=AAAS_palette["red2"]) +
  geom_line(data=best_traj, aes(x=time,y=prev),col=AAAS_palette["black"]) +
  geom_pointrange(data=t0_pointrange,aes(xmin=lower,xmax=upper,x=median),y=0.1, 
                  col=AAAS_palette["blue1"],size=0.25,shape=18) +
  geom_ribbon(data=quants_inc,
              aes(ymin=lower*(0.7/0.12),ymax=upper*(0.7/0.12),x=time),alpha=0.1,fill=AAAS_palette["blue1"]) +
  geom_line(data=best_inc_traj, aes(x=time,y=inc*(0.7/0.12)),col=AAAS_palette["blue1"]) +
  
  coord_cartesian(ylim=c(0,0.7),xlim=as.Date(c("2020-02-01","2020-05-09"))) +
  scale_y_continuous(expand=c(0,0), sec.axis=sec_axis(~./(0.7/0.12),name="Per capita incidence",
                                                      breaks=seq(0,0.12,by=0.03))) +
  ylab("PCR detectable prevalence") +
  xlab("Date") +
  export_theme +
  theme(axis.title.x=element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(face="bold",size=10),
        plot.tag.position = "topleft") +
  facet_wrap(~location,nrow=1) +
  labs(tag="A")

p_gr <- dat_inc %>% 
  ggplot() +
  geom_hline(yintercept=0,col="grey50") +
  geom_ribbon(data=gr_quants,aes(x=time,ymin=lower_gr,ymax=upper_gr),fill=AAAS_palette["purple3"],alpha=0.25) +
  geom_ribbon(data=gr_quants,aes(x=time,ymin=lower_gr_avg,ymax=upper_gr_avg),fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=gr_quants,aes(x=time,y=median_gr),col=AAAS_palette["purple3"]) +
  geom_line(data=gr_quants,aes(x=time,y=median_gr_avg),col=AAAS_palette["green1"]) +
  scale_x_date(limits=range(gr_quants$time), breaks="7 days") +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "seir_FALSE_FALSE"),
              aes(x=t+1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["purple3"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  geom_violin(data=betas_all_models_comb %>% filter(ver == "exp_FALSE_FALSE"),
              aes(x=t-1, y=beta,group=t),scale="width",width=2,fill=AAAS_palette["green1"],
              trim=TRUE,alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),size=0.1) +
  coord_cartesian(ylim=c(-0.55,0.55)) +
  export_theme +
  ylab("Growth rate") +
  xlab("Date") +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,size=7)) +
  facet_wrap(~loc,nrow=1) +
  labs(tag="C")

## Incidence over time with Ct values
p_inc <- ggplot() +
  geom_violin(data=dat_use,aes(x=mean_week_date, y=ct,group=mean_week_date),
              width=4,scale="width",fill="grey70",alpha=0.25,col="grey70",
              draw_quantiles=c(0.025,0.5,0.975)) +
  geom_jitter(data=dat_use,aes(x=mean_week_date, y=ct),size=0.25,height=0,width=1.25) +
  geom_line(data=dat_use %>% group_by(location,mean_week_date) %>% summarize(y=median(ct, na.rm=TRUE)),
            aes(x=mean_week_date,y=y),col=AAAS_palette["grey"],linetype="longdash",size=0.25) +
  geom_ribbon(data=quants_inc,
              aes(ymin=40 - lower/0.004,ymax=40 - upper/0.004,x=time),
              alpha=0.1,fill=AAAS_palette["blue1"]) +
  geom_line(data=best_inc_traj, 
            aes(x=time,y=40 - inc/0.004),
            col=AAAS_palette["blue1"]) +
  
  scale_x_date(limits=range(gr_quants$time), breaks="7 days") +
  scale_y_continuous(trans="reverse",
                     sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                       breaks=seq(0,0.12,by=0.03))) +
  ylab("Ct value") +
  xlab("Date") +
  export_theme +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_wrap(~location,nrow=1) +
  labs(tag="B")

main_p <- p_prev1/p_inc/p_gr

if(FALSE){
  ggsave("figures/fig2.pdf",main_p,width=8,height=6)
  ggsave("figures/fig2.png",main_p,width=8,height=6)
}

betas_comb1 <- betas_all_quants %>% filter(ver == "seir_TRUE_FALSE") %>%
  left_join(gr_quants %>% rename(t_int=t, t=time)) %>% rename(Location = loc)


betas_comb <- betas_all_models_comb %>% filter(ver == "seir_TRUE_FALSE") %>%
  left_join(gr_quants %>% rename(t_int=t, t=time)) %>% rename(Location = loc)

is_growing_res <- betas_comb %>% 
  mutate(is_growing=median_gr>0) %>%
  mutate(is_pos=beta > 0) %>% 
  group_by(Location, t,is_growing) %>% 
  summarize(prob_growing=sum(is_pos)/n()) %>%
  mutate(is_correct = ifelse(is_growing == (prob_growing > 0.5), "Yes", "No"))

ggplot() + 
  geom_abline(slope=1,linetype="dashed") +
  geom_hline(yintercept=0,col=AAAS_palette["grey1"]) +
  geom_violin(data=betas_comb %>% left_join(is_growing_res) %>% rename(`Correct direction` = is_correct),
              aes(x=median_gr, y=beta,group=median_gr,fill=`Correct direction`),
              scale="width",
              width=0.005,trim=TRUE,draw_quantiles = c(00.025,0.5,0.975),alpha=0.5) +
  ylab("Ct growth rate estimate") +
  xlab("SEEIRR posterior median growth rate") +
  ggsci::scale_fill_aaas() +
  coord_cartesian(xlim=c(-0.3,0),ylim=c(-0.6,0.3)) +
  export_theme
