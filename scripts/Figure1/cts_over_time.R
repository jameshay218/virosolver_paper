ct_densities <- dat %>% filter(ct < 40) %>%
  group_by(t) %>%
  mutate(density=density/sum(density))
Ct_median <- ct_densities %>% group_by(t) %>% 
  mutate(cumu_dens=cumsum(density)) %>%
  filter(cumu_dens >= 0.5) %>%
  filter(ct == min(ct))%>%
  mutate(var="Median")
Ct_lower <- ct_densities %>% group_by(t) %>% 
  mutate(cumu_dens=cumsum(density)) %>%
  filter(cumu_dens >= 0.25) %>%
  filter(ct == min(ct))%>%
  mutate(var="Lower")
Ct_upper <- ct_densities %>% group_by(t) %>% 
  mutate(cumu_dens=cumsum(density)) %>%
  filter(cumu_dens >= 0.75) %>%
  filter(ct == min(ct)) %>%
  mutate(var="Upper")

inf_times <- simulate_infection_times(100000, incidence)
add_noise <- function(n, vl, sd){rgumbel(n,vl,sd)}
vls <- simulate_viral_loads(inf_times, 0:200,pars,add_noise = add_noise)

p_ct_time <- ggplot() +
  geom_violin(data=vls%>% filter(t %in% vertvals,obs < 40, obs > 0),aes(x=t,y=obs,group=t),
              fill="#631879FF",alpha=0.1,col="black",draw_quantiles=c(0.025,0.5,0.975), trim=TRUE,scale="width",width=10) +
  geom_jitter(data=vls %>% filter(t %in% vertvals,obs < 40, obs > 0) %>% sample_n(500), aes(x=t,y=obs),width=2,
              size=0.5,col="black") +
  geom_line(data=Ct_lower,aes(x=t,y=ct),col="#631879FF",size=1.25) +
  geom_line(data=Ct_upper,aes(x=t,y=ct),col="#631879FF",size=1.25) +
  geom_line(data=Ct_median,aes(x=t,y=ct),col="#631879FF",size=1.25) +
  
  scale_y_continuous(trans="reverse",limits=c(40,8),expand=c(0,0),breaks=seq(10,40,by=5)) +
  coord_cartesian(ylim=c(49,8)) +
  main_theme +
  labs(y="Ct value", x="Days since start of outbreak") + 
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        #axis.text.x = element_blank(),
        #axis.line.x=element_blank(),
        #axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.length.x = unit(0.2,"cm"),
        #panel.grid.minor.x=element_blank()
        panel.grid.major=element_line(size=0.1,colour="grey40")) + 
  scale_x_continuous(limits=c(25,175), breaks=vertvals) +
  labs(tag="G")
if(TRUE){
  ct_quants <- bind_rows(Ct_lower,Ct_median,Ct_upper)
  for (xval in vertvals) {
    max_y <- ct_quants %>% filter(t == xval,density==max(density)) %>% pull(density)
    
    ph <- ggplot() + 
      geom_histogram(data=ct_densities %>% filter(t == xval), aes(x=ct,weight=density), 
                     breaks=c(-0.5,seq(4.5,40.5,by=1)),col="grey70",fill="grey70",
                     inherit.aes=FALSE) + 
      theme_minimal() + 
      geom_point(data=ct_quants %>% filter(t == xval), aes(x=ct, y=density, shape=var), inherit.aes=FALSE,
                 size=2,color="#631879FF") +
      scale_shape_manual(name="Time Since Infection", values=c("Median"=17, "Lower"=18,"Upper"=18), 
                         labels=c("Median"="Median", "Lower"="25th Percentile", "Upper"="75th Percentile")) +
      scale_x_continuous(breaks=seq(2,38,by=4), minor_breaks=NULL, expand = c(0,0), trans="reverse")  +
      scale_y_continuous(expand=c(0,0)) +
      coord_cartesian(clip = 'off')+
      theme(axis.title=element_blank(), 
            legend.position="none",
            axis.text=element_blank(),
            axis.ticks.length=unit(0, "pt"),
            panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    
    p_ct_time <- p_ct_time + annotation_custom(grob=ggplotGrob(ph), xmin=xval-10, xmax=xval+15, 
                                 ymin=-49.5, ymax=-41)
  }
}
