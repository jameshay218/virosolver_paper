VLcurve <- data.frame(TSI=seq(0,35,by=0.1))

VLcurve$Mean <- viral_load_func(pars, VLcurve$TSI, convert_vl=FALSE)
VLcurve$IL <- qnorm(p=0.025, mean=VLcurve$Mean, sd=pars["obs_sd"])
VLcurve$UL <- qnorm(p=0.975, mean=VLcurve$Mean, sd=pars["obs_sd"])
# VLcurve$IL <- extraDistr::qgumbel(p=0.025, mean=VLcurve$Mean, sigma=pars["obs_sd"])
# VLcurve$UL <- extraDistr:qgumbel(p=0.975, mean=VLcurve$Mean, sigma=pars["obs_sd"])
VLcurve.long <- as.data.frame(pivot_longer(VLcurve, cols=c("Mean","IL","UL"), names_to="measure", values_to="Ct"))

VLcurve_means <- VLcurve.long %>% filter(measure == "Mean")

VL_ribbon <- VLcurve.long %>% filter(measure != "Mean") %>%
  pivot_wider(names_from=measure,values_from=Ct)

tsi_points <- ADFmedians %>% 
  select(-Var) %>%
  filter(TestDay %in% vertvals) %>%
  rename(TSI=MedianAge) %>%
  left_join(VLcurve.long %>% filter(measure=="Mean")) %>%
  drop_na()


p_vl <- ggplot() +
  geom_ribbon(data=VL_ribbon, aes(xmin=IL,xmax=UL, y=TSI), fill=cbbPalette[4], alpha=0.1) +
  geom_line(data=VLcurve_means,aes(x=Ct,y=TSI),orientation="y", color=cbbPalette[4],size=1.25) +
  geom_point(data=tsi_points,aes(x=Ct, y=TSI),shape=19, col="black",size=2) +
  #geom_text_repel(data=tsi_points, 
  geom_text(data=tsi_points %>% filter(TestDay %in% c(50,75,100)), aes(x=Ct+5, y=TSI-1, label=paste0("Test day ", TestDay)),size=3) +
  geom_text(data=tsi_points %>% filter(TestDay %in% c(125)),aes(x=Ct-5, y=TSI, label=paste0("Test day ", TestDay)),size=3) +
  geom_text(data=tsi_points %>% filter(TestDay %in% c(150)),aes(x=Ct-5, y=TSI+1, label=paste0("Test day ", TestDay)),size=3) +
  #geom_segment(data=VLcurve_means,aes(x=20,xend=35,y=15,yend=27),lineend="butt",linejoin="mitre",col="grey30",
  #             arrow=arrow(length=unit(0.25,"cm"),type="closed"),size=1)+
  ylab("Days since infection") +
  xlab("Mean Ct value") +
  main_theme + 
  theme(legend.position = "bottom",legend.direction="horizontal",
        #plot.margin = margin(t=5,r=20,b=5,l=10,unit="pt"),
        panel.grid.major=element_line(size=0.1,colour="grey40")) +
  
  scale_y_continuous(breaks=seq(0,lastday+2,by=5),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(10,40,by=5), trans="reverse") +
  coord_cartesian(#ylim=c(-5,lastday+2),
                  ylim=c(0,lastday+2),
                  xlim=c(42,10)) +
  # theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  #theme(axis.line.y=element_blank(),
  #      axis.ticks.y=element_blank(),
  #      axis.text.y=element_blank(),
  #      axis.title.y=element_blank()) +
  labs(tag="D")

