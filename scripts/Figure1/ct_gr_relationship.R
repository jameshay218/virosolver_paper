skew <- function(values,weights) {
  weights_std <- weights/sum(weights, na.rm=TRUE)
  xbar <- sum(values*weights_std, na.rm=TRUE)
  xi_xbar <- values - xbar
  return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
}

ct_skew <- ct_densities %>% filter(t >= 35) %>% group_by(t) %>% 
  mutate(mean_ct=sum(ct*density),top_part=(ct-mean_ct)^3,bot_part=(ct-mean_ct)^2) %>% 
  group_by(t) %>% summarize(top=sum(top_part*density),bot=(sum(bot_part*density)^1.5),skew=top/bot) 

ggplot(ct_skew) + geom_line(aes(x=t,y=skew))
ph_dat <- ct_skew %>% left_join(GR %>% rename(t=Days)) %>% filter(t > 35)
p_ct <- ggplot() + 
  geom_vline(xintercept = 0,col="grey40",linetype="dashed",size=0.5) +
  geom_path(data=ph_dat,
            aes(x=DailyGR, y=skew), size=1, col=cbbPalette[4]) + 
  labs(x="Daily growth rate", y="Ct distribution skew") + 
  theme_light() +
  scale_x_continuous(name="Growth rate", limits=c(.11,-.11), breaks=seq(-.1,.1,by=.05), 
                     expand=expansion(add=c(0,0)),
                     trans="reverse") +
  #scale_y_continuous(limits=c(10,18),breaks=seq(10,18,by=2)) +
  scale_linetype_manual(name="Distribution Measure", values=c("solid","dotdash")) +
  main_theme +
  theme(legend.position = "bottom",legend.direction="horizontal",
        #plot.margin = margin(t=5,r=20,b=5,l=10,unit="pt"),
        panel.grid.major=element_line(size=0.1,colour="grey40")) +
  geom_point(data=ph_dat[ph_dat$t %in% vertvals,], 
             shape=18,
             aes(x=DailyGR, y=skew),
             col="black", size=3) +
  geom_text_repel(data=ph_dat[ph_dat$t %in% vertvals,], 
                  aes(x=DailyGR, y=skew, label=paste0("Test day ", t)),
                  size=3) +
  scale_shape_manual(name="Test Day", values=21:25) +
  guides(color=guide_legend(title.position = "top", title.hjust = 0.5),
         shape=guide_legend(title.position = "top", title.hjust = 0.5)) +
  labs(tag="F")


