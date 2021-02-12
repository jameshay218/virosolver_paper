coeff <- 30000
p1A <- ggplot(seir_dynamics$seir_outputs)+ 
  geom_line(aes(x=step,y=inc),col="red")+
  geom_line(aes(x=step,y=Rt*coeff),col="forestgreen") +
  scale_y_continuous(expand=c(0,0),limits=c(0,80000),
                     sec.axis=sec_axis(~.*1/coeff, name="Rt")) +
  scale_x_continuous(limits=c(50,200),breaks=seq(0,200,by=25)) +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New infections") +
  labs(tag="A")
p1A


plot_dat <- simulated_viral_loads %>% mutate(week = floor(sampled_time/7)) %>%
  group_by(week) %>%
  mutate(first_day = min(sampled_time,na.rm=TRUE)) %>%
  filter(ct_obs < pars["intercept"]) 
use_days <- plot_dat %>% 
  group_by(first_day) %>%
  tally() %>% 
  filter(n >= 50) %>% 
  pull(first_day)

subset_dat <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  group_by(first_day) %>% sample_n(pmin(50,n()))

p2 <- plot_dat %>% 
  filter(first_day %in% use_days) %>%
  ggplot() + 
  geom_violin(aes(x=first_day, y=ct_obs,group=week),scale="width",
              fill="grey70",alpha=0.75,color=NA) +
  #geom_jitter(aes(x=first_day, y=panther_Ct,group=week),width=2,height=0,
  #            fill="grey70",size=0.25,alpha=0.25) +
  geom_dotplot(data=subset_dat,aes(x=first_day, y=ct_obs,group=first_day),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  geom_smooth(data= . %>% group_by(first_day) %>% summarize(median_ct=median(ct_obs)),
              aes(x=first_day,y=median_ct),col="blue",se=FALSE) +
  scale_y_continuous(trans="reverse",limits=c(45, 5),expand=c(0,0)) +
  scale_x_continuous(limits=c(50,200),breaks=seq(0,200,by=25)) +
  geom_hline(yintercept=40,linetype="dashed") +
  theme_bw() +
  theme_classic()+
  xlab("Day of sample") +
  ylab("Ct value") +
  theme(legend.position="none",
        plot.title=element_blank(),
        panel.grid.minor =element_blank()) +
  labs(tag="B")
p2

p1A/p2

seir_weekly <- seir_dynamics$seir_outputs %>% 
  as_tibble() %>% 
  rename(sampled_time=step) %>%
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
  group_by(week, first_day) %>%
  summarize(median_ct=median(ct_obs),
            skew_ct=moments::skewness(ct_obs),
            n=n()) %>%
  #rename(date=first_day) %>%
  full_join(seir_weekly) %>%
  filter(n >= 50)


pA <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=skew_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=skew_ct),size=2) +
  geom_vline(xintercept=1,linetype="dashed") +
  #scale_y_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  theme_classic() +
  ylab("Skewness of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="C")

pB <- combined_dat %>%
  ggplot() +
  geom_smooth(aes(x=mean_rt,y=median_ct),se=TRUE)+
  geom_point(aes(x=mean_rt,y=median_ct),size=2)+
  geom_vline(xintercept=1,linetype="dashed") +
  theme_classic() +
  #scale_y_continuous(trans="reverse",limits=c(35,28),breaks=seq(28,35,by=1)) +
  scale_y_continuous(trans="reverse") +
  ylab("Median of Ct distribution") +
  xlab("Rt (posterior mean)") +
  labs(tag="D")
pB
pC <- ggplot(combined_dat) +
  geom_point(aes(x=skew_ct,y=median_ct,col=mean_rt),alpha=0.9,size=2) +
  scale_color_gradient2(low="green",mid="blue",high="red",midpoint=1.25,
                        limits=c(0,2.5),
                        guide=guide_colorbar(title="Rt",
                                             barwidth=1,ticks=FALSE,barheight=4))+
  #scale_x_continuous(limits=c(-1.2,0),breaks=seq(-1.2,0,by=0.2)) +
  #scale_y_continuous(trans="reverse",limits=c(35,28),breaks=seq(28,35,by=1)) +
  scale_y_continuous(trans="reverse") +
  theme_classic() +
  xlab("Skewness of Ct distribution") +
  ylab("Median of Ct distribution") +
  theme(legend.position=c(0.2,0.8),
        legend.text=element_text(size=6),
        legend.title=element_text(size=6)) +
  labs(tag="E")


main_p <- p1A/p2/(pA|pB|pC)
main_p
