
dat1 <- ct_densities %>% filter(t %in% seq(25,175,by=25))
dat1 <- dat1 %>% group_by(t) %>% mutate(cumu_dens=cumsum(density))
dat1 <- dat1 %>% mutate(ct_group=as.character(floor(ct/5)))

ct_group_key <- c("0"="(0-5]","1"="(5-10]","2"="(10-15]","3"="(15-20]","4"="(20-25]","5"="(25-30]","6"="(30-35]","7"="(35-40]")
dat1$ct_group <- ct_group_key[dat1$ct_group]
dat1$ct_group <- factor(dat1$ct_group,levels=ct_group_key)
dat2 <- dat1 %>% group_by(ct_group, t) %>% summarize(density=sum(density))

samps <- NULL
for(i in 1:1000){
  samps[[i]] <- dat2 %>% mutate(i=i,sampled=rbinom(n(),60,prob=density))
}
samps <- do.call("bind_rows",samps)
samps_quants <- samps %>% group_by(ct_group,t) %>%
  summarize(lower=quantile(sampled,0.025),
            median=median(sampled),
            upper=quantile(sampled,0.975))


p_normal <- ggplot(dat1) + 
  geom_bar(aes(x=ct,y=density),stat="identity",fill="grey50",col="black",size=0.25) + 
  ylab("Probability density") + xlab("")+
  facet_wrap(~t,nrow=1) + theme_bw()
p_normal_group <- ggplot(dat2) + 
  geom_bar(aes(x=ct_group,y=density),stat="identity",fill="grey50",col="black",size=0.25) + 
  ylab("Probability density") +
  facet_wrap(~t,nrow=1) + theme_bw() + xlab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1)) 

p_observed <- ggplot(samps_quants) + 
  geom_pointrange(aes(x=ct_group,y=median,ymin=lower,ymax=upper),col="blue") + 
  ylab("Number observed") +
  ggtitle("Simulated 60 observations") +
  facet_wrap(~t,nrow=1) + theme_bw() + xlab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1)) 

p_normal_cumu <- ggplot(dat1) + 
  geom_line(aes(x=ct,y=cumu_dens),size=1,col="red")+ 
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  scale_x_continuous(breaks=seq(0,40,by=5)) +
  facet_wrap(~t,nrow=1) + 
  ylab("Cumulative density") + xlab("")+
  geom_hline(yintercept=c(0.25,0.5,0.75),linetype="dashed") + theme_bw()

pA1 <- pA + theme(axis.text.x=element_text(size=8)) + labs(tag="")
main_p <- pA1/p_normal_cumu/p_normal/p_normal_group/p_observed
ggsave("~/Documents/local_data/Liverpool/distributions_over_time_normal.png",height=10,width=12,units="in",dpi=300)
