rt_dat_top <- rt_dat %>%
  mutate(bottom=pmax(1,bottom),
         top=pmax(1,top))
rt_dat_bot <- rt_dat %>%
  mutate(bottom=pmin(1,bottom),
         top=pmin(1,top))

omg <- rt_dat %>% mutate(is_grow=ifelse(mean>1,"Growing","Declining"))
index <- 1
omg$index <- index
for(i in 2:nrow(omg)){
  if(omg$is_grow[i] != omg$is_grow[i-1]) {
    index <- index + 1
  }
  omg$index[i] <- index
}
coeff <- 2000
p1 <- nyt_dat %>% 
  ## Highlight anomalous cases
  mutate(new_cases1=ifelse(new_cases > 3500, lag(new_cases,1), new_cases)) %>%
  ## Get rolling mean
  mutate(roll_mean=rollmean(new_cases1, 7,fill=NA,align="right")) %>%
  ggplot()+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.5,fill=AAAS_palette["grey1"]) +
  #geom_ribbon(data=dat_infections,aes(x=date,ymin=bottom,ymax=top),
  #            fill=AAAS_palette["grey1"],alpha=0.5) +
  #geom_line(data=dat_infections,aes(x=date,y=mean),col=AAAS_palette["grey1"]) +
  geom_hline(yintercept=coeff,linetype="dashed",col=AAAS_palette["grey1"]) +
  geom_ribbon(data=rt_dat_top,aes(x=date,ymin=bottom*coeff,ymax=top*coeff),
              fill=AAAS_palette["red1"],alpha=0.25) +
  geom_ribbon(data=rt_dat_bot,aes(x=date,ymin=bottom*coeff,ymax=top*coeff),
              fill=AAAS_palette["green1"],alpha=0.25) +
  geom_line(data=omg,
            aes(x=date,y=mean*coeff,col=is_grow,group=index)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,5000),
                     sec.axis=sec_axis(~.*1/coeff, name="Rt")) +
  #scale_fill_npg() +
  scale_color_manual(values=c("Growing"=as.character(AAAS_palette["red1"]),"Declining"=as.character(AAAS_palette["green1"])))+
  scale_x_date(limits=as.Date(c("2020-04-01", "2020-09-01"), "%Y-%m-%d"), breaks="7 days",
               expand=c(0,0)) +
  export_theme+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New infections") +
  labs(tag="A")

p1

trajs1_quants_joined <- left_join(trajs1_quants, date_key)

gr_dat <- estimates$estimates$summarised %>% filter(variable == "growth_rate") %>% 
  dplyr::select(date, mean,top,bottom) %>%
  rename(gr_mean=mean,
         gr_lower=bottom,
         gr_upper=top)


gr_dat_top <- gr_dat %>%
  mutate(gr_lower=pmax(0,gr_lower),
         gr_upper=pmax(0,gr_upper))
gr_dat_bot <- gr_dat %>%
  mutate(gr_lower=pmin(0,gr_lower),
         gr_upper=pmin(0,gr_upper))

omg1 <- gr_dat %>% mutate(is_grow=ifelse(gr_mean>0,"R(t), growing","R(t), declining"))
index <- 1
omg1$index <- index
for(i in 2:nrow(omg1)){
  if(omg1$is_grow[i] != omg1$is_grow[i-1]) {
    index <- index + 1
  }
  omg1$index[i] <- index
}


trajs_dat_top <- trajs1_quants_joined %>%
  mutate(lower=pmax(0,lower),
         upper=pmax(0,upper))
trajs_dat_bot <- trajs1_quants_joined %>%
  mutate(lower=pmin(0,lower),
         upper=pmin(0,upper))

omg2 <- trajs1_quants_joined %>% mutate(is_grow=ifelse(median>0,"Ct estimate","Ct estimate"))
index <- 1
omg2$index <- index
for(i in 2:nrow(omg2)){
  if(omg2$is_grow[i] != omg2$is_grow[i-1]) {
    index <- index + 1
  }
  omg2$index[i] <- index
}




p_grs <- ggplot() + 
  geom_ribbon(data=gr_dat_top,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill=AAAS_palette["red1"]) +
  geom_ribbon(data=gr_dat_bot,aes(x=date,ymin=gr_lower,ymax=gr_upper),alpha=0.25,fill=AAAS_palette["green1"]) +
  geom_line(data=omg1,aes(x=date,y=gr_mean,col=is_grow,group=index)) +
  geom_ribbon(data=trajs_dat_top,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) +
  geom_ribbon(data=trajs_dat_bot,aes(x=date,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) +
  geom_line(data=omg2,aes(x=date,y=median,col=is_grow,group=index)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  scale_x_date(expand=c(0,0),limits=as.Date(c("2020-03-09", "2020-09-01"), "%Y-%m-%d"), 
               breaks=as.Date(c("2020-04-01","2020-06-01","2020-08-01"))) +
  geom_hline(yintercept=0,linetype="dashed",col=AAAS_palette["grey1"]) +
  scale_color_manual(values=c("R(t), growing"=as.character(AAAS_palette["red1"]),
                              "R(t), declining"=as.character(AAAS_palette["green1"]),
                              "Ct estimate"=as.character(AAAS_palette["blue1"])))+
  export_theme +
  theme(legend.position=c(0.5,0.9),
        legend.title = element_blank(),
        axis.text.x=element_text(size=6),
        legend.text=element_text(size=6)) +
  xlab("Date") +
  ylab("Growth rate") +
  labs(tag="F")

