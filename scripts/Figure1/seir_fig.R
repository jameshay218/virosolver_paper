#### Setting Testing Days to Highlight in Figures: ####
vertvals <- c(50,75,100,125,150)
GRs <- format(round(GR$GR[GR$Days %in% vertvals],3), nsmall=3)


#### Panel A: ####
DF.a <- data.frame(Days=c(0:199,GR$Days,GR$Days), Value=c(incidence,GR$GR/10+.01,GR$DailyGR/10+.01), 
                   Type=c(rep("Incidence",length(incidence)),rep("35-day growth rate",length(GR$GR)),rep("Daily growth rate", length(GR$DailyGR))))
DF.a$Typef <- factor(DF.a$Type, levels=c("Incidence", "Daily growth rate", "35-day growth rate"))

xlength <- 10
Tangents <- data.frame(Days=vertvals, dailyGR=GR$DailyGR[GR$Days %in% vertvals], incidence=DF.a$Value[DF.a$Days %in% vertvals & DF.a$Type=="Incidence"])
Tangents$slope1 <- (exp(Tangents$dailyGR)-1)
Tangents$slope <- Tangents$slope1*Tangents$incidence
Tangents$x <- Tangents$Days-xlength/2
Tangents$xend <- Tangents$Days+xlength/2
Tangents$y <- Tangents$incidence-Tangents$slope*xlength/2
Tangents$yend <- Tangents$incidence+Tangents$slope*xlength/2
Tangents$y2 <- Tangents$dailyGR/10+.01-Tangents$slope*xlength/2
Tangents$yend2 <- Tangents$dailyGR/10+.01+Tangents$slope*xlength/2


p_seir <- ggplot(data=DF.a %>% filter(Type == "Daily growth rate"), 
                    aes(x=Days,y=Value, group=Typef)) + 
  geom_bar(data=DF.a %>% filter(Type == "Incidence"),aes(fill="Incidence"),stat="identity",col="black",size=0.25) +
  geom_line(size=1,aes(color="Daily growth rate")) + 
  theme_light() +
  scale_color_manual(values=c("Daily growth rate" = "#3B4992FF")) +
  scale_fill_manual(values=c("Incidence"="#808180FF")) +
  #scale_color_aaas() +
  scale_y_continuous(expand=c(0,0),
                     name="Per capita incidence", limits=c(0,0.021),
                     sec.axis=sec_axis(trans=~.*10-.1, name="Growth rate")) + 
  geom_hline(yintercept = 0.01,col="black",linetype="dashed",size=0.5) +
  main_theme +
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.length.x = unit(0.2,"cm"),
        #panel.grid.minor.x=element_blank()
        panel.grid.major=element_line(size=0.1,colour="grey40")
        ) + 
  #geom_vline(xintercept=vertvals, col="grey40",linetype="dotted",size=0.5) + 
  scale_x_continuous(limits=c(25,175), breaks=vertvals) +
  labs(tag="A")
