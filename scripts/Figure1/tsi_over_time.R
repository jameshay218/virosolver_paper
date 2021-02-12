# Summaries of age and VL distributions:
age_res_std <- age_res/apply(age_res, 1, sum, na.rm=TRUE)
age_mean <- apply(age_res_std, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))
age_res_std_csum <- t(apply(age_res_std, 1, FUN=function(res) cumsum(res)))
age_median <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.5]))
age_lower <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.25]))
age_upper <- apply(age_res_std_csum, 1, FUN=function(res) min((1:lastday)[res >= 0.75]))
plot(x=0:199, y=age_mean, type="l", lwd=2, col="blue")
lines(x=0:199, y=age_median, lwd=2, col="green", lty=2)

NumShifts <- length(vertvals)
AgeDF <- age_res[1+vertvals,]
rownames(AgeDF) <- as.character(vertvals)
AgeDF2 <- NULL
sample_n <- 100000
for (i in 1:(dim(AgeDF)[1])) { ## Assigns [sample_n] points to test days & times of infection by coarsening cumulative distributions
  probs <- AgeDF[i,]*sample_n
  cum.probs <- cumsum(probs)
  cum.probs.r <- floor(cum.probs+.5)
  NumPerAge <- diff(c(0,cum.probs.r))
  ages <- rep(1:35, times=NumPerAge)
  cum.ages <- cumsum(ages)
  AgeDF2 <- rbind(AgeDF2, data.frame(TestDay=rep(vertvals[i], length(ages)), TSI=ages))
}
ADFmeans <- apply(age_res, 1, function(x) sum(x*(1:lastday), na.rm=TRUE)/sum(x, na.rm=TRUE))
ADFmeans <- data.frame(TestDay=0:199, MeanAge=ADFmeans)
ADFmedians <- data.frame(TestDay=0:199, MedianAge=age_median, Var="Median")
ADFlower <- data.frame(TestDay=0:199, MedianAge=age_lower, Var="Lower")
ADFupper <- data.frame(TestDay=0:199, MedianAge=age_upper, Var="Upper")

ADFcomb <- data.frame(TestDay=0:199, Lower=age_lower,Upper=age_upper)
ADFcomb1 <- bind_rows(ADFlower,ADFupper)

AgeDFd <- age_res_std_csum[vertvals+1,]
AgeDFd2 <- data.frame(TestDay=rep(vertvals,3),"measure"=rep(c("Median","P25","P75"),each=length(vertvals)),
                      "TSI"=c(apply(AgeDFd,1,function(x) min((1:lastday)[x>=0.5], na.rm=TRUE)),
                              apply(AgeDFd,1,function(x) min((1:lastday)[x>=.25], na.rm=TRUE)),
                              apply(AgeDFd,1,function(x) min((1:lastday)[x>=.75], na.rm=TRUE))))

age_distns <- age_res_std[vertvals+1,]
AgeDFe <- NULL
for (i in 1:NumShifts) {
  AgeDFe <- rbind(AgeDFe, data.frame(Pos=NumShifts-i, TestDay=vertvals[i], TSI=1:35, 
                                     Distn=age_distns[i,]))
}

AgeDFe2 <- AgeDFd2[,c("TestDay","measure","TSI")]
AgeDFe2$TSI <- AgeDFe2$TSI
AgeDFe2 <- merge(AgeDFe2, AgeDFe, by=c("TestDay","TSI"), all.x=TRUE)

p_tsi_time <- ggplot(data=AgeDF2, aes(x=TestDay,y=TSI)) +  
  #geom_ribbon(data=ADFcomb, aes(x=TestDay, ymin=Lower,ymax=Upper), 
  #          linetype=1,fill=cbbPalette[2],inherit.aes = FALSE,alpha=0.05,col="red") +
  #geom_boxplot(aes(group=TestDay),fill="red",alpha=0.1) +
  #geom_dotplot(aes(group=TestDay), binaxis="y",
  #             binwidth=1,stackdir="center",binpositions="all",dotsize=0.21) +
  
  geom_violin(aes(group=TestDay),fill="#008B45FF",col="black",alpha=0.1,
              draw_quantiles=c(0.025,0.5,0.975), 
              trim=TRUE,scale="width",width=10) + 
  geom_jitter(data=AgeDF2%>% sample_frac(0.02),height=0, width=2,size=0.5,col="black") + 
  
  geom_line(data=ADFmedians, aes(x=TestDay, y=MedianAge,group=Var), 
            linetype=1, inherit.aes=FALSE, size=1.25,
            col=cbbPalette[3]) +
  geom_line(data=ADFcomb1, aes(x=TestDay, y=MedianAge,group=Var), 
            linetype=1, inherit.aes=FALSE, size=1.25,
            col=cbbPalette[3]) +
  main_theme +
  scale_y_continuous(breaks=seq(0,lastday+2,by=5),expand=c(0,0)) +
  coord_cartesian(#ylim=c(-5,lastday+2)
                  ylim=c(-8,lastday+2)
                  ) +
  labs(y="Days since infection", x="Days since start") + 
  theme(legend.title=element_blank(), 
        legend.position=c(0.8,0.8),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.length.x = unit(0.2,"cm"),
        #panel.grid.minor.x=element_blank()
        panel.grid.major=element_line(size=0.1,colour="grey40")) + 
  scale_x_continuous(limits=c(25,175), breaks=vertvals) +
  labs(tag="E")

if(TRUE){
for (i in 0:(NumShifts-1)) {
  ph <- ggplot() + 
    geom_bar(data=AgeDFe[AgeDFe$Pos==i,], 
             aes(x=TSI, y=Distn), 
             stat="identity",
             col="grey70",
             fill="grey70") +
    #geom_histogram(data=AgeDFe[AgeDFe$Pos==i,], 
    #               aes(x=TSI, weight=Distn), breaks=seq(0.5,35.5,by=1),
    #               fill="grey70") +
    geom_point(data=AgeDFe2[AgeDFe2$Pos==i,], aes(x=TSI, y=Distn, shape=measure), inherit.aes=FALSE,
               color=cbbPalette[3], size=2) +
    scale_shape_manual(name="Time Since Infection", values=c("Median"=17, "P25"=18,"P75"=18), 
                       labels=c("Median"="Median", "P25"="25th Percentile", "P75"="75th Percentile")) +
    theme_minimal() + 
    scale_x_continuous(breaks=seq(3,33,by=5), minor_breaks=NULL, expand=expansion(add=c(0,0))) +
    scale_y_continuous(expand=c(0,0)) + 
    coord_cartesian(clip = 'off')+
    theme(
      legend.position="none",
      axis.title=element_blank(), 
      axis.text=element_blank(),
      axis.ticks=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  xval <- min(AgeDFe$TestDay[AgeDFe$Pos==i])
  p_tsi_time <- p_tsi_time + annotation_custom(grob=ggplotGrob(ph), xmin=xval-12, xmax=xval+12, 
                                                 ymin=-8.5, ymax=0)
}
}
