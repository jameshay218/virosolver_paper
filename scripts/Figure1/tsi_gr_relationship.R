CtDF <- Ct_res[1+vertvals,]
CtDFf <- NULL
for (i in 1:(dim(CtDF)[1])) { ## Assigns [sample_n] points to VLs and test days by coarsening cumulative distribution
  probs <- CtDF[i,]*sample_n
  cum.probs <- cumsum(probs)
  cum.probs.r <- floor(cum.probs+.5)
  NumPerCt <- diff(c(0,cum.probs.r))
  Cts <- rep(x, times=NumPerCt)
  CtDFf <- rbind(CtDFf, data.frame(TestDay=rep(vertvals[i], length(Cts)), Ct=Cts))
}

DF.h.2.age <- data.frame(Days=0:199, Median=age_mean)
skew <- function(values,weights) {
  weights_std <- weights/sum(weights, na.rm=TRUE)
  xbar <- sum(values*weights_std, na.rm=TRUE)
  xi_xbar <- values - xbar
  return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
}


CtDFf$Val <- -1*CtDFf$Ct/100+.3
x1 <- 1:35
age_skews <- apply(age_res_std, 1, function(row) skew(values=x1, weights=row))
ageSkewDF <- data.frame(Days=0:199, Ageskew=age_skews)
ageSkewDF <- merge(ageSkewDF, GR, by="Days", all.x=FALSE, all.y=FALSE)
ageSkewDF$SkewVal <- ageSkewDF$Ageskew*2/6+1/6
ageSkewDF$Label <- paste0("Skew: ",format(round(ageSkewDF$Ageskew, 2), nsmall=2))

## Age skewness
ageSkewDF2 <- tibble(merge(ageSkewDF[,c("Days","Ageskew","DailyGR")], DF.h.2.age[,c("Days","Median")], 
                           by="Days", all.x=TRUE, all.y=FALSE))
ageSkewDF2$ScaleMed <- ageSkewDF2$Median#*0.16-5.624
ageSkewDF2 <- pivot_longer(ageSkewDF2 %>% dplyr::select(-Median), cols=c(Ageskew,ScaleMed), names_to="Type", values_to="Value")
ageSkewDF2$Type <- ifelse(ageSkewDF2$Type=="Ageskew","Skew","Median")

new_skew_df <- ageSkewDF2 %>% filter(Days > lastday) %>% mutate(var="Time since infection") 

p_tsi_relationship <- ggplot() + 
  geom_vline(xintercept = 0,col="grey40",linetype="dashed",size=0.5) +
  geom_path(data=new_skew_df %>% filter(var == "Time since infection" &
                                          Type == "Median"),
            aes(x=DailyGR, y=Value), size=1.25, col=cbbPalette[3]) + 
  labs(x="Daily growth rate", y="Median days since infection") + 
  theme_light() +
  scale_x_continuous(name="Growth rate", limits=c(.11,-.11), breaks=seq(-.1,.1,by=.05), 
                     expand=expansion(add=c(0,0)),
                     trans="reverse") +
  scale_y_continuous(limits=c(7,20),breaks=seq(8,20,by=2)) +
  scale_linetype_manual(name="Distribution Measure", values=c("solid","dotdash")) +
  scale_color_manual(name="Variable", values=cbbPalette[3:4])  +
  main_theme +
  theme(legend.position = "bottom",legend.direction="horizontal",
        #plot.margin = margin(t=5,r=20,b=5,l=10,unit="pt"),
        panel.grid.major=element_line(size=0.1,colour="grey40")) +
  geom_point(data=ageSkewDF2 %>% filter(Days %in% vertvals &
                                          Type == "Median"), 
             shape=19,
             aes(x=DailyGR, y=Value),
             col="black", size=2) +
  
  #geom_text_repel(data=ageSkewDF2 %>% filter(Days %in% vertvals &
  geom_text(data=ageSkewDF2 %>% filter(Days %in% vertvals &Type == "Median", Days %in% c(50,75,100)), 
             aes(x=DailyGR-0.02, y=Value-0.5, label=paste0("Test day ", Days)),
             size=3) +
  geom_text(data=ageSkewDF2 %>% filter(Days %in% vertvals &Type == "Median", Days %in% c(125,150)), 
            aes(x=DailyGR+0.02, y=Value+0.5, label=paste0("Test day ", Days)),
            size=3) +
  #geom_segment(aes(x=0.08,xend=-0.06,y=10,yend=18),lineend="butt",linejoin="mitre",col="grey30",
  #             arrow=arrow(length=unit(0.25,"cm"),type="closed"),size=1)+
  #scale_shape_manual(name="Test Day", values=21:25) +
  guides(color=guide_legend(title.position = "top", title.hjust = 0.5),
         shape=guide_legend(title.position = "top", title.hjust = 0.5)) +
  labs(tag="B")
  
