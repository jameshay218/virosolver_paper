#### Panel D Alt Heatmap: 
###### Function to generate random viral load trajectories:
rx <- function(pars, obs_model="normal", additional_detect=FALSE, a.range=1:35, corr=.8) {
  ct_mean <- viral_load_func(pars, a.range)
  if (obs_model=="normal") {
    Sigma <- matrix(pars["obs_sd"]^2, nrow=length(a.range), ncol=length(a.range))
    for (i in 1:length(a.range)) {
      for (j in 1:length(a.range)) {
        Sigma[i,j] <- Sigma[i,j]*corr^(abs(i-j))
      }
    }
    xs <- rmvnorm(1, mean=ct_mean, sigma=Sigma)[1,]
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  } else {
    ind <- extraDistr::rgumbel(1, mu=0, sigma=pars["obs_sd"])
    xs <- ct_mean + ind
    xs[xs >= pars["intercept"]] <- NA
    xs[xs < 1] <- 1
  }
  if (additional_detect) {
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    t_switch_int <- ceiling(t_switch)
    t_switch_inc <- t_switch_int-t_switch
    ts <- a.range[a.range >= t_switch_int]
    addl_det <- c(rbinom(1, 1, (1-pars["prob_detect"]*pars["t_unit"])^t_switch_inc),
                  rbinom(length(ts)-1, 1, (1-pars["prob_detect"]*pars["t_unit"])))
    last.day <- min(c(max(a.range)+1,ts[addl_det==0]))-1
    xs[a.range > last.day] <- NA
  }
  xs
}

set.seed(3712)
times <- sort(sample(x=0:199, size=500, replace=TRUE, prob=incidence), decreasing=TRUE)
DF.heat <- data.frame()
for (n in 1:500) {
  Person <- n
  Days <- times[n] + 1:lastday
  Cts <- rx(pars, obs_model="gumbel", additional_detect=TRUE, a.range=1:lastday, corr=0.9)
  DF.heat <- rbind(DF.heat, data.frame(Person, Days, Cts))
}

p_traj <- ggplot(data=DF.heat) + 
  geom_tile(aes(x=Days,y=as.factor(Person),fill=Cts)) + 
  # geom_rect(aes(xmin=101,xmax=181,ymin=224,ymax=251),fill=NA,col="grey10") +
  #scale_fill_viridis_c(limits=c(0,10),breaks=seq(0,10,by=2)) +
  scale_fill_gradientn(name="Ct value",
                       colours = c("white","#00468bff","#ffdd55ff","#ff751a","#ed0000ff","#ed0000ff","#ed0000ff"),
                       breaks=seq(0,40,by=8), labels=seq(0,40,by=8), na.value="transparent",
                       trans="reverse", limits=c(42,-2)) + 
  scale_x_continuous(limits=c(25,175), breaks=seq(50,150,by=25)) +
  main_theme +
  xlab("Days since start of outbreak") +
  theme(legend.position=c(0.9,0.6),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.length.x = unit(0.2,"cm"),
        #panel.grid.minor.x=element_blank()
        panel.grid.major.x=element_line(size=0.1,colour="grey40"),
    
      axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.background=element_blank()) +
  labs(tag="C", y="Individuals")
