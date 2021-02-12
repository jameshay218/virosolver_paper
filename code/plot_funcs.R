
## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))


export_theme <- theme_tufte() +
  theme(
    axis.text.x = element_text(size=7,family="sans"),
    axis.text.y=element_text(size=7,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),

    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),

    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    plot.tag = element_text(family="sans",size=10,face="bold"),

    ## Legends
    legend.title=element_text(size=8,family="sans",face="italic"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans",face="bold"),
    #strip.background=element_rect(fill="#f0f0f0")
    strip.background=element_blank()
    )



AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="#1B1919FF")


theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(0 - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(0, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    #theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
    #                                         data = xtick, size = size, 
    #                                         color = color)
    
    #Add labels to the x-ticks
    #theme.list[[nlist + 4*k-2]] <- annotate("text", 
    #                                        x = ticks_x[k], 
    #                                        y = ygeo - 2.5*epsilon,
    #                                        size = textsize,
    #                                        label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}

plot_nursing_home <- function(chain, obs_dat, plot_pars=c("viral_peak","R0","t0"), loc="NH 1"){
  
  
  ## Plot trajectories
  p_trajectories <- ggplot(posterior_dat) +
    geom_line(aes(x=t,y=40-prob_infection/0.004,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat,aes(x=t,y=40-prob_infection/0.004,col="MAP"),size=0.5) +
    geom_violin(data=obs_dat,aes(x=t, y=ct,group=t),
                width=4,scale="width",fill="grey70",alpha=0.25,col="grey70",
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_jitter(data=obs_dat,aes(x=t, y=ct),size=0.25,height=0,width=1.25) +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                         breaks=seq(0,0.12,by=0.03))) +
    xlab("Time") +
    ylab("Probability of infection") +
    theme_classic() + +
    scale_color_manual(values=c("Posterior draw"="gray50","MAP"="green","Ground truth"="blue")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position="bottom") +
    facet_wrap(~t, nrow=1)
  
  ggplot() +
    #geom_vline(data=betas_all_models_comb %>% dplyr::select(loc, t) %>% distinct() %>%
    #             rename(location=loc),
    #           aes(xintercept=t),col=AAAS_palette["grey1"],linetype="dashed",size=0.25) +
   
    
    geom_line(data=dat_use %>% group_by(location,mean_week_date) %>% summarize(y=median(ct, na.rm=TRUE)),
              aes(x=mean_week_date,y=y),col=AAAS_palette["grey"],linetype="longdash",size=0.25) +
    
    geom_ribbon(data=quants_inc,
                aes(ymin=40 - lower/0.004,ymax=40 - upper/0.004,x=time),
                alpha=0.1,fill=AAAS_palette["blue1"]) +
    geom_line(data=best_inc_traj, 
              aes(x=time,y=40 - inc/0.004),
              col=AAAS_palette["blue1"]) +
    
    scale_x_date(limits=range(gr_quants$time), breaks="7 days") +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                         breaks=seq(0,0.12,by=0.03))) +
    ylab("Ct value") +
    xlab("Date") +
    export_theme
  
  ## Plot 
  
  ## Plot densities
  p_densities <- chain[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_density(aes(x=value,fill=chain),alpha=0.25) +
    ggsci::scale_fill_aaas() +
    facet_wrap(t~name,scales="free",nrow=1) +
    export_theme
  
}


plot_distribution_fits_bwh <- function(chain, obs_dat,MODEL_FUNC, nsamps=100,pos_only=TRUE,date_key){
  
  best_pars <- get_best_pars(chain)
  best_dat <- MODEL_FUNC(best_pars)
  
  ## Generate posterior draws for Ct distribution prediction
  samps <- sample(unique(chain$sampno),nsamps)
  all_res <- NULL
  for(i in seq_along(samps)){
    samp <- samps[i]
    tmp_pars <- lazymcmc::get_index_pars(chain, samp)
    all_res[[i]] <- MODEL_FUNC(tmp_pars) %>% mutate(sampno=i)
  }
  posterior_dat <- do.call("bind_rows",all_res)
  
  ## Make sure only plotting detectable, and get label
  obs_dat1 <- obs_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    mutate(obs_t=paste0("Sample day: ", t))
  
  ## Get number of observations per time point
  obs_tally <- obs_dat1 %>% group_by(t) %>% tally()
  
  ## Re-scale densities to only detectable Ct distribution densities
  total_density <- posterior_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    group_by(t,sampno) %>%
    summarize(total_dens=sum(density)) %>%
    left_join(obs_tally)
  
  ## Get expected number of observations for each Ct,
  ## also simulate observations
  summary_posterior_dat <- posterior_dat %>%
    filter(ct < best_pars["intercept"]) %>%
    left_join(total_density) %>%
    group_by(t, sampno) %>%
    mutate(density=density/total_dens) %>%
    ungroup() %>%
    mutate(expectation=density*n) %>%
    ungroup() %>%
    mutate(sim_obs=rbinom(n(),n,density)) %>%
    group_by(ct, t)
  
  summary_expectation <- summary_posterior_dat %>%
    group_by(ct, t) %>%
    ## Quantiles on expectations
    summarize(lower_expec=quantile(expectation,0.025),
              median_expec=quantile(expectation,0.5),
              upper_expec=quantile(expectation,0.975))
  
  ## Quantiles on observations
  summary_obs <- summary_posterior_dat %>%
    group_by(ct, t) %>%
    summarize(lower_obs=quantile(sim_obs,0.025),
              median_obs=quantile(sim_obs,0.5),
              upper_obs=quantile(sim_obs,0.975))
  
  summary_obs <- summary_obs %>% left_join(date_key)
  summary_expectation <- summary_expectation %>% left_join(date_key)
  obs_dat1 <- obs_dat1  %>% left_join(date_key)
  
  p1 <- ggplot(obs_dat1) +
    geom_histogram(aes(x=ct,y=..count..),binwidth=1,fill="grey70",col="grey20",boundary=0) +
    geom_ribbon(data=summary_obs,aes(x=ct+0.5,ymin=lower_obs,ymax=upper_obs),fill="blue",alpha=0.25)+
    geom_ribbon(data=summary_expectation,aes(x=ct+0.5,ymin=lower_expec,ymax=upper_expec),fill="blue",alpha=0.5)+
    geom_line(data=summary_expectation,aes(x=ct+0.5,y=median_expec),col="blue") +
    #geom_line(data=best_dat%>%group_by(t) %>%filter(ct < best_pars["intercept"]), aes(x=ct+0.5,y=density),col="green") +
    scale_x_continuous(trans="reverse",expand=c(0,0),limits=c(41,5),breaks=seq(0,40,by=5)) +
    #scale_y_continuous(expand=c(0,0),limits=c(0,10),breaks=seq(0,10,by=2)) +
    scale_y_continuous(breaks=c(0,2,seq(5,20,by=5))) +
    coord_cartesian(xlim=c(0,39)) +
    coord_flip() +
    xlab("Ct value") +
    ylab("Count") +
    facet_wrap(~format(date,"%m-%d"),nrow=1,scales="free_x") +
    export_theme +
    theme(panel.grid.major = element_line(color="grey80",size=0.25),
          panel.grid.minor = element_line(color="grey80",size=0.25),
          panel.spacing = unit(0.5, "lines"),
          axis.text.x=element_text(size=5),
          plot.margin = margin(0,0,0,0, "cm"),
          axis.text.y=element_text(size=5),
          strip.text=element_text(size=6),
          strip.background = element_blank()) +
    labs(tag="A")
  
  ## Get predictions for prob undetectable
  summary_prop_detectable <- posterior_dat %>% filter(ct == best_pars["intercept"]) %>%
    group_by(t) %>%
    mutate(density=1-density) %>%
    summarize(lower=quantile(density,0.025),
              median=quantile(density,0.5),
              upper=quantile(density,0.975))
  summary_prop_detectable <- summary_prop_detectable %>% left_join(date_key)
  p2 <-  ggplot(obs_dat1 %>%
                  group_by(t) %>%
                  mutate(is_detectable=ct < best_pars["intercept"]) %>%
                  summarize(prop_detectable=sum(is_detectable)/n()))
  if(!pos_only){
    p2 <- p2 +
      geom_point(aes(y=prop_detectable,x=date-0.5,col="Data"),size=3,shape=18) +
      geom_point(data=summary_prop_detectable,aes(x=date+0.5,y=median,col="Posterior median & 95% CI"),size=1) +
      geom_errorbar(data=summary_prop_detectable,aes(x=date,ymin=lower,ymax=upper,col="Posterior median & 95% CI"),width=0.5)
  } else{
    p2 <- p2 +
      geom_point(data=summary_prop_detectable,aes(x=date,y=median,col="Posterior median & 95% CI"),size=1) +
      geom_errorbar(data=summary_prop_detectable,aes(x=date,ymin=lower,ymax=upper,col="Posterior median & 95% CI"),width=0.5)
  }
  p2 <- p2 +
    #geom_point(data=best_dat %>%filter(ct==best_pars["intercept"]) %>%mutate(density = 1-density),aes(x=0.5,y=density,col="MAP"),size=1) +
    scale_color_manual(values=c("Data"="grey40",
                                "Posterior median & 95% CI"="blue",
                                "MAP"="green")) +
    guides(color=guide_legend(title=NULL)) +
    #facet_wrap(~date,nrow=2) +
    ylab("Proportion detectable") +
    scale_x_date(breaks="1 month",expand=c(0.01,0.01)) +
    xlab("Sample time") +
    export_theme +
    theme(panel.grid.major = element_line(color="grey80",size=0.25),
          panel.grid.minor = element_line(color="grey80",size=0.25),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          strip.text=element_text(size=7),
          #axis.text.x=element_blank(),
          legend.position="none") +
    labs(tag="B")
  list(p1, p2)
}

## Plot predicted trajectories for Figure 2
plot_predictions_fig2 <- function(use_loc, preds_all, best_prob_dat, obs_dat,n_samp=50,ymax=0.12,yscale=0.004){
  samps <- sample(unique(preds_all$sampno),n_samp)
  p_predictions <- ggplot(preds_all %>% filter(loc == use_loc,sampno %in% samps)) + 
    geom_line(aes(x=t,y=prob_infection,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat%>% filter(loc == use_loc),aes(x=t,y=prob_infection,col="MAP"),size=0.5) +
    scale_y_continuous(breaks=seq(0,ymax,by=0.03),expand=c(0,0),limits=c(0,ymax)) +
    xlab("Date") +
    ylab("Per capita incidence") +
    export_theme + 
    scale_color_manual(values=c("Posterior draw"="#EE0000FF","MAP"="#008B45FF")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position=c(0.1,0.8),
          plot.margin = margin(0,0,0,0, "cm")) +
    facet_wrap(~samp_t, nrow=1)+
    labs(tag="C") + 
    coord_cartesian(xlim=as.Date(c("2020-02-25","2020-05-09"))) +
    theme(legend.position=c(0.8,0.9),panel.grid.major=element_line(size=0.1,color="grey40"))+
          #axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank()) +
    facet_wrap(~samp_t,ncol=1)
  p_predictions
}

## Plot predicted trajectories for Figure 2 individual times
plot_predictions_fig2_indiv <- function(use_loc, preds_all, best_prob_dat, obs_dat,n_samp=50,ymax=0.12,yscale=0.004,timepoint){
  preds_all <- preds_all %>% filter(samp_t == timepoint)
  best_prob_dat <- best_prob_dat %>% filter(samp_t == timepoint)
  obs_dat <- obs_dat %>% filter(samp_t == timepoint)
  
  samps <- sample(unique(preds_all$sampno),n_samp)
  p_predictions <- ggplot(preds_all %>% filter(loc == use_loc,sampno %in% samps)) + 
    geom_line(aes(x=t,y=prob_infection,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat%>% filter(loc == use_loc),aes(x=t,y=prob_infection,col="MAP"),size=1) +
    geom_vline(aes(xintercept=samp_t),linetype="dashed") +
    scale_y_continuous(breaks=seq(0,ymax,by=0.03),expand=c(0,0),limits=c(0,ymax)) +
    xlab("Date") +
    ylab("Per capita incidence") +
    export_theme + 
    scale_color_manual(values=c("Posterior draw"="#EE0000FF","MAP"="#5F559BFF")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position=c(0.1,0.8),
          plot.margin = margin(0,0,0.15,0, "cm")) +
    coord_cartesian(xlim=as.Date(c("2020-03-01","2020-05-09"))) +
    scale_x_date(expand=c(0,0)) +
    theme(legend.position=c(0.8,0.8),panel.grid.major=element_line(size=0.1,color="grey40"))
    #axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank()) 
  p_predictions
}


## Plot predicted trajectories for Figure 2
plot_predictions_fig2_old <- function(use_loc, preds_all, best_prob_dat, obs_dat,n_samp=50,ymax=0.12,yscale=0.004){
  samps <- sample(unique(preds_all$sampno),n_samp)
  p_predictions <- ggplot(preds_all %>% filter(loc == use_loc,sampno %in% samps)) + 
    geom_line(aes(x=t,y=40-prob_infection/yscale,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat%>% filter(loc == use_loc),aes(x=t,y=40-prob_infection/yscale,col="MAP"),size=0.5) +
    #geom_violin(data=obs_dat%>% filter(location == use_loc),aes(x=mean_week_date, y=ct,group=mean_week_date),
    #            width=8,scale="width",fill=AAAS_palette[1],alpha=0.25,col="grey10",
    #            draw_quantiles=c(0.025,0.5,0.975)) +
    #geom_jitter(data=obs_dat%>% filter(location == use_loc),aes(x=mean_week_date, y=ct),size=0.25,height=0,width=1.25,col=AAAS_palette[1]) +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*yscale) - .*yscale,name="Per capita incidence",
                                         breaks=seq(0,ymax,by=0.03))) +
    xlab("Date") +
    ylab("Ct value") +
    export_theme + 
    scale_color_manual(values=c("Posterior draw"="#EE0000FF","MAP"="#008B45FF")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position=c(0.1,0.8),
          plot.margin = margin(0,0,0,0, "cm")) +
    facet_wrap(~samp_t, nrow=1)+
    labs(tag="C") + 
    coord_cartesian(xlim=as.Date(c("2020-02-25","2020-05-09"))) +
    theme(legend.position=c(0.8,0.9),panel.grid.major=element_line(size=0.1,color="grey40"))+
    #axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank()) +
    facet_wrap(~samp_t,ncol=1)
  p_predictions
}

## Plot predicted trajectories for Figure 2
plot_predictions_supp_fig2 <- function(dont_use_loc, preds_all, best_prob_dat, obs_dat,n_samp=50,ymax=0.12,yscale=0.004){
  samps <- sample(unique(preds_all$sampno),n_samp)
  time_key <- preds_all %>% filter(loc != dont_use_loc) %>% dplyr::select(loc,samp_t) %>% distinct() %>% group_by(loc) %>% mutate(timepoint=1:n()) %>% rename(location=loc)
  p_predictions <- ggplot(preds_all %>% filter(loc != dont_use_loc,sampno %in% samps) %>% rename(location=loc) %>% left_join(time_key)) + 
    geom_line(aes(x=t,y=40-prob_infection/yscale,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat%>% filter(loc != dont_use_loc)%>% rename(location=loc)%>% left_join(time_key),aes(x=t,y=40-prob_infection/yscale,col="MAP"),size=0.5) +
    geom_violin(data=obs_dat%>% filter(location != dont_use_loc)%>% left_join(time_key),aes(x=mean_week_date, y=ct,group=mean_week_date),
                width=8,scale="width",fill=AAAS_palette[1],alpha=0.25,col="grey10",
                draw_quantiles=c(0.025,0.5,0.975)) +
    #geom_jitter(data=obs_dat%>% filter(location != dont_use_loc)%>% left_join(time_key),aes(x=mean_week_date, y=ct),size=0.25,height=0,width=1.25,col=AAAS_palette[1]) +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*yscale) - .*yscale,name="Per capita incidence",
                                         breaks=seq(0,ymax,by=0.03))) +
    xlab("Date") +
    ylab("Ct value") +
    export_theme + 
    scale_color_manual(values=c("Posterior draw"="#EE0000FF","MAP"="#008B45FF")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position=c(0.1,0.8),
          plot.margin = margin(0,0,0,0, "cm")) +
    facet_wrap(paste0(location, ", timepoint ", timepoint)~.,ncol=3,dir="v") +
    coord_cartesian(xlim=as.Date(c("2020-02-25","2020-05-09"))) +
    theme(legend.position=c(0.8,0.9),panel.grid.major=element_line(size=0.1,color="grey40"))
          #axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank())
  p_predictions
}

## Plot inset distribution fits for Figure 2
inset_theme <-  theme(panel.grid.major = element_line(color="grey80",size=0.25),
                      #axis.text.x=element_text(size=6),
                      #axis.title.x=element_text(size=6),
                      #axis.text.y=element_text(size=6),
                      #axis.title.y=element_text(size=6),
                      panel.grid.minor = element_blank(),
                      #plot.background=element_rect(fill="white",colour=NA),
                      plot.background = element_blank(),
                      panel.background = element_rect(fill="white",colour=NA),
                      plot.margin = margin(0,0,0.15,0, "cm"))

plot_inset_fig2 <- function(use_loc, timepoint, obs_dat, summary_posterior_dat, summary_prob_detectable=NULL, ymax=10,ymax2=0.6){
  p_inset <- ggplot(obs_dat %>% filter(location == use_loc,samp_t==timepoint)) +
    geom_histogram(aes(x=ct),binwidth=1,fill="grey90",col="black",boundary=0) +
    geom_ribbon(data=summary_posterior_dat%>% filter(loc == use_loc,samp_t==timepoint),aes(x=ct+0.5,ymin=lower_count,ymax=upper_count),fill="blue",alpha=0.15)+
    geom_ribbon(data=summary_posterior_dat %>% filter(loc == use_loc,samp_t==timepoint),aes(x=ct+0.5,ymin=lower,ymax=upper),fill="blue",alpha=0.6)+
    geom_line(data=summary_posterior_dat %>% filter(loc == use_loc,samp_t==timepoint),
              aes(x=ct+0.5,y=median)) +
    scale_x_continuous(trans="reverse",expand=c(0,0),limits=c(39,5),breaks=seq(0,40,by=10)) +
    coord_cartesian(xlim=c(39,5),ylim=c(0,ymax)) +
    scale_y_continuous(expand=c(0,0)) +
    #scale_y_continuous(expand=c(0,0),breaks=seq(0,1,by=0.5),limits=c(0,1.01)) +
    export_theme +
    inset_theme +
    xlab("Ct value") +
    ylab("Count") 
  
  if(!is.null(summary_prob_detectable)){
    p_inset2 <-  ggplot(summary_prob_detect  %>% filter(loc == use_loc, samp_t==timepoint)) + 
      geom_errorbar(aes(x=1.2,ymin=lower,ymax=upper),col="blue",width=0.1) + 
      geom_point(aes(x=1.2,y=med),col="blue",size=1) + 
      geom_point(aes(x=0.8,y=dat_detectable),shape=8,col="grey40",size=2) + 
      scale_x_continuous(limits=c(0.6,1.4),breaks=c(0.8,1.2),labels=c("Observed","Predicted")) + 
      scale_y_continuous(expand=c(0,0),limits=c(0,ymax2))  +
      export_theme +
      inset_theme +
      theme(axis.text.x=element_blank())+
      #      axis.title.y=element_text(size=5),
      #      axis.title.x=element_blank(),
      #      plot.margin = margin(0,0,0,0, "cm")) +
      xlab("") +
      ylab("PCR prev") 
    return(list(p_inset,p_inset2))
  } else {
   return(p_inset)
  }
}

plot_inset_fig2_prev_combined <- function(use_loc, summary_prob_detectable, ymax=10,ymax2=0.6){
  use_dat <- summary_prob_detectable  %>% filter(loc == use_loc)
  use_dat$t_plot <- 1:nrow(use_dat)
  xlabels <- unique(use_dat$samp_t)
  
  p_inset2 <-  ggplot(use_dat) + 
    geom_errorbar(aes(x=t_plot+0.1,ymin=lower,ymax=upper),col="blue",width=0.1) + 
    geom_point(aes(x=t_plot+0.1,y=med,shape="Predicted"),col="blue",size=1) + 
    geom_point(aes(x=t_plot-0.1,y=dat_detectable,shape="Observed"),col="grey20",size=2) + 
    scale_x_continuous(limits=c(0.5,3.5),breaks=seq(1,3),labels=xlabels) + 
    scale_y_continuous(expand=c(0,0),limits=c(0,ymax2))  +
    scale_shape_manual(values=c("Predicted"=19,"Observed"=8)) +
    scale_color_manual(values=c("Predicted"="blue","Observed"="grey20")) +
    guides(shape=guide_legend(title=NULL,override.aes=list(color=c("grey20","blue")))) +
    export_theme +
    inset_theme +
    theme(legend.position=c(0.5,0.8),axis.text.x=element_text(size=6)) +
    xlab("") +
    ylab("PCR detectable prevalence")
  p_inset2
}


