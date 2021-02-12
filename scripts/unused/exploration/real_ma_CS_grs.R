########################################
## 1. Headers
########################################
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
devtools::load_all("~/Documents/GitHub/virosolver")

## Arguments for this run
set.seed(1)
n_samp <- 1000

runname_seir <- "real_ma_seir_pt"
runname_exp <- "real_ma_exp"

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- paste0("~/Documents/GitHub/virosolver_paper/mcmc_chains/5.real_ma_single_timepoint/")
plot_wd <- paste0("~/Documents/GitHub/virosolver_paper/plots/5.real_ma_single_timepoint/")
setwd(main_wd)

## MCMC parameters for Ct model fits
mcmcPars_ct_seir <- c("iterations"=100000,"popt"=0.44,"opt_freq"=5000,
                 "thin"=100,"adaptive_period"=50000,"save_block"=1000)
mcmcPars_ct_exp <- c("iterations"=200000,"popt"=0.44,"opt_freq"=5000,
                 "thin"=100,"adaptive_period"=100000,"save_block"=1000)

## Code for plotting
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

########################################
## 2. Model parameters and simulation settings
########################################
## GP model parameters for fitting
parTab_exp <- read.csv("pars/partab_exp_pos_model.csv")
parTab_seir <- read.csv("pars/partab_seir_model.csv")

########################################
## 3. Read in MA data
########################################
obs_dat_all <- read_csv("~/Documents/GitHub/ct_inference_preprint/data/BWH_COVID_Cts_deid_20200403-20200831.csv") %>%
  mutate(id=1:n())
obs_dat1 <- obs_dat_all

obs_dat1 <-  obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  mutate(first_day = as.numeric(first_day)) %>%
  mutate(first_day = first_day - min(first_day) + 35) %>% ## Start 35 days before first sample
  arrange(first_day) %>%
  rename(t = first_day, ct=panther_Ct)

obs_times <- unique(obs_dat1$t)

obs_dat_all <- obs_dat_all %>% 
  filter(platform=="Panther" &
           first_pos %in% c(1,0)) %>%
  filter(coll_date > "2020-04-15") %>% ## After biased symptomatic sampling time
  rename(date=coll_date) %>%
  left_join(epi_calendar) %>%
  dplyr::select(first_day,  panther_Ct, id) %>%
  arrange(first_day) %>%
  rename(date = first_day, ct=panther_Ct)

comb_dat <- left_join(obs_dat1, obs_dat_all)
date_key <- distinct(comb_dat %>% dplyr::select(t, date))

date_min_date <- min(date_key$date)
date_min_t <- min(date_key$t)

date_max_date <- max(date_key$date)
date_max_t <- max(date_key$t)

integer_seq_times <- seq(0, date_max_t)
date_seq_times <- seq(date_min_date-date_min_t, date_max_date,by="1 day")
date_key <- tibble(t=integer_seq_times,date=date_seq_times)

res_seir <- NULL
res_exp <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  
  ## Read in SEIR chains
  chainwd_tmp <- paste0(chainwd,"/",runname_seir,"/",timepoint)
  chain_seir <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,
                                 #mcmcPars_ct_seir["adaptive_period"],
                                 30000,
                                 multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain_seir <- as.data.frame(chain_seir)
  #chain_seir$sampno <- 1:nrow(chain_seir)
  res_seir[[i]] <- chain_seir
  
  ## Read in exp chains
  chainwd_tmp <- paste0(chainwd,"/",runname_exp,"/",timepoint)
  chain_exp <- load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct_exp["adaptive_period"],multi=TRUE,chainNo=TRUE)$chain
  chain_exp <- as.data.frame(chain_exp)
  chain_exp$sampno <- 1:nrow(chain_exp)
  res_exp[[i]] <- chain_exp
}

grs_daily <- NULL
grs_average <- NULL

prop_grow_daily <- NULL
prop_grow_average <- NULL

ps <- NULL

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use_seir <- runname_seir
  runname_use_exp <- runname_exp
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  obs_dat_seir <- obs_dat_exp <- obs_dat_use
  
  ## Observation times
  obs_dat_exp <- obs_dat_exp %>% mutate(t = t - min(t), t = t + 35)
  
  ages_exp <- 1:max(obs_dat_exp$t)
  times_exp <- 0:max(obs_dat_exp$t)
  ages_seir <- 1:max(obs_dat_seir$t)
  times_seir <- 0:max(obs_dat_seir$t)
  
  chain_seir <- res_seir[[i]]
  chain_comb_seir <- chain_seir
  chain_comb_seir$sampno <- 1:nrow(chain_comb_seir)
  chain1_seir <- chain_seir
  chain_comb_seir <- chain_comb_seir[,colnames(chain_comb_seir) != "chain"]
  
  ## Get daily growth rate
  samps <- sample(unique(chain_comb_seir$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times_seir))
  for(j in seq_along(samps)){
    trajs[j,] <- pmax(solveSEIRModel_rlsoda_wrapper(get_index_pars(chain_comb_seir, samps[j]),times_seir),0.0000001)
  }
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1[trajs1 < -0.5] <- -0.5
  trajs1[trajs1 > 0.5] <- 0.5
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower95","lower50","median","upper50","upper95","t")
  
  x <- trajs1[,ncol(trajs1)]
  
  omg <- cbind(seq(0.05,0.95,by=0.01), t(sapply(seq(0.05,0.95,by=0.01), function(y) quantile(x, c(0.5-y/2,0.5+y/2)))))
  beta_max <- 0.5
  dats_all <- NULL
  
  overall_circle_top <- data.frame(y=c(0,seq(0,beta_max,by=0.001)),x=c(0,sqrt(beta_max^2-seq(0,beta_max,by=0.001)^2)))
  overall_circle_bot <- data.frame(y=c(0,seq(-beta_max,0,by=0.001)),x=c(0,sqrt(beta_max^2-seq(-beta_max,0,by=0.001)^2)))
  beta_max1 <- beta_max*0.8
  small_overall_circle_top <- data.frame(y=c(0,seq(0,beta_max1,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(0,beta_max1,by=0.001)^2)))
  small_overall_circle_bot <- data.frame(y=c(0,seq(-beta_max1,0,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(-beta_max1,0,by=0.001)^2)))
  
  
  p <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="red",alpha=0.1) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill="green",alpha=0.1) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1)
  for(j in 1:nrow(omg)){
    alpha <- omg[j,1]
    lower <- omg[j,2]
    upper <- omg[j, 3]
    y <- seq(lower,upper,by=0.01)
    tmp <- data.frame(y=y,x=sqrt(beta_max^2 - y^2),alpha=1-alpha,quant=i)
    tmp <- bind_rows(tmp, data.frame(y=0,x=0,alpha=1-alpha,quant=i))
    dats_all[[i]] <- tmp
    p <- p + geom_polygon(data=tmp,aes(x=x,y=y),alpha=0.01,fill="blue")
  }
  
  med_segment <- quantile(x, 0.5)
  ps[[i]] <- p + scale_y_continuous(limits=c(-0.5,0.5)) +
    geom_segment(data=data.frame(y=0,yend=med_segment,x=0,xend=sqrt((beta_max^2) - med_segment^2)),
                 aes(x=x,y=y,xend=xend,yend=yend),
                 arrow=arrow(length=unit(0.1,"npc")),col="black") +
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) + 
    geom_hline(yintercept = 0,linetype="dashed",size=0.25,color="grey40") +
    coord_flip() +
    theme_void() +
    theme(axis.text=element_blank())
  
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + 
    coord_cartesian(ylim=c(-0.5,0.5))
  
  chain_exp <- res_exp[[i]]
  chain_comb_exp <- chain_exp
  chain_comb_exp$sampno <- 1:nrow(chain_comb_exp)
  chain1_exp <- chain_exp
  chain_comb_exp <- chain_comb_exp[,colnames(chain_comb_exp) != "chain"]
  
  
  grs_daily[[i]] <- trajs1_quants[nrow(trajs1_quants),]
  grs_average[[i]] <- c(quantile(chain_comb_exp$beta, c(0.025,0.25,0.5,0.75,0.975)),timepoint)
  
  prop_grow_daily[[i]] <- sum(trajs1[,ncol(trajs1)] > 0)/nrow(trajs1)
  prop_grow_average[[i]] <- sum(chain_comb_exp$beta > 0)/nrow(chain_comb_exp)
}

grs_daily_dat <- do.call("bind_rows",grs_daily)
colnames(grs_daily_dat) <-c("lower95","lower50","median","upper50","upper95","t")
grs_daily_dat$ver <- "Daily"

grs_average_dat <- do.call("bind_rows",grs_average)
colnames(grs_average_dat) <- c("lower95","lower50","median","upper50","upper95","t")
grs_average_dat$ver <- "Average"

grs_all <- bind_rows(grs_daily_dat, grs_average_dat)

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
      xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
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
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                             data = xtick, size = size, 
                                             color = color)
    
    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))
    
    
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

plot_growth_dial <- function(grs_all, t, beta_max=0.5){
  overall_circle_top <- data.frame(y=c(0,seq(0,beta_max,by=0.001)),x=c(0,sqrt(beta_max^2-seq(0,beta_max,by=0.001)^2)))
  overall_circle_bot <- data.frame(y=c(0,seq(-beta_max,0,by=0.001)),x=c(0,sqrt(beta_max^2-seq(-beta_max,0,by=0.001)^2)))
  beta_max1 <- beta_max*0.8
  small_overall_circle_top <- data.frame(y=c(0,seq(0,beta_max1,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(0,beta_max1,by=0.001)^2)))
  small_overall_circle_bot <- data.frame(y=c(0,seq(-beta_max1,0,by=0.001)),x=c(0,sqrt(beta_max1^2-seq(-beta_max1,0,by=0.001)^2)))
  
  
  ## Daily GRs
  daily_circle <- grs_all[grs_all$t == t & grs_all$ver == "Daily",]
  
  ## If positive upper quantile
  upper_est95 <- daily_circle$upper95[1]
  lower_est95 <- daily_circle$lower95[1]
  upper_est50 <- daily_circle$upper50[1]
  lower_est50 <- daily_circle$lower50[1]
  median_est <- daily_circle$median[1]
  
  ## Posterior dist
  y <- seq(lower_est95,upper_est95,by=0.0001)
  daily_circle_dat <- data.frame(y=y,x=sqrt(beta_max^2 - y^2))
  
  daily_circle_growth <- daily_circle_dat[daily_circle_dat$y > 0,]
  daily_circle_growth <- bind_rows(daily_circle_growth,data.frame(x=0,y=0))
  daily_circle_decline <- daily_circle_dat[daily_circle_dat$y < 0,]
  daily_circle_decline <- bind_rows(daily_circle_decline,data.frame(x=0,y=0))
  
  
  
  ## Posterior dist
  y <- seq(lower_est50,upper_est50,by=0.0001)
  daily_circle_dat <- data.frame(y=y,x=sqrt(beta_max^2 - y^2))
  
  daily_circle_growth1 <- daily_circle_dat[daily_circle_dat$y > 0,]
  daily_circle_growth1 <- bind_rows(daily_circle_growth1,data.frame(x=0,y=0))
  daily_circle_decline1 <- daily_circle_dat[daily_circle_dat$y < 0,]
  daily_circle_decline1 <- bind_rows(daily_circle_decline1,data.frame(x=0,y=0))
  
  
  ## Quantiles
  #daily_quantile_lines <- data.frame(x=c(0,0,0),
  #                                   xend=c(sqrt(beta_max^2-upper_est^2), sqrt(beta_max^2-lower_est^2), sqrt(beta_max^2-median_est^2)),
  #                                   y=c(0,0,0),
  #                                   yend=c(upper_est,lower_est,median_est))
  
  daily_quantile_lines <- data.frame(x=c(0),
                                     xend=c(sqrt(beta_max^2-median_est^2)),
                                     y=c(0),
                                     yend=c(median_est))
  
  ## Average GRs
  average_circle <- grs_all[grs_all$t == t & grs_all$ver == "Average",]
  
  ## If positive upper quantile
  upper_est95 <- average_circle$upper95[1]
  lower_est95 <- average_circle$lower95[1]
  upper_est50 <- average_circle$upper50[1]
  lower_est50 <- average_circle$lower50[1]
  median_est <- average_circle$median[1]
  
  ## Posterior dist
  y <- seq(lower_est95,upper_est95,by=0.0001)
  average_circle_dat <- data.frame(y=c(0,y),x=c(0,sqrt(beta_max^2 - y^2)))
  
  average_circle_growth <- average_circle_dat[average_circle_dat$y > 0,]
  average_circle_growth <- bind_rows(average_circle_growth,data.frame(x=0,y=0))
  average_circle_decline <- average_circle_dat[average_circle_dat$y < 0,]
  average_circle_decline <- bind_rows(average_circle_decline,data.frame(x=0,y=0))
  
  ## Quantiles
  #average_quantile_lines <- data.frame(x=c(0,0,0),
  #                                     xend=c(sqrt(beta_max^2-upper_est^2), sqrt(beta_max^2-lower_est^2), sqrt(beta_max^2-median_est^2)),
 #                                      y=c(0,0,0),
 #                                      yend=c(upper_est,lower_est,median_est))
  
  p1 <- ggplot() +
    geom_polygon(data=overall_circle_top, aes(x=x,y=y),fill="red",alpha=0.1) +
    geom_polygon(data=overall_circle_bot, aes(x=x,y=y),fill="green",alpha=0.1) +
    geom_polygon(data=small_overall_circle_top, aes(x=x,y=y),fill="white",alpha=1) +
    geom_polygon(data=small_overall_circle_bot, aes(x=x,y=y),fill="white",alpha=1) +
    #geom_polygon(data=overall_circle_top, aes(x=-x,y=y),fill="white",col="grey40") +
    #geom_polygon(data=overall_circle_bot, aes(x=-x,y=y),fill="white",col="grey40") +
    
    geom_polygon(data=daily_circle_growth,aes(x=x,y=y),alpha=0.1,fill="blue") +
    geom_polygon(data=daily_circle_decline,aes(x=x,y=y),alpha=0.1,fill="blue") +
    geom_polygon(data=daily_circle_growth1,aes(x=x,y=y),alpha=0.25,fill="blue") +
    geom_polygon(data=daily_circle_decline1,aes(x=x,y=y),alpha=0.25,fill="blue") +
    
    #geom_polygon(data=average_circle_growth,aes(x=-x,y=y),alpha=0.9,fill="red") +
    #geom_polygon(data=average_circle_decline,aes(x=-x,y=y),alpha=0.9,fill="blue") +
    
    geom_segment(data=daily_quantile_lines,col="blue",
                 aes(x=x,y=y,xend=xend,yend=yend), col="black", arrow=arrow(length=unit(0.1,"npc"))) +
    #geom_segment(data=average_quantile_lines,
    #             aes(x=x,y=y,xend=-xend,yend=yend), col="black", arrow=arrow(length=unit(0.03,"npc"))) +
    #geom_hline(yintercept=0) +
    geom_hline(yintercept = 0,linetype="dashed",size=0.25,color="grey40") +
    #geom_text(data=data.frame(x=-0.25,y=0.5,label="35-day average"),aes(x=x,y=y,label=label)) +
    #geom_text(data=data.frame(x=0.25,y=0.5,label="Recent daily"),aes(x=x,y=y,label=label)) +
    
    #geom_text(data=data.frame(x=-0.4,y=0.4,label="Growth"),aes(x=x,y=y,label=label)) +
    #geom_text(data=data.frame(x=0.4,y=0.4,label="Growth"),aes(x=x,y=y,label=label)) +
    #geom_text(data=data.frame(x=-0.4,y=-0.4,label="Decline"),aes(x=x,y=y,label=label)) +
    #geom_text(data=data.frame(x=0.4,y=-0.4,label="Decline"),aes(x=x,y=y,label=label)) +
    
    coord_cartesian(xlim=c(-beta_max,beta_max), ylim=c(-beta_max,beta_max)) +
    coord_flip() +
    
    #theme_geometry(seq(-beta_max,beta_max,by=0.001),seq(-beta_max,beta_max,by=0.001),
    #               ticks=11,xlab="",ylab="") +
    theme_void() +
    theme(axis.text=element_blank())
  p1
  
}
#ps <- NULL
#ps_all <- ggplot()
#for(i in seq_along(obs_times)){
#  ps[[i]] <- plot_growth_dial(grs_all, obs_times[i])
#  ps_all <- ps_all | ps[[i]]
#}
#pdf("test.pdf",width=25,height=1)
#plot(ps_all)
#dev.off()

ps_all <- ggplot()
for(i in seq_along(obs_times)){
  #ps[[i]] <- plot_growth_dial(grs_all, obs_times[i])
  ps_all <- ps_all | ps[[i]]
}
pdf("test.pdf",width=25,height=1)
plot(ps_all)
dev.off()


plot_growth_dial(grs_all, 42)


daily_growth <- do.call("c",prop_grow_daily)
average_growth <- do.call("c",prop_grow_average)
growth_dat <- data.frame(timepoint=obs_times,growth=daily_growth,decline=1-daily_growth,ver="Daily")
average_dat <- data.frame(timepoint=obs_times,growth=average_growth,decline=1-average_growth,ver="Average")
all_growth_dat <- bind_rows(growth_dat, average_dat)
all_growth_dat$ver <- factor(all_growth_dat$ver, levels=c("Daily","Average"))

p1 <- ggplot(all_growth_dat) + 
  geom_segment(x=0.25,xend=0.25,y=0,yend=0.5,aes(size=decline,alpha=decline),
               col="red",linejoin="mitre",
               arrow=arrow(length=unit(0.30,"cm"),ends="first",type="closed")) +
  geom_segment(x=1.75,xend=1.75,y=0,yend=0.5,aes(size=growth,alpha=growth),
               col="green",linejoin="mitre",
               arrow=arrow(length=unit(0.30,"cm"),ends="last",type="closed")) +
  scale_size_continuous(limits=c(0,1),range=c(0,2)) +
  scale_x_continuous(limits=c(-1,3)) +
  scale_y_continuous(limits=c(-0.1,0.6)) +
  facet_grid(ver~timepoint) +
  theme_void() +
  theme(legend.position="none",
        panel.background = element_rect())
p1


x1 <- data.frame(x=0,y=0,xend=sqrt(0.5^2 - (daily_growth-0.5)^2), yend=(daily_growth-0.5),timepoint=obs_times,ver="Daily")
x2 <- data.frame(x=1,y=0,xend=sqrt(0.5^2 - (average_growth-0.5)^2)+1, yend=(average_growth-0.5),timepoint=obs_times,ver="Average")
x_all <- bind_rows(x1, x2)
ggplot(x_all) +
  geom_segment(aes(x=x,xend=xend,y=y,yend=yend, col=ver,group=timepoint), arrow=arrow(length=unit(0.3,"npc"))) +
  scale_y_continuous(limits=c(-0.5,0.5)) +
  #scale_x_continuous(limits=c(0,1)) +
  facet_wrap(~timepoint,nrow=1) +
  theme_void()

