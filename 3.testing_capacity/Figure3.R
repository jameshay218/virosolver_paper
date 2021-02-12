########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(ggpubr)
library(data.table)

# setwd("~/Documents/GitHub/virosolver_paper")
plot_out <- "../virosolver_paper/figures"

## Code for plotting
source("code/plot_funcs.R")
export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                     axis.text.y=element_text(size=7),
                                     axis.title.x=element_text(size=8),
                                     axis.title.y=element_text(size=8),
                                     strip.text.x = element_text(size=7),
                                     legend.background = element_rect(fill="white",color=NA),
                                     legend.title=element_blank(),
                                     legend.text=element_text(size=5),
                                     legend.spacing.x=unit(0.05,"cm"),
                                     legend.spacing.y=unit(0.05,"cm"),
                                     plot.margin=margin(0,0,0,0,"pt"),
                                     panel.grid.major=element_line(size=0.1,color="grey40"))


## Setting Directories:
main_wd <- "../virosolver_paper/"
chainwd <- "../virosolver_paper/mcmc_chains/ReportSims/Symptom"
plot_wd <- "../virosolver_paper/plots/ReportSims/Symptom"
data_wd <- "../virosolver_paper/data/ReportSims/Symptom"
results_wd <- "../virosolver_paper/results/ReportSims/Symptom"
data_out <- "../virosolver_paper/data"

### Note: to generate the figures without re-compiling the data, 
###        run the above headers and then go to line 219


#########################################
## 2. Symptom-Based Reporting Comparisons
#########################################

## Loading SEIR dynamics used for these simulations:
load(file=paste0(data_wd,"/SEIR_dynamics.Rda"))

## Set number of simulations that were run:
SimNo <- 100

## Testing Schemes (from script to run simulations):
days <- 36
t_7 <- 1:(which.max(seir_dynamics$incidence)-2*7)
t_8 <- 1:(which.max(seir_dynamics$incidence)+2*7)

lastdays <- c(max(t_7),max(t_8))
names(lastdays) <- c("t_7","t_8")

days_R1 <- c(0)
days_R2 <- c(7,0)
days_R3 <- c(14,7,0)
days_R4 <- c(14,7,0)
days_R5 <- c(14,7,0)


## Combining results from the simulations:
Tests <- NULL
Res <- NULL
SEIR <- NULL
missing <- NULL
for (Sim in 1:SimNo) {
  for (Day in c(7,8)) {
    lastday <- lastdays[paste0("t_",Day)]
    for (Rv in c("R1","R2","R3","R4","R5")) {
      if (file.exists(paste0(results_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,".Rda"))) {
        load(paste0(results_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,".Rda"))
        R_Ests_Full$lastday <- ifelse(R_Ests_Full$t==lastday,1,0)
        Res <- bind_rows(Res, R_Ests_Full)
        SEIR <- bind_rows(SEIR, True_SEIR_sim)
        Test_in <- Cts_Full %>% dplyr::filter(sampled_time %in% (lastday - get(paste0("days_",Rv)))) %>%
          group_by(TestDay,TestProbs,sampled_time) %>% 
          summarize(Tests=length(ct_obs), Cases=sum(ct_obs < 40), MedianCt=median(ifelse(ct_obs<40,ct_obs, NA), na.rm=TRUE))
        Test_in$TestProbs <- factor(Rv, levels=c("R1","R2","R3","R4","R5"))
        Tests <- bind_rows(Tests, Test_in)
      } else {
        missing <- c(missing,paste0(Sim,"_",Day,Rv))
      }
    }
    for (pv in c("A","B","C","D","E")) {
      if (file.exists(paste0(results_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))) {
        load(paste0(results_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))
        Ests_Full$lastday <- ifelse(Ests_Full$t==lastday,1,0)
        Res <- bind_rows(Res, 
                         Ests_Full %>% dplyr::select(Sim,TestDay,TestProbs,variable,median,lower_95,lower_50,upper_50,upper_95,t,lastday))
        Test_in <- Cases_Full %>% dplyr::select(TestDay,TestProbs,confirmed_time,confirm) %>%
          rename(sampled_time=confirmed_time, Cases=confirm)
        Tests <- bind_rows(Tests, Test_in)
      } else {
        missing <- c(missing,paste0(Sim,"_",Day,pv))
      }
    }
  }
}

## List of simulation numbers/scenarios without results:
missing

## Getting results from SEIR dynamics
SEIR <- SEIR %>% dplyr::select(t,variable,median,TestDay,TestProbs)
SEIR <- distinct(SEIR)
True_Vals <- SEIR %>% filter(variable=="R") %>% dplyr::select(t,TestDay,median) %>% rename(Truth=median)

## Focus on Rt estimation results from last day for each sample:
Res_last <- as_tibble(Res %>% filter(lastday==1,variable=="R") %>% dplyr::select(-c(lastday,variable)) %>% 
                        left_join(True_Vals))
Res_last$Cvg <- ifelse(Res_last$Truth >= Res_last$lower_95 & Res_last$Truth <= Res_last$upper_95,1,0)
Res_last$Err <- Res_last$median - Res_last$Truth
Res_last$SqErr <- Res_last$Err^2
Res_last$CIW <- Res_last$upper_95-Res_last$lower_95

## Summarizing results across simulations:
MyRes <- function(x) {
  tibble(lower=quantile(x, .025, na.rm=TRUE),
         Q25=quantile(x, .25, na.rm=TRUE),
         median=quantile(x, .5, na.rm=TRUE),
         Q75=quantile(x, .75, na.rm=TRUE),
         upper=quantile(x, .975, na.rm=TRUE),
         mean=mean(x, na.rm=TRUE),
         sd=sd(x, na.rm=TRUE))
}
Res_Summ <- as_tibble(Res_last) %>% group_by(t,TestDay,TestProbs) %>% 
  summarize(MyRes(median),Cvg=mean(Cvg),MSE=mean(SqErr),MeanWidth=mean(CIW))
Res_Summ <- bind_rows(Res_Summ, SEIR %>% filter(variable=="R") %>% dplyr::select(!variable))
Res_Summ$mean <- ifelse(is.na(Res_Summ$mean), Res_Summ$median, Res_Summ$mean)
Res_Summ$label <- ifelse(Res_Summ$TestDay==7,"Growing Epidemic",
                         ifelse(Res_Summ$TestDay==8,"Declining Epidemic",NA))
Res_Summ$label_factor <- factor(Res_Summ$label,
                                levels=c("Growing Epidemic","Declining Epidemic"))
Res_Summ$type <- ifelse(Res_Summ$TestProbs %in% c("R1","R2","R3","R4","R5"),"Surveillance Ct Values",
                        ifelse(Res_Summ$TestProbs=="T","Simulation Truth","Reported Case Counts"))
Res_Summ$type2 <- ifelse(Res_Summ$TestProbs %in% c("A","R1","R2","R3"),"Flat Testing",
                         ifelse(Res_Summ$TestProbs %in% c("B","C","R4"),"Rising Testing",
                                ifelse(Res_Summ$TestProbs %in% c("D","E","R5"),"Falling Testing",
                                       ifelse(Res_Summ$TestProbs=="T","Simulation Truth",NA))))
Res_Summ$type2_factor <- factor(Res_Summ$type2,
                                levels=c("Flat Testing","Rising Testing","Falling Testing"))
Res_Summ$type3 <- ifelse(Res_Summ$TestProbs=="R1","One Surveillance Ct Sample",
                         ifelse(Res_Summ$TestProbs=="R2","Two Surveillance Ct Samples",
                                ifelse(Res_Summ$TestProbs %in% c("R3","R4","R5"),"Three Surveillance Ct Samples",
                                       ifelse(Res_Summ$TestProbs=="T","Simulation Truth",
                                              ifelse(Res_Summ$TestProbs %in% c("A","B","C","D","E"),"Reported Case Counts",NA)))))
Res_Summ$type3_factor <- factor(Res_Summ$type3,
                                levels=c("Reported Case Counts","One Sample","Two Samples","Three Samples"))
Res_Summ <- Res_Summ %>% dplyr::filter(TestProbs !="C",TestProbs !="E")
Res_Summ$xVal <- ifelse(Res_Summ$TestProbs=="A",-3,
                        ifelse(Res_Summ$TestProbs=="B",-2,
                               ifelse(Res_Summ$TestProbs=="D",-1,
                                      ifelse(Res_Summ$TestProbs=="T",0,
                                             ifelse(Res_Summ$TestProbs=="R1",1,
                                                    ifelse(Res_Summ$TestProbs=="R2",2,
                                                           ifelse(Res_Summ$TestProbs=="R3",3,
                                                                  ifelse(Res_Summ$TestProbs=="R4",4,
                                                                         ifelse(Res_Summ$TestProbs=="R5",5,NA)))))))))

Res_last2 <- Res_last %>% dplyr::filter(TestProbs !="C",TestProbs !="E")
Res_last2$xVal <- ifelse(Res_last2$TestProbs=="A",-3,
                        ifelse(Res_last2$TestProbs=="B",-2,
                               ifelse(Res_last2$TestProbs=="D",-1,
                                      ifelse(Res_last2$TestProbs=="T",0,
                                             ifelse(Res_last2$TestProbs=="R1",1,
                                                    ifelse(Res_last2$TestProbs=="R2",2,
                                                           ifelse(Res_last2$TestProbs=="R3",3,
                                                                  ifelse(Res_last2$TestProbs=="R4",4,
                                                                         ifelse(Res_last2$TestProbs=="R5",5,NA)))))))))
Res_last2$label <- ifelse(Res_last2$TestDay==7,"Growing Epidemic",
                         ifelse(Res_last2$TestDay==8,"Declining Epidemic",NA))
Res_last2$label_factor <- factor(Res_last2$label,
                                levels=c("Growing Epidemic","Declining Epidemic"))

Tests_Summ <- Tests %>% group_by(TestDay,TestProbs,sampled_time) %>% summarize(Cases=mean(Cases, na.rm=TRUE),
                                                                               Tests=mean(Tests, na.rm=TRUE))
Tests_Summ$Negatives <- Tests_Summ$Tests - Tests_Summ$Cases
Tests_Summ$label <- ifelse(Tests_Summ$TestDay==7,"Growing Epidemic",
                           ifelse(Tests_Summ$TestDay==8,"Declining Epidemic",NA))
Tests_Summ$label_factor <- factor(Tests_Summ$label,
                                  levels=c("Growing Epidemic","Declining Epidemic"))
Tests_Summ$type <- ifelse(Tests_Summ$TestProbs %in% c("R1","R2","R3","R4","R5"),"Surveillance Ct Samples",
                          "Reported Case Counts")
Tests_Summ$type2 <- ifelse(Tests_Summ$TestProbs %in% c("A","R1","R2","R3"),"Flat Testing",
                           ifelse(Tests_Summ$TestProbs %in% c("B","C","R4"),"Rising Testing",
                                  ifelse(Tests_Summ$TestProbs %in% c("D","E","R5"),"Falling Testing",NA)))
Tests_Summ$type2_factor <- factor(Tests_Summ$type2,
                                  levels=c("Flat Testing","Rising Testing","Falling Testing"))
Tests_Summ$type3 <- ifelse(Tests_Summ$TestProbs=="R1","One Surveillance Ct Sample",
                           ifelse(Tests_Summ$TestProbs=="R2","Two Surveillance Ct Samples",
                                  ifelse(Tests_Summ$TestProbs %in% c("R3","R4","R5"),"Three Surveillance Ct Samples",
                                         "Reported Case Counts")))
Tests_Summ$type3_factor <- factor(Tests_Summ$type3,
                                  levels=c("Reported Case Counts","One Surveillance Ct Sample",
                                           "Two Surveillance Ct Samples","Three Surveillance Ct Samples"))
Tests_Summ <- Tests_Summ %>% filter(TestProbs !="C", TestProbs !="E")

Tests_Summ$xVal <- Tests_Summ$sampled_time + 
  ifelse(Tests_Summ$type=="Reported Case Counts",0,
         ifelse(Tests_Summ$TestProbs=="R1",-2,
                ifelse(Tests_Summ$TestProbs=="R2",-1,
                       ifelse(Tests_Summ$TestProbs=="R3",0,
                              ifelse(Tests_Summ$TestProbs=="R4",1,
                                     ifelse(Tests_Summ$TestProbs=="R5",2,0))))))

Tests_Summ <- Tests_Summ %>% dplyr::select(-c("Tests")) %>% 
  pivot_longer(cols=c(Cases,Negatives)) %>% filter(!is.nan(value))

## Saving Data for Fig 3:
save(list=c("Tests_Summ","Res_Summ","Res_last2"), file=paste0(data_out,"/Fig3_data.Rda"))


############################
## Plotting for Figure 3: ##
############################
## Load data:
load(file=paste0(data_out,"/Fig3_data.Rda"))

sample_times <- Tests_Summ %>% ungroup() %>% group_by(label_factor) %>% 
  summarize(sampled_time=max(sampled_time)) %>% pull(sampled_time)
sample_time_dat <- data.frame(samp_times=c(sample_times[1],sample_times[1],sample_times[2]),
                              label_factor=c("Growing Epidemic","Declining Epidemic","Declining Epidemic"),
                              lty=c("A","A","B"))
sample_time_dat$label_factor <- factor(sample_time_dat$label_factor,levels=c("Growing Epidemic","Declining Epidemic"))


# conv_label_factor <- c("Symptom Testing, Flat"="Symptom Testing\nFlat", 
#                        "Symptom Testing, Rising"="Symptom Testing\nRising", 
#                        "Symptom Testing, Falling"="Symptom Testing\nFalling",
#                        "Random Sampling"="Random Sampling")
# Tests_Summ_C$type_factor2 <- conv_label_factor[as.character(Tests_Summ_C$type_factor)]
# Tests_Summ_C$type_factor2 <- factor(Tests_Summ_C$type_factor2, levels=conv_label_factor)

## Panel A: Plot of positive tests by day and testing scheme:
p.A <- ggplot(Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 110 & name=="Cases"),
              aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                  size=type)) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs != "A" & TestProbs != "D" & name=="Cases")) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs=="A"), alpha=.75) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs=="D"), alpha=.75) +
  geom_vline(data=sample_time_dat,aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor), 
             scales="free_x") +
  scale_y_continuous(name="Observed positive tests per day",breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(30,110,by=10)) +
  scale_linetype_manual(values=c("dashed","dotted")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample","Two Samples, Flat","Three Samples, Flat",
                             "Three Samples, Rising","Three Samples, Falling")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample","Two Samples, Flat","Three Samples, Flat",
                              "Three Samples, Rising","Three Samples, Falling")) +
  guides(linetype="none",size="none") +
  export_theme +
  theme(strip.text.y = element_text(size=6)) +
  labs(tag="A")
p.A

p.A.left <- ggplot(Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 70 & name=="Cases" & TestDay=="7"),
              aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                  size=type)) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs != "A" & TestProbs != "D" & name=="Cases" & TestDay=="7")) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs=="A" & TestDay=="7"), alpha=.75) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs=="D" & TestDay=="7"), alpha=.75) +
  geom_vline(data=sample_time_dat %>% filter(label_factor=="Growing Epidemic"),aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor),
             scales="free_x") +
  scale_y_continuous(name="Observed positive tests per day",
                     limits=c(0,500),
                     breaks=seq(0,500,by=100))+ #,breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(32,60,by=7)) +
  scale_linetype_manual(values=c("dashed","dotted")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample: 3000 Tests","Two Samples: 1500 Tests Each",
                             "Three Samples: 1000 Tests Each",
                             "Three Samples: 500, 1000, 1500 Tests","Three Samples: 1500, 1000, 500 Tests")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample: 3000 Tests","Two Samples: 1500 Tests Each",
                              "Three Samples: 1000 Tests Each",
                              "Three Samples: 500, 1000, 1500 Tests","Three Samples: 1500, 1000, 500 Tests")) +
  guides(linetype="none",size="none") +
  export_theme +
  theme(strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.tag=element_text(face="bold"))+
  labs(tag="A")
p.A.left

p.A.right <- ggplot(Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 100 & name=="Cases" & TestDay=="8"),
                   aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                       size=type)) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs != "A" & TestProbs != "D" & name=="Cases" & TestDay=="8")) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs=="A" & TestDay=="8"), alpha=.75) +
  geom_col(data=Tests_Summ %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs=="D" & TestDay=="8"), alpha=.75) +
  geom_vline(data=sample_time_dat %>% filter(label_factor=="Declining Epidemic"),aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor),
             scales="free_x") +
  scale_y_continuous(name="")+ #,breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(32, 88,by=7)) +
  scale_linetype_manual(values=c("dashed","dotted")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample: 3000 Tests","Two Samples: 1500 Tests Each",
                             "Three Samples: 1000 Tests Each",
                             "Three Samples: 500, 1000, 1500 Tests","Three Samples: 1500, 1000, 500 Tests")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple3",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample: 3000 Tests","Two Samples: 1500 Tests Each",
                              "Three Samples: 1000 Tests Each",
                              "Three Samples: 500, 1000, 1500 Tests","Three Samples: 1500, 1000, 500 Tests")) +
  guides(linetype="none",size="none") +
  export_theme +
  theme(strip.text.y = element_text(size=6))
p.A.right
p.A.left + p.A.right + plot_layout(nrow=1,ncol=2, guides="collect")

## Panel B: Plot of estimated Rt values, median +/- 1 SD across the simulations

p.B.meds <- ggplot(Res_Summ,
              aes(x=xVal, y=median, ymin=Q25, ymax=Q75, 
                  color=TestProbs, fill=TestProbs, shape=TestProbs, linetype=TestProbs,
                  size=TestProbs)) + 
  geom_hline(yintercept=1,col="black",size=0.25) +
  geom_hline(data=Res_Summ %>% filter(type=="Simulation Truth"),aes(yintercept=median,linetype="T"),
             col="grey40",size=0.25) +
  geom_point() + 
  geom_errorbar(aes(size="L"),width=0.2,size=0.25) + 
  export_theme + 
  theme(axis.line.x=element_blank(),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black",size=0.5),
        #legend.position=c(0.75,0.7),
        legend.position="right",
        legend.text=element_text(size=5)) +
  scale_x_discrete(name="Testing scheme and estimation method", breaks=NULL, labels=NULL) +
  scale_y_continuous(name=expression("Estimated "*R[t]))+
                     # expand=c(0,0))+
                     # limits=c(0.1,3.6),
                     # trans="log10", breaks=c(1/3,1/2,2/3,1,3/2,2,3),
                     # labels=c("0.33","0.50","0.67","1.0","1.50","2.0","3.0")) +
                     # limits=c(0,4),
                     # breaks=seq(0,4,by=.5)) +
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                              "Three Surveillance Ct Samples, Flat Testing",
                              "Three Surveillance Ct Samples, Rising Testing",
                              "Three Surveillance Ct Samples, Falling Testing")) +
  scale_size_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5","L"), 
                    values=c(1,1,1,2,1,1,1,1,1,0.5)/1.2,
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                             "Three Surveillance Ct Samples, Flat Testing",
                             "Three Surveillance Ct Samples, Rising Testing",
                             "Three Surveillance Ct Samples, Falling Testing",NA)) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                                 "Three Surveillance Ct Samples, Flat Testing",
                                 "Three Surveillance Ct Samples, Rising Testing",
                                 "Three Surveillance Ct Samples, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple3","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                              "Three Surveillance Ct Samples, Flat Testing",
                              "Three Surveillance Ct Samples, Rising Testing",
                              "Three Surveillance Ct Samples, Falling Testing")) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple3","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                             "Three Surveillance Ct Samples, Flat Testing",
                             "Three Surveillance Ct Samples, Rising Testing",
                             "Three Surveillance Ct Samples, Falling Testing")) +
  guides(size="none") +
  facet_wrap(~label_factor, scales="free_y") +
  labs(tag="B")
p.B.meds

## Panel B alternate: Plot of estimated Rt values from each simulation

p.B.dots <- ggplot(data=Res_last2, aes(x=xVal, y=median, 
                                   color=TestProbs, fill=TestProbs,shape=TestProbs)) + 
  geom_hline(yintercept=1,col="black",size=0.25) +
  geom_hline(data=Res_Summ %>% filter(type=="Simulation Truth"),aes(yintercept=median),
             linetype="dashed",col="black",size=0.25) +
  geom_point(data=Res_Summ %>% filter(type=="Simulation Truth"),
             aes(x=xVal,y=median,shape=TestProbs,color=TestProbs),
             size=1.5) +
  geom_jitter(height=0, width=.2, size=.5) +
  export_theme+
  facet_wrap(~label_factor, scales="free_y") +
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"),
                     values=c(16,16,16,8,16,16,16,16,16),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing",
                              "Three Surveillance Ct Samples, Flat Testing",
                              "Three Surveillance Ct Samples, Rising Testing",
                              "Three Surveillance Ct Samples, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple3","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                              "Three Surveillance Ct Samples, Flat Testing",
                              "Three Surveillance Ct Samples, Rising Testing",
                              "Three Surveillance Ct Samples, Falling Testing")) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple3","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample", "Two Surveillance Ct Samples, Flat Testing", 
                             "Three Surveillance Ct Samples, Flat Testing",
                             "Three Surveillance Ct Samples, Rising Testing",
                             "Three Surveillance Ct Samples, Falling Testing")) +
  scale_x_discrete(name="Testing scheme and estimation method", breaks=NULL, labels=NULL) +
  scale_y_continuous(name=expression("Median estimated "*R[t]))+
  guides(linetype="none") +
  theme(panel.grid.major.y=element_blank(),
        plot.tag=element_text(face="bold"))+
  labs(tag="B")
p.B.dots

ggsave(filename=paste0(plot_out,"/Figure3.pdf"),
       plot=(p.A.left+p.A.right+plot_layout(nrow=1,ncol=2,guides="collect"))/p.B.dots,
       width=7, height=6, units="in", dpi=300)

