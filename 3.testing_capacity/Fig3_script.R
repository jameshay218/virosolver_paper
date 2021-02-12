########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(ggpubr)
# library(data.table)

setwd("..")
data_out <- "../virosolver_paper/data"
plot_out <- "../virosolver_paper/figures"

## Code for plotting
source("code/plot_funcs.R")
export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                     axis.text.y=element_text(size=7),
                                     axis.title.x=element_text(size=8),
                                     axis.title.y=element_text(size=8))


### Note: to generate the figures without re-compiling the data, 
###        run the above headers and then go to line 267

##################################################################
## 2. Comparison with Changing Random-Sample Testing Probabilities
##################################################################

## Setting Directories:
main_wd <- "../virosolver_paper/"
chainwd <- "../virosolver_paper/mcmc_chains/ReportSims/ChangingTests"
plot_wd <- "../virosolver_paper/plots/ReportSims/ChangingTests"
data_wd <- "../virosolver_paper/data/ReportSims/ChangingTests/virosolver"

## Loading SEIR dynamics used for these simulations:
load(file=paste0(data_wd,"/SEIR_dynamics.Rda"))

## Set number of simulations that were run:
SimNo <- 100

## Testing Schemes:
days <- 36
t_4 <- seq(from=98-days+1, by=1, length.out=days) ## the start days depend on those used in the simulation
t_5 <- seq(from=126-days+1, by=1, length.out=days)
t_6 <- seq(from=154-days+1, by=1, length.out=days)

p_E <- rep(0.002, days)
p_F <- seq(from=0.001, to=0.003, length.out=days)
p_G <- exp(seq(from=log(0.0005), to=log(0.0041), length.out=days))
p_H <- seq(from=0.003, to=0.001, length.out=days)

## Combine results from the simulations:
Res_Full <- NULL
Pos_Counts_Full <- NULL
for (Sim in 1:SimNo) {
  for (time in 4:6) {
    for (type in c("E","F","G","H")) {
      if(file.exists(paste0(data_wd,"/Res_Comb_Sim_popno",Sim,"_",time,type,".csv"))) {
        Res_Comb <- read_csv(paste0(data_wd,"/Res_Comb_Sim_popno",Sim,"_",time,type,".csv"),
                           col_types="ddcdcdcdd")
        Res_Comb$Sim <- Sim
        if (is.null(Res_Full)) {
          Res_Full <- bind_rows(Res_Full, Res_Comb)
        } else if (dim(Res_Full %>% filter(Sim==Sim, TestDay==time))[1]==0) {
          Res_Full <- bind_rows(Res_Full, Res_Comb)
        } else {
          Res_Full <- bind_rows(Res_Full, Res_Comb %>% dplyr::filter(TestProbs==type))
        }
      }
      if (file.exists(paste0(data_wd,"/Pos_Counts_Sim_popno",Sim,"_",time,type,".csv"))) {
        Pos_Counts <- read_csv(paste0(data_wd,"/Pos_Counts_Sim_popno",Sim,"_",time,type,".csv"),
                               col_types="ddccddd")
        Pos_Counts$Sim <- Sim
        if (is.null(Pos_Counts_Full)) {
          Pos_Counts_Full <- bind_rows(Pos_Counts_Full, Pos_Counts)
        } else if (dim(Pos_Counts_Full %>% filter(Sim==Sim, TestDay==time))[1]==0) {
          Pos_Counts_Full <- bind_rows(Pos_Counts_Full, Pos_Counts)
        } else {
          Pos_Counts_Full <- bind_rows(Pos_Counts_Full, Pos_Counts %>% dplyr::filter(TestProbs==type))
        }
      }
    }
  }
}

## Summarize results from simulations for Panels A and B, with changing random-sample testing rates:
MyRes <- function(x) {
  tibble(lower=quantile(x, .025, na.rm=TRUE), 
         median=quantile(x, .5, na.rm=TRUE), 
         upper=quantile(x, .975, na.rm=TRUE),
         mean=mean(x, na.rm=TRUE),
         sd=sd(x, na.rm=TRUE))
}
Res_Summ <- Res_Full %>% dplyr::filter(variable %in% c("Est_Pos","GR35_ObsCases","GR35")) %>% 
  dplyr::filter(TestProbs != "G") %>% dplyr::filter(TestProbs != "Cases") %>% 
  group_by(t,TestDay,TestProbs,variable) %>% summarize(MyRes(median))
Pos_Counts_Summ <- Pos_Counts_Full %>% group_by(TestDay,TestProbs,sampled_time) %>%
  summarize(Tests=mean(Tests), Cases=mean(Positive))
Res_Summ$xVal <- Res_Summ$t + ifelse(Res_Summ$variable=="GR35",10,
                                               ifelse(Res_Summ$TestProbs=="E",-8,
                                                      ifelse(Res_Summ$TestProbs=="F",-4,
                                                             ifelse(Res_Summ$TestProbs=="G",0,
                                                                    ifelse(Res_Summ$TestProbs=="H",4,
                                                                           ifelse(Res_Summ$TestProbs=="Cases",8,0))))))+ifelse(Res_Summ$variable=="Est_Pos",1,
                                                                                                                                ifelse(Res_Summ$variable=="GR35_ObsCases",2,0))
Res_Summ$label <- ifelse(Res_Summ$TestDay==4,"Pre-Peak",
                                ifelse(Res_Summ$TestDay==5,"During Peak",
                                       ifelse(Res_Summ$TestDay==6,"Post-Peak",NA)))
Res_Summ$label_factor <- factor(Res_Summ$label,
                                       levels=c("Pre-Peak","During Peak","Post-Peak"))
Pos_Counts_Summ$label <- ifelse(Pos_Counts_Summ$TestDay==4,"Pre-Peak",
                                ifelse(Pos_Counts_Summ$TestDay==5,"During Peak",
                                       ifelse(Pos_Counts_Summ$TestDay==6,"Post-Peak",NA)))
Pos_Counts_Summ$label_factor <- factor(Pos_Counts_Summ$label,
                                       levels=c("Pre-Peak","During Peak","Post-Peak"))
Pos_Counts_Summ2 <- Pos_Counts_Summ %>% filter(TestProbs != "True Incidence")

Res_Summ$xVal <- 0+ifelse(Res_Summ$TestProbs=="E",-3,
                          ifelse(Res_Summ$TestProbs=="F",-2,
                                 ifelse(Res_Summ$TestProbs=="H",-1,0)))*ifelse(Res_Summ$variable=="Est_Pos",-1,1)
Res_Summ$Mthd <- ifelse(Res_Summ$TestProbs=="SEIR","T",
                        ifelse(Res_Summ$TestProbs=="E","Flat_",
                               ifelse(Res_Summ$TestProbs=="F","Rise_",
                                      ifelse(Res_Summ$TestProbs=="H","Fall_",NA))))
Res_Summ$Mthd <- ifelse(Res_Summ$variable=="Est_Pos",paste0(Res_Summ$Mthd,"PCR"),
                        ifelse(Res_Summ$variable=="GR35_ObsCases",paste0(Res_Summ$Mthd,"Cases"),Res_Summ$Mthd))

Pos_Counts_Summ3 <- bind_rows(Pos_Counts_Summ2 %>% select(!Tests))
Pos_Counts_Summ3 <- Pos_Counts_Summ3 %>% dplyr::filter(TestProbs !="G")
Pos_Counts_Summ3$type_factor <- factor(ifelse(Pos_Counts_Summ3$TestProbs=="E","Flat Testing",
                                              ifelse(Pos_Counts_Summ3$TestProbs=="F","Rising Testing",
                                                     ifelse(Pos_Counts_Summ3$TestProbs=="H","Falling Testing",NA))),
                                       levels=c("Flat Testing","Rising Testing","Falling Testing"))

## Saving Data for Figs 3A and 3B:
Pos_Counts_Summ_A <- Pos_Counts_Summ3
save(Pos_Counts_Summ_A, file=paste0(data_out,"/Fig3A_data.Rda"))
Res_Summ_B <- Res_Summ
save(Res_Summ_B, file=paste0(data_out,"/Fig3B_data.Rda"))

#########################################
## 3. Symptom-Based Reporting Comparisons
#########################################

## Setting Directories:
chainwd <- "../virosolver_paper/mcmc_chains/ReportSims/ChangingTests"
plot_wd <- "../virosolver_paper/plots/ReportSims/Symptom"
data_wd <- "../virosolver_paper/data/ReportSims/Symptom"
data_out <- "../virosolver_paper/data"

## Loading SEIR dynamics used for these simulations:
load(file=paste0(data_wd,"/SEIR_dynamics.Rda"))

## Set number of simulations that were run:
SimNo <- 100

## Symptom-Based Testing Schemes:
days <- 36

t_7 <- 1:(92-2*7) ## depends on what was used in simulations
t_8 <- 1:(92+2*7)
t_9 <- 1:(92+6*7)

lastdays <- c(max(t_7),max(t_8),max(t_9))
names(lastdays) <- c("t_7","t_8","t_9")

p_A <- rep(.1, length.out=days)
p_B <- seq(from=.1,to=.2, length.out=days)
p_C <- exp(seq(from=log(.1),to=log(.25),length.out=days))
p_D <- seq(from=.1,to=.01, length.out=days)
p_E <- 0.2-exp(seq(from=log(.1),to=log(.195),length.out=days))

p_R <- c(0.003,0.003,0.003)
days_R1 <- c(0) # On last day
days_R2 <- c(7,0) # On last day and week before
days_R3 <- c(14,7,0) # On last day, week before and 2 weeks before

## Combining results from the simulations
Tests <- NULL
Res <- NULL
SEIR <- NULL
for (Sim in 1:SimNo) {
  for (Day in c(7,8,9)) {
    lastday <- lastdays[paste0("t_",Day)]
    for (Rv in c("R1","R2","R3")) {
      if (file.exists(paste0(data_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,".Rda"))) {
        load(paste0(data_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,".Rda"))
        R_Ests_Full$lastday <- ifelse(R_Ests_Full$t==lastday,1,0)
        Res <- bind_rows(Res, R_Ests_Full)
        SEIR <- bind_rows(SEIR, True_SEIR_sim)
        Test_in <- Cts_Full %>% dplyr::filter(sampled_time %in% (lastday - get(paste0("days_",Rv)))) %>%
          group_by(TestDay,TestProbs,sampled_time) %>% 
          summarize(Tests=length(ct_obs), Cases=sum(ct_obs < 40), MedianCt=median(ifelse(ct_obs<40,ct_obs, NA), na.rm=TRUE))
        Test_in$TestProbs <- factor(Rv, levels=c("R1","R2","R3"))
        Tests <- bind_rows(Tests, Test_in)
      }
    }
    for (pv in c("A","B","C","D","E")) {
      if (file.exists(paste0(data_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))) {
        load(paste0(data_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))
        Ests_Full$lastday <- ifelse(Ests_Full$t==lastday,1,0)
        Res <- bind_rows(Res, 
                         Ests_Full %>% dplyr::select(Sim,TestDay,TestProbs,variable,median,lower_95,lower_50,upper_50,upper_95,t,lastday))
        Test_in <- Cases_Full %>% dplyr::select(TestDay,TestProbs,confirmed_time,confirm) %>%
          rename(sampled_time=confirmed_time, Cases=confirm)
        Tests <- bind_rows(Tests, Test_in)
      }
    }
  }
}

## Getting results from SEIR dynamics
SEIR <- SEIR %>% dplyr::select(t,variable,median,TestDay,TestProbs)
SEIR <- distinct(SEIR)

## Summarizing results across simulations:
MyRes <- function(x) {
  tibble(lower=quantile(x, .025, na.rm=TRUE), 
         median=quantile(x, .5, na.rm=TRUE), 
         upper=quantile(x, .975, na.rm=TRUE),
         mean=mean(x, na.rm=TRUE),
         sd=sd(x, na.rm=TRUE))
}
Res_Summ <- as_tibble(Res) %>% filter(variable=="R", lastday==1) %>%
  group_by(t,TestDay,TestProbs) %>% summarize(MyRes(median))
Res_Summ <- bind_rows(Res_Summ, SEIR %>% filter(variable=="R") %>% select(!variable))
Res_Summ$mean <- ifelse(is.na(Res_Summ$mean), Res_Summ$median, Res_Summ$mean)
Res_Summ$label <- ifelse(Res_Summ$TestDay==7,"Growing Epidemic",
                                ifelse(Res_Summ$TestDay==8,"Declining Epidemic",
                                       ifelse(Res_Summ$TestDay==9,"Ending Epidemic",NA)))
Res_Summ$label_factor <- factor(Res_Summ$label,
                                       levels=c("Growing Epidemic","Declining Epidemic","Ending Epidemic"))
Res_Summ$type <- ifelse(Res_Summ$TestProbs %in% c("R1","R2","R3"),"Random Sampling",
                        ifelse(Res_Summ$TestProbs=="T","Simulation Truth","Symptom-Based Testing"))
Res_Summ <- Res_Summ %>% dplyr::filter(TestProbs !="C",TestProbs !="E")
Res_Summ$xVal <- ifelse(Res_Summ$TestProbs=="A",-3,
                        ifelse(Res_Summ$TestProbs=="B",-2,
                               ifelse(Res_Summ$TestProbs=="D",-1,
                                      ifelse(Res_Summ$TestProbs=="T",0,
                                             ifelse(Res_Summ$TestProbs=="R1",1,
                                                    ifelse(Res_Summ$TestProbs=="R2",2,
                                                           ifelse(Res_Summ$TestProbs=="R3",3,NA)))))))

Tests_Summ <- Tests %>% group_by(TestDay,TestProbs,sampled_time) %>% summarize(Cases=mean(Cases, na.rm=TRUE))
Tests_Summ$label <- ifelse(Tests_Summ$TestDay==7,"Growing Epidemic",
                           ifelse(Tests_Summ$TestDay==8,"Declining Epidemic",
                                  ifelse(Tests_Summ$TestDay==9,"Ending Epidemic",NA)))
Tests_Summ$label_factor <- factor(Tests_Summ$label,
                                  levels=c("Growing Epidemic","Declining Epidemic","Ending Epidemic"))
Tests_R <- Tests_Summ %>% dplyr::filter(TestProbs %in% c("R1","R2","R3")) %>% 
  group_by(TestDay,sampled_time,label,label_factor) %>% 
  summarize(Cases=mean(Cases, na.rm=TRUE))
Tests_R_times <- Tests_Summ %>% filter(TestProbs=="A") %>% ungroup %>% select(TestDay,sampled_time,label,label_factor)
Tests_R2 <- full_join(Tests_R, Tests_R_times)
Tests_R2$TestProbs <- "PCR"
Tests_R2$Cases <- ifelse(is.na(Tests_R2$Cases), 0, Tests_R2$Cases)
Tests_Summ2 <- bind_rows(Tests_Summ %>% 
                           dplyr::filter(TestProbs %in% c("A","B","D")), 
                         Tests_R2)
Tests_Summ2$type_factor <- factor(ifelse(Tests_Summ2$TestProbs=="A","Symptom Testing, Flat",
                                         ifelse(Tests_Summ2$TestProbs=="B","Symptom Testing, Rising",
                                                ifelse(Tests_Summ2$TestProbs=="D","Symptom Testing, Falling",
                                                       ifelse(Tests_Summ2$TestProbs=="PCR","Random Sampling",NA)))),
                                  levels=c("Symptom Testing, Flat","Symptom Testing, Rising",
                                           "Symptom Testing, Falling","Random Sampling"))

## Saving Data for Figs 3C and 3D:
Tests_Summ_C <- Tests_Summ2 %>% dplyr::filter(TestDay != 9)
save(Tests_Summ_C, file=paste0(data_out,"/Fig3C_data.Rda"))
Res_Summ_D <- Res_Summ %>% dplyr::filter(TestDay != 9)
save(Res_Summ_D, file=paste0(data_out,"/Fig3D_data.Rda"))


############################
## Plotting for Figure 3: ##
############################
## Load data for Panels A and B:
load(file=paste0(data_out,"/Fig3A_data.Rda"))
load(file=paste0(data_out,"/Fig3B_data.Rda"))

## Panel A: Plot of Positive Test Results by sampling scheme:
p.A <- ggplot(Pos_Counts_Summ_A %>% dplyr::filter(TestProbs !="G"), 
              aes(x=sampled_time, y=Cases, fill=TestProbs)) +
  geom_col() +
  facet_grid(rows=vars(type_factor), cols=vars(label_factor), scales="free_x") +
  scale_y_continuous(name="Positive Tests Per Day") +
  scale_x_continuous(name="Test Day") +
  scale_fill_manual(name="Testing Scheme",
                    breaks=c("E","F","H"), values=unname(AAAS_palette[c("red2","teal1","purple1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing")) +
  guides(fill="none") +
  export_theme

## Panel B: Plot of estimated growth rates, with mean and +/- 1 SD across simulations:
p.B <- ggplot(data=Res_Summ_B,
              aes(x=xVal, y=mean, ymin=mean-sd, ymax=mean+sd, color=Mthd, shape=Mthd, size=Mthd, linetype=Mthd)) + 
  geom_point() +
  geom_errorbar(aes(size="L")) +
  scale_x_continuous(name="Testing Scheme and Estimation Method", breaks=NULL, labels=NULL) +
  scale_y_continuous(name="Estimated 35-Day Growth Rate") +
  scale_shape_manual(name="Testing Scheme and Estimation Method", 
                     breaks=c("Flat_Cases","Rise_Cases","Fall_Cases","T","Flat_PCR","Rise_PCR","Fall_PCR"),
                     values=c(15,15,15,19,16,16,16), 
                     labels=c("Case Counts: Flat Testing","Case Counts: Rising Testing","Case Counts: Falling Testing",
                              "Simulation Truth",
                              "PCR Ct Values: Flat Testing","PCR Ct Values: Rising Testing","PCR Ct Values: Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method", 
                     breaks=c("Flat_Cases","Rise_Cases","Fall_Cases","T","Flat_PCR","Rise_PCR","Fall_PCR"),
                     values=unname(AAAS_palette[c("red2","teal1","purple1","black","red2","teal1","purple1")]), 
                     labels=c("Case Counts: Flat Testing","Case Counts: Rising Testing","Case Counts: Falling Testing",
                              "Simulation Truth",
                              "PCR Ct Values: Flat Testing","PCR Ct Values: Rising Testing","PCR Ct Values: Falling Testing")) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("Flat_Cases","Rise_Cases","Fall_Cases","T","Flat_PCR","Rise_PCR","Fall_PCR"),
                        values=c(1,1,1,0,1,1,1),
                        labels=c("Case Counts: Flat Testing","Case Counts: Rising Testing","Case Counts: Falling Testing",
                                 "Simulation Truth",
                                 "PCR Ct Values: Flat Testing","PCR Ct Values: Rising Testing","PCR Ct Values: Falling Testing")) +
  scale_size_manual(name="Testing Scheme and Estimation Method", 
                    breaks=c("Flat_Cases","Rise_Cases","Fall_Cases","T","Flat_PCR","Rise_PCR","Fall_PCR","L"),
                    values=c(2,2,2,4,2,2,2,1), 
                    labels=c("Case Counts: Flat Testing","Case Counts: Rising Testing","Case Counts: Falling Testing",
                             "Simulation Truth",
                             "PCR Ct Values: Flat Testing","PCR Ct Values: Rising Testing","PCR Ct Values: Falling Testing",NA)) +
  guides(size="none")+
  export_theme +
  facet_wrap(facets=vars(label_factor), nrow=1, ncol=3, scales="free_x")


## Load data for Panels C and D:
load(file=paste0(data_out,"/Fig3C_data.Rda"))
load(file=paste0(data_out,"/Fig3D_data.Rda"))

## Panel C: Plot of positive tests by day and testing scheme:
p.C <- ggplot(Tests_Summ_C,
              aes(x=sampled_time, y=Cases, fill=TestProbs)) +
  geom_col() +
  facet_grid(rows=vars(type_factor), cols=vars(label_factor), scales="free_x") +
  scale_y_continuous(name="Positive Tests Per Day") +
  scale_x_continuous(name="Test Day") +
  scale_fill_manual(name="Testing Scheme",
                    breaks=c("A","B","D","PCR"), 
                    values=unname(AAAS_palette[c("red1","green1","purple3","blue1")]),
                    labels=c("Symptom Testing, Flat","Symptom Testing, Rising",
                             "Symptom Testing, Falling","Random Sampling")) +
  guides(fill="none") +
  export_theme

## Panel D: Plot of estimated Rt values, mean +/- 1 SD across the simulations
p.D <- ggplot(Res_Summ_D,
              aes(x=xVal, y=mean, ymin=mean-sd, ymax=mean+sd, 
                  color=TestProbs, shape=TestProbs, alpha=TestProbs, linetype=TestProbs,
                  size=TestProbs)) + 
  geom_point() + geom_errorbar(aes(size="L")) + 
  export_theme + 
  scale_x_discrete(name="Testing Scheme and Estimation Method", breaks=NULL, labels=NULL) +
  scale_y_continuous(name=expression("Estimated "*R[t])) +
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3"), 
                     values=c(16,16,16,19,15,15,15),
                     labels=c("Symptom-Based Testing, Flat","Symptom-Based Testing, Rising",
                              "Symptom-Based Testing, Falling","Simulation Truth",
                              "Random Sampling: One Sample", "Random Sampling: Two Samples", 
                              "Random Sampling: Three Samples")) +
  scale_size_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","L"), 
                    values=c(2,2,2,4,2,2,2,1),
                    labels=c("Symptom-Based Testing, Flat","Symptom-Based Testing, Rising",
                             "Symptom-Based Testing, Falling","Simulation Truth",
                             "Random Sampling: One Sample", "Random Sampling: Two Samples", 
                             "Random Sampling: Three Samples",NA)) +
  guides(size="none") +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3"), 
                        values=c(1,1,1,0,1,1,1),
                        labels=c("Symptom-Based Testing, Flat","Symptom-Based Testing, Rising",
                                 "Symptom-Based Testing, Falling","Simulation Truth",
                                 "Random Sampling: One Sample", "Random Sampling: Two Samples", 
                                 "Random Sampling: Three Samples")) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3"), 
                     values=c(1,1,1,1,0.33,0.67,1),
                     labels=c("Symptom-Based Testing, Flat","Symptom-Based Testing, Rising",
                              "Symptom-Based Testing, Falling","Simulation Truth",
                              "Random Sampling: One Sample", "Random Sampling: Two Samples", 
                              "Random Sampling: Three Samples")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3"), 
                     values=unname(AAAS_palette[c("red1","green1","purple3","black",rep("blue1",3))]),
                     labels=c("Symptom-Based Testing, Flat","Symptom-Based Testing, Rising",
                              "Symptom-Based Testing, Falling","Simulation Truth",
                              "Random Sampling: One Sample", "Random Sampling: Two Samples", 
                              "Random Sampling: Three Samples")) +
  facet_wrap(~label_factor, scales="free")

## Combined plot with both scenarios:
ggsave(filename=paste0(plot_out,"/Fig3_Testing_Changes.png"),
       plot=(p.A+labs(tag="A)") | p.B+labs(tag="B)")) / (p.C+labs(tag="C)") | p.D+labs(tag="D)")), 
       width=15, height=12, units="in", dpi=300)

