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
results_wd <- "../virosolver_paper/results/ReportSims/Symptom_SS"
data_out <- "../virosolver_paper/data"

## Values of Testing Probability
test_probs <- c("1e-04"=100, "2e-04"=200, "5e-04"=500,
                "0001"=1000, "0002"=2000, "0003"=3000)

### Note: to generate the figures without re-compiling the data, 
###        run the above headers and then go to line 264


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

control_table <- expand.grid(pop_no=1:100,strategy=c("R1","R2","R3","R4","R5"),
                             testing_day=as.character(c(7,8)),total_prob=c(0.0001,0.0002,0.0005,0.001,0.002,0.003)) %>% 
  arrange(pop_no, testing_day, strategy)
control_table$run_name <- 1:nrow(control_table)


## Combining results from the simulations:
Tests <- NULL
Res <- NULL
SEIR <- NULL
missing <- NULL
missing_runs <- NULL
for (Sim in 1:SimNo) {
  for (Day in c(7,8)) {
    print(paste0(Sim,", ",Day))
    lastday_num <- lastdays[paste0("t_",Day)]
    for (Rv in c("R1","R2","R3","R4","R5")) {
      for (tprob in names(test_probs)) {
        if (file.exists(paste0(results_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,tprob,".Rda"))) {
          load(paste0(results_wd,"/virosolver/Sim_viro_pop",Sim,"_",Day,Rv,tprob,".Rda"))
          R_Ests_Full <- R_Ests_Full %>% mutate(lastday=if_else(t==lastday_num,1,0),
                                                analysis="Ct",
                                                size=test_probs[tprob])
          Res <- bind_rows(Res, R_Ests_Full)
          SEEIRR_Res <- SEEIRR_Res %>% mutate(lastday=if_else(t==lastday_num,1,0),
                                              analysis="Pos_SEEIRR",
                                              size=test_probs[tprob]) %>% select(-c(model))
          Res <- bind_rows(Res, SEEIRR_Res)
          SEIR_Res <- SEIR_Res %>% mutate(lastday=if_else(t==lastday_num,1,0),
                                              analysis="Pos_SEIR",
                                              size=test_probs[tprob]) %>% select(-c(model))
          Res <- bind_rows(Res, SEIR_Res)
          SEIR <- bind_rows(SEIR, True_SEIR_sim)
          Test_in <- Cts_Full %>% dplyr::filter(sampled_time %in% (lastday_num - get(paste0("days_",Rv)))) %>%
            group_by(TestDay,TestProbs,sampled_time) %>% 
            summarize(Tests=length(ct_obs), Cases=sum(ct_obs < 40), MedianCt=median(ifelse(ct_obs<40,ct_obs, NA), na.rm=TRUE)) %>%
            mutate(size=test_probs[tprob])
          Test_in$TestProbs <- factor(Rv, levels=c("R1","R2","R3","R4","R5"))
          Tests <- bind_rows(Tests, Test_in)
        } else {
          missing <- c(missing,paste0(Sim,"_",Day,Rv,tprob))
          missing_runs <- c(missing_runs, control_table$run_name[control_table$pop_no==Sim & control_table$strategy==Rv & 
                                                                   control_table$testing_day==as.character(Day) & 
                                                                   control_table$total_prob==test_probs[tprob]/1000000])
        }
      }
    }
    for (pv in c("A","B","C","D","E")) {
      if (file.exists(paste0(results_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))) {
        load(paste0(results_wd,"/EpiNow2/Sim_EN2_pop",Sim,"_",Day,pv,".Rda"))
        Ests_Full <- Ests_Full %>% mutate(lastday=if_else(t==lastday_num,1,0),
                                          analysis="Cases") %>%
          dplyr::select(Sim,TestDay,TestProbs,variable,median,lower_95,lower_50,upper_50,upper_95,t,lastday,analysis)
        Res <- bind_rows(Res, Ests_Full)
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
missing_runs

## Saving Combined Results Before Processing
save(list=c("Tests","Res","SEIR","missing","missing_runs"), file=paste0(data_out,"/Fig3_data_preproc.Rda"))

## Getting results from SEIR dynamics
SEIR <- SEIR %>% dplyr::select(t,variable,median,TestDay,TestProbs)
SEIR <- distinct(SEIR)
True_Vals <- SEIR %>% filter(variable=="R") %>% dplyr::select(t,TestDay,median) %>% rename(Truth=median)

## Focus on Rt estimation results from last day for each sample:
Res_last <- Res %>%
  filter(lastday==1, variable=="R") %>% dplyr::select(-c(lastday,variable)) %>%
  left_join(True_Vals) %>%
  mutate(Cvg=if_else(Truth >= lower_95 & Truth <= upper_95, 1, 0),
         Err=median-Truth, SqErr=(median-Truth)^2, CIW=upper_95-lower_95,
         CorrSign=if_else(Truth > 1 & median > 1, 1, if_else(Truth < 1 & median < 1, 1, 0)), 
         CorrSignCI=if_else(Truth > 1 & lower_95 > 1, 1,
                            if_else(Truth < 1 & upper_95 < 1, 1, 0)))

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
Res_Summ <- as_tibble(Res_last) %>% group_by(t,TestDay,TestProbs,analysis,size) %>% 
  summarize(MyRes(median),Cvg=mean(Cvg),MSE=mean(SqErr),MeanWidth=mean(CIW),
            EstimateSign=mean(CorrSign), CrIPower=mean(CorrSignCI))
Res_Summ <- bind_rows(Res_Summ, SEIR %>% filter(variable=="R") %>% dplyr::select(!variable) %>% mutate(analysis="Truth"))
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
Res_Summ$xVal <- sapply(Res_Summ$TestProbs,
                        FUN=function(x) switch(x,
                                               A=-2,
                                               B=-1.5,
                                               D=-1,
                                               T=0,
                                               R1=1,
                                               R2=1.5,
                                               R3=2,
                                               R4=2.5,
                                               R5=3)) +
  if_else(substr(Res_Summ$analysis,1,3)=="Pos",3,0)

Res_last2 <- Res_last %>% dplyr::filter(TestProbs !="C",TestProbs !="E")
Res_last2$xVal <- sapply(as.character(Res_last2$TestProbs),
                         FUN=function(x) switch(x,
                                                A=-2,
                                                B=-1.5,
                                                D=-1,
                                                T=0,
                                                R1=1,
                                                R2=1.5,
                                                R3=2,
                                                R4=2.5,
                                                R5=3)) +
  if_else(substr(Res_last2$analysis,1,3)=="Pos",3,0)

Res_last2$label <- ifelse(Res_last2$TestDay==7,"Growing Epidemic",
                         ifelse(Res_last2$TestDay==8,"Declining Epidemic",NA))
Res_last2$label_factor <- factor(Res_last2$label,
                                levels=c("Growing Epidemic","Declining Epidemic"))

Tests_Summ <- Tests %>% group_by(TestDay,TestProbs,sampled_time,size) %>% 
  summarize(Cases=mean(Cases, na.rm=TRUE), Tests=mean(Tests, na.rm=TRUE))
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

## Panel A: Plot of positive tests by day and testing scheme:
Test_Data_A <- Tests_Summ %>% filter(type3=="Reported Case Counts" | size==3000)
p.A <- ggplot(Test_Data_A %>% 
                filter(sampled_time > 30 & sampled_time < 110 & name=="Cases"),
              aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                  size=type)) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs != "A" & TestProbs != "D" & name=="Cases")) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs=="A"), alpha=.75) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 110 & TestProbs=="D"), alpha=.75) +
  geom_vline(data=sample_time_dat,aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor), 
             scales="free_x") +
  scale_y_continuous(name="Observed positive tests per day",breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(30,110,by=10)) +
  scale_linetype_manual(values=c("dashed","dotted")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample Period","Two Sample Periods, Flat","Three Sample Periods, Flat",
                             "Three Sample Periods, Rising","Three Sample Periods, Falling")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample Period","Two Sample Periods, Flat","Three Sample Periods, Flat",
                              "Three Sample Periods, Rising","Three Sample Periods, Falling")) +
  guides(linetype="none",size="none") +
  export_theme +
  theme(strip.text.y = element_text(size=6)) +
  labs(tag="A")
p.A

p.A.left <- ggplot(Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 70 & name=="Cases" & TestDay=="7"),
              aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                  size=type)) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs != "A" & TestProbs != "D" & name=="Cases" & TestDay=="7")) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs=="A" & TestDay=="7"), alpha=.75) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 70 & TestProbs=="D" & TestDay=="7"), alpha=.75) +
  geom_vline(data=sample_time_dat %>% filter(label_factor=="Growing Epidemic"),aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor),
             scales="free_x") +
  scale_y_continuous(name="Observed positive tests per day",
                     limits=c(0,500),
                     breaks=seq(0,500,by=100))+ #,breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(32,60,by=7)) +
  scale_linetype_manual(values=c("dashed","dotted"),
                        labels=c("Analysis 1: Epidemic Day 60",
                                 "Analysis 2: Epidemic Day 88")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample Period: 3000 Tests","Two Sample Periods: 1500 Tests Each",
                             "Three Sample Periods: 1000 Tests Each",
                             "Three Sample Periods: 500, 1000, 1500 Tests","Three Sample Periods: 1500, 1000, 500 Tests")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample Period: 3000 Tests","Two Sample Periods: 1500 Tests Each",
                              "Three Sample Periods: 1000 Tests Each",
                              "Three Sample Periods: 500, 1000, 1500 Tests","Three Sample Periods: 1500, 1000, 500 Tests")) +
  guides(linetype="none",size="none") +
  export_theme +
  theme(strip.text.y = element_blank(),
        strip.background = element_blank(),
        plot.tag=element_text(face="bold"))+
  labs(tag="A")
p.A.left

p.A.right <- ggplot(Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 100 & name=="Cases" & TestDay=="8"),
                   aes(x=xVal, y=value, fill=TestProbs, col=TestProbs, 
                       size=type)) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs != "A" & TestProbs != "D" & name=="Cases" & TestDay=="8")) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs=="A" & TestDay=="8"), alpha=.75) +
  geom_col(data=Test_Data_A %>% filter(sampled_time > 30 & sampled_time < 100 & TestProbs=="D" & TestDay=="8"), alpha=.75) +
  geom_vline(data=sample_time_dat %>% filter(label_factor=="Declining Epidemic"),aes(xintercept=samp_times,linetype=lty)) +
  facet_grid(rows=vars(type), cols=vars(label_factor),
             scales="free_x") +
  scale_y_continuous(name="")+ #,breaks=seq(0,1500,by=500),expand=c(0,0)) +
  scale_x_continuous(name="Test day",breaks=seq(32, 88,by=7)) +
  scale_linetype_manual(values=c("dashed","dotted"),
                        labels=c("Analysis 1: Epidemic Day 60",
                                 "Analysis 2: Epidemic Day 88")) +
  scale_size_manual(breaks=c("Reported Case Counts","Surveillance Ct Samples"),
                    values=c(0.05,0.2)) +
  scale_fill_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                    values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Flat Testing","Rising Testing","Falling Testing",
                             "One Sample Period: 3000 Tests","Two Sample Periods: 1500 Tests Each",
                             "Three Sample Periods: 1000 Tests Each",
                             "Three Sample Periods: 500, 1000, 1500 Tests","Three Sample Periods: 1500, 1000, 500 Tests")) +
  scale_color_manual(breaks=c("A","B","D","R1","R2","R3","R4","R5"),
                     values=unname(AAAS_palette[c("grey1","purple2","purple1",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Flat Testing","Rising Testing","Falling Testing",
                              "One Sample Period: 3000 Tests","Two Sample Periods: 1500 Tests Each",
                              "Three Sample Periods: 1000 Tests Each",
                              "Three Sample Periods: 500, 1000, 1500 Tests","Three Sample Periods: 1500, 1000, 500 Tests")) +
  guides(size="none") +
  export_theme +
  theme(strip.text.y = element_text(size=6))
p.A.right
p.A.left + p.A.right + plot_layout(nrow=1,ncol=2, guides="collect")

## Panel B original: Plot of estimated Rt values, median +/- 1 SD across the simulations
Res_Data_B <- Res_Summ %>% filter(TestProbs %in% c("A","B","D","T") | (size==3000 & analysis=="Ct"))
p.B.meds <- ggplot(Res_Data_B,
              aes(x=xVal, y=median, ymin=Q25, ymax=Q75, 
                  color=TestProbs, fill=TestProbs, shape=TestProbs, linetype=TestProbs,
                  size=TestProbs)) + 
  geom_hline(yintercept=1,col="black",size=0.25) +
  geom_hline(data=Res_Data_B %>% filter(type=="Simulation Truth"),aes(yintercept=median,linetype="T"),
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
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_size_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5","L"), 
                    values=c(1,1,1,2,1,1,1,1,1,0.5)/1.2,
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing",NA)) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing")) +
  guides(size="none") +
  facet_wrap(~label_factor, scales="free_y") +
  labs(tag="B")
p.B.meds

## Panel B used: Plot of estimated Rt values from each simulation
Res_Data_B_2 <- Res_last2 %>% filter(TestProbs %in% c("A","B","D","T") | (size==3000 & analysis=="Ct"))
p.B.dots <- ggplot(data=Res_Data_B_2, aes(x=xVal, y=median, 
                                   color=TestProbs, fill=TestProbs,shape=TestProbs)) + 
  geom_hline(yintercept=1,col="black",size=0.25) +
  geom_hline(data=Res_Data_B %>% filter(type=="Simulation Truth"),aes(yintercept=median),
             linetype="dashed",col="black",size=0.25) +
  geom_point(data=Res_Data_B %>% filter(type=="Simulation Truth"),
             aes(x=xVal,y=median,shape=TestProbs,color=TestProbs,fill=TestProbs),
             size=1.5) +
  geom_jitter(height=0, width=.2, size=.5) +
  export_theme+
  facet_wrap(~label_factor, scales="free_y") +
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"),
                     values=c(16,16,16,8,16,16,16,16,16),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_x_discrete(name="Testing scheme and estimation method",
                      breaks=NULL, labels=NULL) +
  scale_y_continuous(name=expression("Median estimated "*R[t]))+
  guides(linetype="none") +
  theme(panel.grid.major.y=element_blank(),
        plot.tag=element_text(face="bold"))+
  labs(tag="B")
p.B.dots

ggsave(filename=paste0(plot_out,"/Figure3.pdf"),
       plot=(p.A.left+p.A.right+plot_layout(nrow=1,ncol=2,guides="collect"))/p.B.dots,
       width=7.5, height=6, units="in", dpi=300)

## Panel B alternate: Plot of estimated Rt values from each simulation

#for (Pos_type in c("Pos_SEIR","Pos_SEEIRR")) {
for (Pos_type in c("Pos_SEIR")) {
Res_Data_B_3 <- Res_last2 %>% filter(TestProbs %in% c("A","B","D","T") | (size==3000)) %>% 
  filter(analysis %in% c("Ct",Pos_type,"Cases"))
p.B.dots.alt <- ggplot(data=Res_Data_B_3, aes(x=xVal, y=median, 
                                          color=TestProbs, fill=TestProbs,shape=TestProbs,
                                          alpha=if_else(analysis==Pos_type,"Positivity","Other"))) + 
  geom_hline(yintercept=1,col="black",size=0.25) +
  geom_hline(data=Res_Data_B %>% filter(type=="Simulation Truth"),aes(yintercept=median),
             linetype="dashed",col="black",size=0.25) +
  geom_point(data=Res_Data_B %>% filter(type=="Simulation Truth") %>% mutate(analysis="Truth"),
             # aes(x=xVal,y=median,
             #     shape=TestProbs,color=TestProbs,alpha="Other"),
             inherit.aes=TRUE,
             size=1.5) +
  geom_jitter(height=0, width=.1, size=.5) +
  export_theme+
  facet_wrap(~label_factor, scales="free_y") +
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"),
                     values=c(16,16,16,8,16,16,16,16,16),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Positivity","Other"),
                     values=c(0.4,1),
                     labels=c("Surveillance Sample: Positivity Only",""), guide=NULL) +
  scale_x_discrete(name="Testing scheme and estimation method", breaks=NULL, labels=NULL) +
  scale_y_continuous(name=expression("Median estimated "*R[t])) +
  # geom_text(data=tibble(x=c(2,5,2,5)-.2,y=c(1.7,2.7,0.6,0.465),
  #                       label_factor=factor(c("Growing Epidemic","Growing Epidemic",
  #                                      "Declining Epidemic","Declining Epidemic"),
  #                                      levels=c("Growing Epidemic","Declining Epidemic")),
  #                       text_val=c("{","}","}","}")),
  #           aes(x=x,y=y,label=text_val), inherit.aes=FALSE,
  #           angle=90, size=14, vjust=0.5) +
  geom_text(data=tibble(x=c(2,5,2,5),y=c(1.35,3.65,0.7,0.63),
                        label_factor=factor(c("Growing Epidemic","Growing Epidemic",
                                              "Declining Epidemic","Declining Epidemic"),
                                            levels=c("Growing Epidemic","Declining Epidemic")),
                        text_val=c("Ct Values","Positivity\nOnly",
                                   "Ct Values","Positivity\nOnly")),
            aes(x=x,y=y,label=text_val), inherit.aes=FALSE,
            hjust=0.5, vjust=1, size=8/ggplot2:::.pt) +
  guides(linetype="none") +
  theme(panel.grid.major.y=element_blank(),
        plot.tag=element_text(face="bold"))+
  labs(tag="B")
p.B.dots.alt

ggsave(filename=paste0(plot_out,"/Figure3_wPos_",substr(Pos_type, 5, nchar(Pos_type)),".pdf"),
       plot=(p.A.left+p.A.right+plot_layout(nrow=1,ncol=2,guides="collect"))/p.B.dots.alt,
       width=7.5, height=6, units="in", dpi=300)

ggsave(filename=paste0(plot_out,"/Figure3_wPos_",substr(Pos_type, 5, nchar(Pos_type)),".png"),
       plot=(p.A.left+p.A.right+plot_layout(nrow=1,ncol=2,guides="collect"))/p.B.dots.alt,
       width=7.5, height=6, units="in", dpi=300)

### Comparison of Sample Size and Ct vs. Positivity Only:
Res_Data_C <- Res_Summ %>% filter(analysis %in% c("Ct",Pos_type)) %>%
  mutate(size_pos = if_else(size>100,10,0) + if_else(size>200,10,0) + if_else(size>500,10,0) + 
           if_else(size>1000,10,0) + if_else(size>2000,10,0),
         analysis=if_else(analysis==Pos_type,"Positivity",analysis)) %>% 
  mutate(xVal2=xVal + size_pos,
         SampleLabel=factor(TestProbs,
                            levels=c("R1","R2","R3","R4","R5"),
                            labels=c("One","Two","Three, Flat Testing",
                                     "Three, Rising Testing","Three, Falling Testing")))
p.C.meds <- ggplot(Res_Data_C,
                   aes(x=size_pos+ifelse(analysis=="Positivity",2,-2), y=median, ymin=Q25, ymax=Q75, 
                       color=TestProbs, fill=TestProbs, shape=TestProbs,
                       alpha=analysis)) + 
  geom_hline(yintercept=1,col="black") +
  geom_hline(data=Res_Summ %>% filter(type=="Simulation Truth") %>% ungroup() %>% select(t,label_factor,median) %>%
               full_join(Res_Data_C %>% ungroup() %>% select(label_factor,TestProbs) %>% distinct()),
             mapping=aes(yintercept=median,linetype="T"),
             col="grey40") +
  geom_point() + 
  geom_errorbar(width=0.2) + 
  export_theme + 
  labs(title="Sample Periods")+
  theme(panel.border=element_rect(fill=NA),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black"),
        legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=-45,hjust=0),
        plot.tag=element_text(face="bold")) +
  scale_x_continuous(name="Sample Size",
                     breaks=seq(0,50,by=10),
                     labels=as.character(unname(test_probs))) +
  scale_y_continuous(name=expression("Median and IQR of Median Posterior "*R[t]*" Estimates"))+
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing"),
                    guide=FALSE) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Ct","Positivity"), values=c(1, 0.5),
                     labels=c("Ct Values","Positivity Only")) +
  guides(size="none") +
  facet_grid(rows=vars(label_factor), cols=vars(SampleLabel), scales="free_y")
p.C.meds

p.D.CIW <- ggplot(Res_Data_C,
                   aes(x=size_pos+ifelse(analysis=="Positivity",2,-2), y=MeanWidth, 
                       color=TestProbs, fill=TestProbs, shape=TestProbs,
                       alpha=analysis)) + 
  geom_point() +
  export_theme + 
  
  labs(title="Sample Periods")+
  theme(panel.border=element_rect(fill=NA),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black"),
        legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=-45,hjust=0),
        plot.tag=element_text(face="bold")) +
  scale_x_continuous(name="Sample Size",
                     breaks=seq(0,50,by=10),
                     labels=as.character(unname(test_probs))) +
  scale_y_continuous(name=expression("Mean Width of 95% Credible Intervals for "*R[t]))+
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing"),
                    guide=FALSE) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Ct","Positivity"), values=c(1, 0.5),
                     labels=c("Ct Values","Positivity Only")) +
  guides(size="none") +
  facet_grid(rows=vars(label_factor), cols=vars(SampleLabel), scales="free_y")
p.D.CIW

p.E.MSE <- ggplot(Res_Data_C,
                  aes(x=size_pos+ifelse(analysis=="Positivity",2,-2), y=MSE, 
                      color=TestProbs, fill=TestProbs, shape=TestProbs,
                      alpha=analysis)) + 
  geom_point() +
  export_theme + 
  
  labs(title="Sample Periods")+
  theme(panel.border=element_rect(fill=NA),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black"),
        legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=-45,hjust=0),
        plot.tag=element_text(face="bold")) +
  scale_x_continuous(name="Sample Size",
                     breaks=seq(0,50,by=10),
                     labels=as.character(unname(test_probs))) +
  scale_y_continuous(name=expression("Mean Squared Error of Median Posterior "*R[t]*" Estimates"))+
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing"),
                    guide=FALSE) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Ct","Positivity"), values=c(1, 0.5),
                     labels=c("Ct Values","Positivity Only")) +
  guides(size="none") +
  facet_grid(rows=vars(label_factor), cols=vars(SampleLabel), scales="free_y")
p.E.MSE

# ggsave(filename=paste0(plot_out,"/Figure3_SScomp_",substr(Pos_type, 5, nchar(Pos_type)),".pdf"),
#        plot=p.C.meds+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
#          p.D.CIW+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position="none")+
#          p.E.MSE+theme(legend.position="none")+
#          plot_layout(ncol=1,nrow=3,guides="collect"),
#        width=7.6, height=9, units="in", dpi=300)

p.F.EstSign <- ggplot(Res_Data_C,
                  aes(x=size_pos+ifelse(analysis=="Positivity",2,-2), y=EstimateSign*100, 
                      color=TestProbs, fill=TestProbs, shape=TestProbs,
                      alpha=analysis)) + 
  geom_point() +
  export_theme + 
  labs(title="Sample Periods")+
  theme(panel.border=element_rect(fill=NA),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black"),
        legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=-45,hjust=0),
        plot.tag=element_text(face="bold")) +
  scale_x_continuous(name="Sample Size",
                     breaks=seq(0,50,by=10),
                     labels=as.character(unname(test_probs))) +
  scale_y_continuous(name=expression("Median Posterior "*R[t]*" Estimates with Correct Direction (%)"),
                     limits=c(0,100))+
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing"),
                    guide=FALSE) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Ct","Positivity"), values=c(1, 0.5),
                     labels=c("Ct Values","Positivity Only")) +
  guides(size="none") +
  facet_grid(rows=vars(label_factor), cols=vars(SampleLabel), scales="fixed")
p.F.EstSign

p.G.CrISign <- ggplot(Res_Data_C,
                      aes(x=size_pos+ifelse(analysis=="Positivity",2,-2), y=CrIPower*100, 
                          color=TestProbs, fill=TestProbs, shape=TestProbs,
                          alpha=analysis)) + 
  geom_point() +
  export_theme + 
  labs(title="Sample Periods")+
  theme(panel.border=element_rect(fill=NA),
        panel.grid.major=element_blank(),
        axis.line.y=element_line(color="black"),
        legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=-45,hjust=0),
        plot.tag=element_text(face="bold")) +
  scale_x_continuous(name="Sample Size",
                     breaks=seq(0,50,by=10),
                     labels=as.character(unname(test_probs))) +
  scale_y_continuous(name=expression(R[t]*" 95% Credible Intervals Entirely in Correct Direction (%)"),
                     limits=c(0,100))+
  scale_shape_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=c(22,24,25,8,22,22,22,24,25),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_linetype_manual(name="Testing Scheme and Estimation Method",
                        breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                        values=c(1,1,1,2,1,1,1,1,1),
                        labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                                 "Reported Case Counts, Falling Testing","Simulation Truth",
                                 "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                                 "Three Surveillance Ct Sample Periods, Flat Testing",
                                 "Three Surveillance Ct Sample Periods, Rising Testing",
                                 "Three Surveillance Ct Sample Periods, Falling Testing")) +
  scale_color_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                     values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                  "red1","red2","teal1","blue1","green1")]),
                     labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                              "Reported Case Counts, Falling Testing","Simulation Truth",
                              "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                              "Three Surveillance Ct Sample Periods, Flat Testing",
                              "Three Surveillance Ct Sample Periods, Rising Testing",
                              "Three Surveillance Ct Sample Periods, Falling Testing"),
                     guide=FALSE) +
  scale_fill_manual(name="Testing Scheme and Estimation Method",
                    breaks=c("A","B","D","T","R1","R2","R3","R4","R5"), 
                    values=unname(AAAS_palette[c("grey1","purple2","purple1","black",
                                                 "red1","red2","teal1","blue1","green1")]),
                    labels=c("Reported Case Counts, Flat Testing","Reported Case Counts, Rising Testing",
                             "Reported Case Counts, Falling Testing","Simulation Truth",
                             "One Surveillance Ct Sample Period", "Two Surveillance Ct Sample Periods, Flat Testing", 
                             "Three Surveillance Ct Sample Periods, Flat Testing",
                             "Three Surveillance Ct Sample Periods, Rising Testing",
                             "Three Surveillance Ct Sample Periods, Falling Testing"),
                    guide=FALSE) +
  scale_alpha_manual(name="Testing Scheme and Estimation Method",
                     breaks=c("Ct","Positivity"), values=c(1, 0.5),
                     labels=c("Ct Values","Positivity Only")) +
  guides(size="none") +
  facet_grid(rows=vars(label_factor), cols=vars(SampleLabel), scales="fixed")
p.G.CrISign

# ggsave(filename=paste0(plot_out,"/Figure3_SScomp2_",substr(Pos_type, 5, nchar(Pos_type)),".pdf"),
#        plot=p.F.EstSign+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
#          p.G.CrISign+theme(legend.position="none")+
#          plot_layout(ncol=1,nrow=2,guides="collect"),
#        width=7.6, height=6.5, dpi=300)

ggsave(filename=paste0(plot_out,"/Figure3_SScompFull_",substr(Pos_type, 5, nchar(Pos_type)),".pdf"),
       plot=p.C.meds+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                           axis.title.y=element_text(size=6),
                           strip.text.y=element_text(size=6))+
         labs(tag="A")+
         p.E.MSE+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       strip.text.y=element_text(size=6),
                       strip.text.x=element_blank(),legend.position="none")+
         labs(tag="B",title=NULL)+
         p.D.CIW+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       strip.text.y=element_text(size=6),
                       strip.text.x=element_blank(),legend.position="none")+
         labs(tag="C",title=NULL)+
         p.G.CrISign+theme(legend.position="none",strip.text.y=element_text(size=6),
                           strip.text.x=element_blank(),
                           axis.title.y=element_text(size=6))+labs(tag="D",title=NULL)+
         plot_layout(ncol=1,nrow=4,guides="collect"),
       width=7.6, height=10, dpi=300)
ggsave(filename=paste0(plot_out,"/Figure3_SScompFull_",substr(Pos_type, 5, nchar(Pos_type)),".png"),
       plot=p.C.meds+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                           axis.title.y=element_text(size=6),
                           strip.text.y=element_text(size=6))+
         labs(tag="A")+
         p.E.MSE+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       strip.text.y=element_text(size=6),
                       strip.text.x=element_blank(),legend.position="none")+
         labs(tag="B",title=NULL)+
         p.D.CIW+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                       axis.title.y=element_text(size=6),
                       strip.text.y=element_text(size=6),
                       strip.text.x=element_blank(),legend.position="none")+
         labs(tag="C",title=NULL)+
         p.G.CrISign+theme(legend.position="none",strip.text.y=element_text(size=6),
                           strip.text.x=element_blank(),
                           axis.title.y=element_text(size=6))+labs(tag="D",title=NULL)+
         plot_layout(ncol=1,nrow=4,guides="collect"),
       width=7.6, height=10, dpi=300)
}

