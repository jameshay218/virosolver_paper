## Pre-script 2: Rt estimations based on MA case data
library(tidyverse)
library(patchwork)
library(EpiEstim)
library(zoo)
library(ggsci)
library(EpiNow2)

main_wd <- "~/Documents/GitHub/virosolver_paper/"
setwd(main_wd)

rerun <- FALSE

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

## NYT data
nyt_dat <- read_csv("data/us-counties.csv")
nyt_dat <- nyt_dat %>% 
  filter(state=="Massachusetts",
         county == "Suffolk") %>%
  group_by(state) %>%
  mutate(new_cases=cases-lag(cases, 1)) %>%
  filter(new_cases >= 0)

## Plot case counts
p1 <- nyt_dat %>% 
  ## Highlight anomalous cases
  #mutate(new_cases1=ifelse(new_cases > 3500, lag(new_cases,1), new_cases)) %>%
  ## Get rolling mean
  mutate(roll_mean=rollmean(new_cases, 7,fill=NA,align="right")) %>%
  ggplot()+ 
  geom_bar(aes(x=date,y=new_cases),stat="identity",alpha=0.25) +
  #geom_ribbon(aes(x=date,ymax=roll_mean,ymin=0),fill="grey70",alpha=0.5,col="grey20") +
  #scale_fill_npg() +
  scale_fill_manual(values=c("grey40","darkred")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,1000)) +
  scale_x_date(limits=as.Date(c("2020-04-01", "2020-12-01"), "%Y-%m-%d"), breaks="7 days",
               expand=c(0,0)) +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("Date") +
  ylab("New infections") +
  labs(tag="A")
p1

## Rt estimation on these case counts
rt_dat <- nyt_dat %>% 
  ungroup() %>%
  dplyr::select(date, new_cases) %>%
  drop_na() %>%
  rename(I=new_cases,
         dates=date)

reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(3), 0.5))
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 30

generation_time <- list(mean = EpiNow2::covid_generation_times[1, ]$mean,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = EpiNow2::covid_generation_times[1, ]$sd,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)

incubation_period <- list(mean = EpiNow2::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow2::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)
reported_cases <- EpiNow2::example_confirmed[1:50]
reported_cases <- rt_dat %>% rename(confirm=I, date=dates)

if(rerun){
  estimates <- EpiNow2::epinow(reported_cases = reported_cases, 
                               generation_time = generation_time,
                               delays = list(incubation_period, reporting_delay),
                               horizon = 7, samples = 4000, warmup = 1000, 
                               cores = 4, chains = 4, verbose = TRUE, 
                               adapt_delta = 0.95)
  saveRDS(estimates,"results/ma_county_rt_fit.RData")
}
  