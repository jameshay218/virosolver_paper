## Pre-script 2: Rt estimations based on MA case data
library(tidyverse)
library(patchwork)
library(EpiEstim)
library(zoo)
library(ggsci)
library(EpiNow2)

main_wd <- "~/Documents/GitHub/virosolver_paper/"
setwd(main_wd)

rerun <- TRUE

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

nyt_dat <- read_csv("~/Documents/local_data/BWH/covid-19-raw-data-11-2-2021.csv")  %>% 
  filter(County == "Suffolk") %>% 
  mutate(date=lubridate::mdy(Date)) %>% 
  rename(new_cases=`New Confirmed Cases`)

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
  #scale_x_date(limits=as.Date(c("2020-04-01", "2020-12-01"), "%Y-%m-%d"), breaks="7 days",
  #             expand=c(0,0)) +
  theme_classic()+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        #axis.text.x=element_text(angle=45,hjust=1),
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank(),
        #axis.ticks.x=element_blank()
        ) +
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
         dates=date) %>%
  filter(I >= 0,
         dates >= "2021-06-01")

reporting_delay <- list(mean = convert_to_logmean(4, 1),
                        mean_sd = 0.1,
                        sd = convert_to_logsd(4, 1),
                        sd_sd = 0.1,
                        max = 15)


generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

reported_cases <- EpiNow2::example_confirmed[1:50]
reported_cases <- rt_dat %>% rename(confirm=I, date=dates)

if(rerun){
  estimates <- EpiNow2::epinow(reported_cases = reported_cases, 
                               generation_time = generation_time,
                               delays = delay_opts(incubation_period, reporting_delay), 
                               rt = rt_opts(prior = list(mean = 2, sd = 0.2)),
                               gp = gp_opts(basis_prop = 0.2),
                               horizon = 14, 
                               verbose = TRUE,
                               return_output=TRUE,
                               stan = stan_opts())
  saveRDS(estimates,"results/ma_county_rt_fit.RData")
}
  