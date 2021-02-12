## Code for plotting
source("code/plot_funcs.R")
export_theme <- export_theme + theme(axis.text.x=element_text(size=7),
                                     axis.text.y=element_text(size=7),
                                     axis.title.x=element_text(size=8),
                                     axis.title.y=element_text(size=8))

## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(chainwd2)) dir.create(chainwd2,recursive = TRUE)


AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="1B1919FF")


## Which incidence model to use?
ct_versions <- c("exp","exp","seir","seir","seir","gp","gp")
## Use all Cts or just positive?
use_pos <- c(FALSE,TRUE,FALSE,TRUE,
             FALSE,FALSE,FALSE)
## Consider cross-sections independently or together?
cumu_data <- c(FALSE,FALSE,FALSE,FALSE,
               TRUE,FALSE,TRUE)
## Parameter table
par_table_files <- c("pars/nursing_homes/partab_exp_model.csv","pars/nursing_homes/partab_exp_pos_model.csv",
                     "pars/nursing_homes/partab_seir_model.csv","pars/nursing_homes/partab_seir_model.csv","pars/nursing_homes/partab_seir_model.csv",
                     "pars/nursing_homes/partab_gp_model.csv","pars/nursing_homes/partab_gp_model.csv")
## Prior functions
prior_funcs <- c(prior_func_hinge_exp, prior_func_hinge_exp,
                 prior_func_hinge_seir, prior_func_hinge_seir, prior_func_hinge_seir,
                 prior_func_hinge_gp, prior_func_hinge_gp)
## Incidence functions
inc_funcs <- c(exponential_growth_model, exponential_growth_model,
               solveSEIRModel_rlsoda_wrapper, solveSEIRModel_rlsoda_wrapper, solveSEIRModel_rlsoda_wrapper,
               gaussian_process_model, gaussian_process_model)
## Using cap on infection ages or to the start?
age_maxs <- c(35, 35,NA, NA,
              NA, 35, NA)

## Put together
analysis_control <- tibble(version=ct_versions, use_pos=use_pos, cumu_data=cumu_data,
                           par_table_file=par_table_files, prior_func=prior_funcs,
                           inc_func=inc_funcs, age_max=age_maxs)
analysis_control <- analysis_control[c(1,3),]

## Manage MCMC runs and parallel runs
nchains <- 3
n_clusters <- 8
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

## MCMC parameters for SEEIRR model fit
## Change these MCMC control parameters to run longer/shorter chains
mcmcPars1 <- c("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=10000,"save_block"=1000)
mcmcPars2 <- c("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
               "thin"=10,"adaptive_period"=30000,"save_block"=1000)

## MCMC parameters for Ct model fits if using parallel tempering branch
n_temperatures <- 10
mcmcPars_ct_pt <- list("iterations"=50000,"popt"=0.44,"opt_freq"=1000,
                       "thin"=1,"adaptive_period"=50000,"save_block"=1000,
                       "temperature" = seq(1,101,length.out=n_temperatures),
                       "parallel_tempering_iter" = 5,"max_adaptive_period" = 50000, 
                       "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)


## Duration of epidemic in days
times <- 0:365
t<-times
## Date to consider day 0
min_date <- as.Date("2020-03-01")

## Read in nursing home data Cts
## Note that data are grouped by week of collection
dat_wide <- read_csv(file = "data/nursing_home_data.csv")
dat_wide$collection_date <- as.Date(dat_wide$collection_date,format="%m/%d/%y")

dat_wide <- dat_wide %>% left_join(
  dat_wide %>% filter(is_nursinghome==1) %>% 
  group_by(location, week) %>% 
  count(collection_date) %>% 
  filter(n == max(n)) %>% 
  rename(mean_week_date=collection_date) %>% 
  ungroup() %>% 
  dplyr::select(-n))

## 3rd sampling time for nursing home 3 is not random collection
dat_wide <- dat_wide %>% filter(!(location == "NH 4" & mean_week_date == "2020-04-29"))

## Look at distribution properties
dat_properties <- dat_wide %>% filter(is_nursinghome==1) %>% group_by(location, mean_week_date) %>%
  filter(!is.na(N2) & N2 < 40) %>%
  summarize(ct_median=median(N2,na.rm=TRUE),
            ct_skew=moments::skewness(N2,na.rm=TRUE)) %>%
  rename(Location=location)

p_median <- dat_properties %>% 
  ggplot() +
  geom_point(aes(x=mean_week_date,y=ct_median,group=Location,col=Location),size=2) +
  geom_line(aes(x=mean_week_date,y=ct_median,group=Location,col=Location),size=1) +
  scale_color_viridis_d() +
  scale_y_continuous(trans="reverse") +
  export_theme  +
  ylab("Median of N2 Cts") +
  xlab("Swab collection date") +
  labs(tag="D") +
  theme(legend.position=c(0.9,0.9))
p_skew <- dat_properties %>% 
  ggplot() +
  geom_point(aes(x=mean_week_date,y=ct_skew,group=Location,col=Location),size=2) +
  geom_line(aes(x=mean_week_date,y=ct_skew,group=Location,col=Location),size=1) +
  scale_color_viridis_d() +
  export_theme  +
  ylab("Skewness of N2 Cts") +
  xlab("Swab collection date") +
  labs(tag="E") +
  theme(legend.position="none")

p_nh_relationships <- (p_median/p_skew)

## Look at sampling times
p_samples <- dat_wide %>% filter(is_nursinghome==1) %>%
  mutate(mean_week_date1 = as.factor(mean_week_date)) %>%
  rename(Week=mean_week_date1) %>%
  ggplot() +
  geom_histogram(aes(x=collection_date,fill=Week),binwidth=1) +
  geom_vline(aes(xintercept=mean_week_date),linetype="dashed",size=0.5,col="grey30") +
  facet_wrap(~location, nrow=1) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Swab collection date") +
  ylab("Count") +
  ggsci::scale_fill_aaas() +
  export_theme +
  theme(axis.text.x=element_text(hjust=1,angle=45)) +
  labs(tag="A")

## Look at Ct values
p_ct_N2 <- dat_wide %>% filter(is_nursinghome==1) %>% 
  ggplot() + 
  geom_histogram(aes(x=N2),binwidth=1,fill="grey70",col="grey40") +
  facet_wrap(~location,nrow=1) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Count") +
  xlab("Ct") +
  ggtitle("N2 Ct value") +
  export_theme +
  theme(plot.title=element_text(hjust=0)) +
  #geom_vline(xintercept=33,linetype="dashed", col="red") +
  labs(tag="B")


p_ct_N1 <- dat_wide %>% filter(is_nursinghome==1) %>% 
  ggplot() + 
  geom_histogram(aes(x=N1),binwidth=1,fill="grey70",col="grey40") +
  facet_wrap(~location,nrow=1) +
  ggtitle("N1 Ct value") +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Count") +
  xlab("Ct") +
  export_theme +
  theme(plot.title=element_text(hjust=0)) +
  #geom_vline(xintercept=33,linetype="dashed", col="red") +
  labs(tag="C")


p_lhs <- p_samples/p_ct_N2/p_ct_N1
fig_nh_data <- (p_lhs | p_nh_relationships ) + plot_layout(widths=c(2,1))

ggsave("figures/supplement/nh_data.pdf",p_lhs,height=6,width=8)
ggsave("figures/supplement/nh_data.png",p_lhs,height=6,width=8)

## Subset by only nursing homes
dat_use <- dat_wide %>% filter(is_nursinghome==1) %>%
  dplyr::select(location, mean_week_date, N2) %>%
  rename(ct=N2)
unique_data_combs <- dat_use %>% dplyr::select(location, mean_week_date) %>% distinct()

##########################################
# 1. Look at raw data for prevalence
##########################################
## Look at prevalence over time
prev_time <- dat_wide %>%
  filter(is_nursinghome == 1) %>%
  filter(result %in% c("POS","NEG")) %>%
  group_by(location, result, week) %>%
  tally() %>%
  pivot_wider(values_from=n, names_from=result,values_fill=list(n=0)) %>%
  group_by(location, week) %>%
  mutate(prev=POS/(POS+NEG),
         lower_confint=prop.test(POS, POS+NEG)$conf.int[1],
         upper_confint=prop.test(POS,POS+NEG)$conf.int[2])

p_prev <- prev_time %>% ggplot() +
  geom_point(aes(x=week,y=prev),size=0.25) +
  geom_errorbar(aes(x=week,ymin=lower_confint,ymax=upper_confint),size=0.25,width=0.1) +
  facet_wrap(~location,ncol=2) +
  ylab("Prevalence (95% binomial confidence intervals)") +
  xlab("Week") +
  ggtitle("Prevalence over time by location") +
  export_theme +
  theme(legend.position=c(0.8,0.1),
        axis.text.x=element_text(angle=0,hjust=0.5),
        panel.grid.major = element_line(color="grey70",size=0.1)) +
  scale_y_continuous(limits=c(0,1))
p_prev
