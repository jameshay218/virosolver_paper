require(tidyverse)
require(ggplot2)
require(ggthemes)
require(ggsci)
require(ggpubr)
require(patchwork)
require(triangle)
require(ggrepel)
require(ggpubr)
require(extraDistr)
require(mvtnorm)
library(moments)

devtools::load_all("~/Documents/GitHub/virosolver")
devtools::document("~/Documents/GitHub/virosolver")



viral_load_pars <- c(beta = 0.1, tshift = 0, desired_mode = 3, viral_peak = 19.7, 
                     obs_sd = 5, sd_mod = 0.789, sd_mod_wane = 14, true_0 = 40, 
                     intercept = 40, LOD = 3, incu = 5, t_switch = 13.3, level_switch = 38, 
                     wane_rate2 = 1000, prob_detect = 0.103, t_unit = 1)

seir_pars["I0"] <- 1e-6
sir_pars <- c(R0 = 5, infectious = 5, I0 = 1e-06, t0 = 0)
times <- seq(0,500,by=1)
all_dat <- NULL
R0s <- seq(1,5,by=0.2)
GRs <- round(seq(-0.3,0.3,by=0.025),3)
res <- calculate_gr_ct_relationships_R0s(sir_pars, viral_load_pars, times, R0s=R0s)
virosolver::compare_tsi_dist_at_gr(res[[6]], GRs, 0.01)

