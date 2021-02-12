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

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- pal_aaas(palette = c("default"), alpha = 1)(10)

#setwd("/Users/leekennedy-shaffer/Documents/Harvard/Research/Infectious Diseases/COVID-19/PCR Cts/gitversion/ct_inference")
setwd("~/Documents/GitHub/virosolver_paper")

main_theme <- theme_classic() +
  theme(axis.ticks.length.x = unit(0.2,"cm"),
        text=element_text(family="sans"),
        plot.margin = margin(-1,-1,-1,-1,unit="cm"),
    panel.grid.minor.x=element_blank())
source("code/plot_funcs.R")
source("code/linelist_sim_funcs.R")

R0 <- 2
#### Viral Load Kinetics Parameters: ####
vl_pars <- read.csv("pars/massachusetts/partab_fitted_bwh.csv")
pars <- vl_pars$values
names(pars) <- vl_pars$names


source("scripts/Figure1/generate_sim_data.R")

## Simulation plot
source("scripts/Figure1/seir_fig.R")
## TSI median against growth rate
source("scripts/Figure1/tsi_gr_relationship.R")
## Individual trajectories
source("scripts/Figure1/vl_trajectories.R")
## TSI over time
source("scripts/Figure1/tsi_over_time.R")
## VL model
source("scripts/Figure1/vl_model.R")
## 
source("scripts/Figure1/cts_over_time.R")
## Ct growth rate relationship
source("scripts/Figure1/ct_gr_relationship.R")
## Simulation Rt relationships
source("scripts/Figure1/simulation_rt_relationships.R")


p_lhs <- (p_seir+p_traj+p_tsi_time+p_ct_time) + plot_layout(ncol=1)
p_rhs <- (p_tsi_relationship+p_vl+p_ct+p3)+ plot_layout(ncol=1)

main_p1 <- ((p_lhs) | (p_rhs)) + plot_layout(widths=c(2.5,1),ncol=2)

main_p1


ggsave(filename="figures/fig1_jh3.png",
       plot=main_p1,
       width=10, height=10, units="in", dpi=600)
ggsave(filename="figures/fig1_jh3.pdf",
       plot=main_p1,
       width=10, height=10, units="in", dpi=600)
