## James Hay 31/12/2020
##  - This script runs the analyses to investigate transmission dynamics in 4 Massachusetts care homes
##  - First, we look at the prevalence of positive samples at 3 sampling rounds
##  - Next, we fit a modified SEIR model (SEEIRR) to these point prevalence estimates
##  - Assuming that these fits demonstrate "ground truth", we then fit the Ct model to each individual cross section
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(data.table)
library(patchwork)
library(fitdistrplus)
library(deSolve)
library(lazymcmc) ## devtools::install_github("jameshay218/lazymcmc")
library(doParallel)
## Where all Git repos are saved
GIT_WD <- "~/Documents/Github"
devtools::load_all(paste0(GIT_WD,"/virosolver"))

## CHANGE TO MAIN WD
## Important to set this to the full file path, as on L205 the foreach loop
## must move to the correct working directory to source the model functions
main_wd <- "~/Documents/GitHub/virosolver_paper/"
chainwd <- "~/Documents/GitHub/virosolver_paper/mcmc_chains/1.nursing_home" ## Where to save MCMC chains
chainwd2 <- "~/Documents/GitHub/virosolver_paper/mcmc_chains/1.nursing_home_ct/"
plot_wd <- "~/Documents/GitHub/virosolver_paper/plots/1.nursing_home_ct"
setwd(main_wd)

## IMPORTANT - change this flag to TRUE if running the MCMC for the first time
rerun_mcmc_seeirr <- FALSE
rerun_mcmc_ct <- TRUE

## Clean data and set up model parameters. Also some data plots.
source("1.nursing_homes/1.nursing_home_analyses_headers.R")
## Fit the SEEIRR compartmental model to point prevalence
source("1.nursing_homes/1.nursing_home_analyses_SEEIRR.R")
## Fit the various Ct models to the nursing home cross sections
source("1.nursing_homes/1.nursing_home_analyses_ct_model.R")
## Generate Figure 2
source("1.nursing_homes/Figure2.R")
## Generate corresponding supplementary figures
source("1.nursing_homes/SuppFigure_fits.R")



