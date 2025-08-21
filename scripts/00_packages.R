
library(tidyverse)
library(broom)
library(brms)
library(posterior)
library(bayesplot)
library(multinma)
library(devtools)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Get set_agd_regression from development arm of multinma
# remotes::install_local("C:/multinma-feature-set_agd_regression")

