
rm(list = ls())
gc()

options(scipen = 999)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
