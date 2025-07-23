
# This script is run within a high performance computing environment
# The associated slurm file used to run the script is has the same file name
source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import tidied data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Set placebo as global reference treatment
b_set_ref <- a_imp_df %>% 
  mutate(
    class = factor(class),
    class = relevel(class, ref = "placebo")
  )

# Fit hierarchical logistic regression via brms
# Random effects for treatment class across trials
b_fit <- brm(
  formula = arm_attr | trials(arm_n) ~ class + (1 + class | trial_id),
  data = b_set_ref,
  family = binomial(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "rstan",
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

b_fit <- readRDS("outputs/brm_fit.rds")