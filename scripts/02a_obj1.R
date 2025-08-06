
# This script is run within a high performance computing environment
# The slurm file used to run the script has the same file name
# Fits a hierarchical logistic regression model for main effects
# Random effects for treatment, placebo as global reference

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import tidied data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Set placebo as global reference treatment
b_set_ref <- a_imp_df %>% 
  mutate(
    class = factor(class),
    class = relevel(class, ref = "placebo"),
    arm_compl = arm_n - arm_attr,
    across(arm_n:arm_compl, as.integer)
  )

# Remove IPD trial with 2 attrition events to see if it resolves
# divergent transition
b_ipd_rm <- b_set_ref %>% filter(source == "ipd", trial_id != "NCT01769378")

# Fit hierarchical logistic regression via brms
# Random effects for treatment class across trials
b_fit <- brm(
  # formula = arm_attr | trials(arm_n) ~ class + (1 + class | trial_id),
  formula = arm_attr | trials(arm_n) ~ class + (1 | trial_id) + (0 + class | trial_id),
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

saveRDS(b_fit, "output/obj1/brm_fit.rds")

# Fit with just the IPD trials
b_fit_ipd <- brm(
  formula = arm_attr | trials(arm_n) ~ class + (1 | trial_id) + (0 + class | trial_id),
  data = b_ipd_rm,
  family = binomial(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "rstan",
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

saveRDS(b_fit_ipd, "output/obj1/brm_fit_ipd.rds")
