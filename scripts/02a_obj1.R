
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Rstan config
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Fit hierarchical logistic regression via brms
# Random effects for treatment class across trials
b_fit <- brm(
  formula = arm_attr | trials(arm_n) ~ class_short + (1 + class_short | trial_id),
  data = a_imp_df,
  family = binomial(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "rstan",
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# Save (to prevent having to rerun each time)
b_fit <- readRDS("processed_data/brm_fit.rds")
saveRDS(b_fit, "processed_data/brm_fit.rds")