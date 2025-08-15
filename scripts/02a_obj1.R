
# This script fits a hierarchical logistic regression model of the effects
# of treatment on odds of attrition
# Random effects for treatment, placebo as global reference

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import tidied data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Inspect glp1 trials
a_insp <- a_imp_df %>% 
  group_by(trial_id) %>% 
  filter(any(class == "glp1"))

# Set placebo as global reference treatment
# This brings the unadjusted estimates in line with mean relative effects
b_set_ref <- a_imp_df %>% 
  mutate(
    class = factor(class),
    class = relevel(class, ref = "placebo"),
    arm_compl = arm_n - arm_attr,
    across(arm_n:arm_compl, as.integer)
  )

# Fit hierarchical logistic regression via brms
# Random effects for treatment class across trials
b_fit <- brm(
  formula = arm_attr | trials(arm_n) ~ class + (class | trial_id),
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

# Posterior predictive checks --------------------------------------------------

# Clear memory after fitting model
gc()

# Import (or run continuously using model above)
b_fit <- readRDS("output/obj1/brm_fit.rds")

# Extract posterior draws and summarise
c_draws <- as_draws_df(b_fit)
c_draws_sum <- summarise_draws(c_draws)

# Convergence summary
c_draws_cnvg <- c_draws_sum %>% 
  filter(grepl("^b_", variable)) %>% 
  mutate(variable = gsub("^b_|^b_class_short", "", variable)) %>% 
  select(variable, rhat:ess_tail)

write_csv(c_draws_cnvg, "output/obj1/sum_cnvg.csv")

# Tidy draws
# Separate parameters and terms
c_draws_split <- c_draws %>% 
  select(b_Intercept, contains("b_class")) %>% 
  as_tibble() %>% 
  pivot_longer(-b_Intercept, names_to = "class", values_to = "draw") %>% 
  rename(intercept = b_Intercept) %>% 
  mutate(class = gsub("b_class", "", class))

# Summarise mean and quantiles for credible effect
c_draws_class_sum <- c_draws_split %>% 
  group_by(class) %>% 
  reframe(
    mean_log = mean(draw),
    median_log = median(draw),
    upper_log = quantile(draw, 0.975),
    lower_log = quantile(draw, 0.025),
    mean_or = mean(exp(draw)),
    median_or = median(exp(draw)),
    upper_or = quantile(exp(draw), 0.975),
    lower_or = quantile(exp(draw), 0.025)
  )

write_csv(c_draws_class_sum, "output/obj1/sum_res.csv")

# Plots ------------------------------------------------------------------------

# Predictive check on arm-level attrition counts
plot_ppc_trial <- pp_check(b_fit, type = "intervals", ndraws = 803)

ggsave(
  "output/obj1/ppc_arm_level.jpg",
  plot_ppc_trial,
  width = 6,
  height = 4,
  units = "in"
)

# Plot log-odds
plot_log <- c_draws_class_sum %>% 
  ggplot(aes(x = mean_log, xmin = lower_log, xmax = upper_log, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean log-odds (95% credible intervals)", y = "Treatment class")

plot_log

ggsave(
  "output/obj1/plot_log_odds.png",
  plot_log,
  width = 8,
  height = 4,
  units = "in"
)

# Plot odds ratio
plot_odds <- c_draws_class_sum %>% 
  ggplot(aes(x = mean_or, xmin = lower_or, xmax = upper_or, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

plot_odds

ggsave(
  "output/obj1/plot_odds.png",
  plot_odds,
  width = 8,
  height = 4,
  units = "in"
)
