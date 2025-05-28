
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
# Commented out after saving model fit
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

# Posterior checks -------------------------------------------------------------

# Model summary
summary(b_fit)

# Trace plots
trace_draws <- as_draws_df(b_fit)
plot_trace <- mcmc_trace(trace_draws, pars = vars(starts_with("b_")))

ggsave(
  "output/obj1/trace_treatment_effects.jpg",
  plot_trace,
  width = 12,
  height = 4,
  units = "in"
)


# Summarise draws, including diagnostics (can import after saving)
# sum_diags <- summarise_draws(trace_draws)
sum_diags <- readRDS("processed_data/brm_sum.rds")
saveRDS(sum_diags, "processed_data/brm_sum.rds")

# Overall count recovery
plot_rec_overall <- pp_check(b_fit, ndraws = 803)

ggsave(
  "output/obj1/plot_recovery_overall.jpg",
  plot_rec_overall,
  width = 8,
  height = 4,
  units = "in"
)

# Arm-level count recovery
plot_rec_trial <- pp_check(b_fit, type = "intervals", ndraws = 803)

ggsave(
  "output/obj1/plot_recovery_trials.jpg",
  plot_rec_trial,
  width = 8,
  height = 4,
  units = "in"
)

# Treatment effect estimates ---------------------------------------------------

# Extract posterior draws
# Summarise mean and quantiles for credible effect

c_draws_df <- as_draws_df(b_fit) %>% 
  pivot_longer(everything(), names_to = "class", values_to = "draw") %>% 
  group_by(class) %>% 
  reframe(
    mean_log = mean(draw),
    upper_log = quantile(draw, 0.975),
    lower_log = quantile(draw, 0.025),
    mean_or = mean(exp(draw)),
    upper_or = quantile(exp(draw), 0.975),
    lower_or = quantile(exp(draw), 0.025)
  ) %>% 
  filter(grepl("^b_", class)) %>% 
  mutate(class = gsub("b_class_short", "", class))

# Plot log-odds
plot_log <- c_draws_df %>% 
  filter(!class == "b_Intercept") %>% 
  ggplot(aes(x = mean_log, xmin = lower_log, xmax = upper_log, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean log-odds (95% credible intervals)", y = "Treatment class")

ggsave(
  "output/obj1/plot_log_odds.png",
  plot_log,
  width = 8,
  height = 4,
  units = "in"
)

# Plot odds ratio
plot_odds <- c_draws_df %>% 
  filter(!class == "b_Intercept") %>% 
  ggplot(aes(x = mean_or, xmin = lower_or, xmax = upper_or, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

ggsave(
  "output/obj1/plot_odds.png",
  plot_odds,
  width = 8,
  height = 4,
  units = "in"
)

# Plot inverse odds
plot_odds_inverse <- c_draws_df %>% 
  filter(!class == "b_Intercept") %>% 
  ggplot(aes(x = 1/mean_or, xmin = 1/lower_or, xmax = 1/upper_or, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

ggsave(
  "output/obj1/plot_odds_inverse.png",
  plot_odds_inverse,
  width = 8,
  height = 4,
  units = "in"
)
