
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import fitted model
a_fit <- readRDS("processed_data/brm_fit.rds")
a_draws <- as_draws_df(a_fit)
a_draws_sum <- summarise_draws(a_draws)

# Model summary
summary(a_fit)

# Convergence summary
a_draws_cnvg <- a_draws_sum %>% 
  filter(grepl("^b_", variable)) %>% 
  mutate(variable = gsub("^b_|^b_class_short", "", variable)) %>% 
  select(variable, rhat:ess_tail)

# Tidy draws summary
# Summarise mean and quantiles for credible effect
a_draws_split <- a_draws %>% 
  pivot_longer(everything(), names_to = "class", values_to = "draw") %>% 
  group_by(class) %>% 
  reframe(
    median_log = median(draw),
    upper_log = quantile(draw, 0.975),
    lower_log = quantile(draw, 0.025),
    mean_or = median(exp(draw)),
    upper_or = quantile(exp(draw), 0.975),
    lower_or = quantile(exp(draw), 0.025)
  ) %>% 
  filter(grepl("^b_", class)) %>% 
  mutate(class = gsub("b_class_short", "", class))
  
# Save summaries
write_csv(a_draws_cnvg, "output/obj1/sum_cnvg.csv")
write_csv(a_draws_split, "output/obj1/sum_res.csv")

# Trace plots
plot_trace <- mcmc_trace(a_draws, pars = vars(starts_with("b_")))

ggsave(
  "output/obj1/trace_treatment_effects.jpg",
  plot_trace,
  width = 12,
  height = 8,
  units = "in"
)


# Overall count recovery
plot_rec_overall <- pp_check(a_fit, ndraws = 803)

ggsave(
  "output/obj1/plot_recovery_overall.jpg",
  plot_rec_overall,
  width = 6,
  height = 4,
  units = "in"
)

# Arm-level count recovery
plot_rec_trial <- pp_check(a_fit, type = "intervals", ndraws = 803)

ggsave(
  "output/obj1/plot_recovery_trials.jpg",
  plot_rec_trial,
  width = 6,
  height = 4,
  units = "in"
)

# Treatment effect estimates ---------------------------------------------------

# Plot log-odds
plot_log <- a_draws_split %>% 
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
plot_odds <- a_draws_split %>% 
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
plot_odds_inverse <- a_draws_split %>% 
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
