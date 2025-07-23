
# Objective 1 model checks (convergence, ess, predictive check)
# Plotting pooled odds of attrition versus placebo

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model
a_fit <- readRDS("output/obj1/brm_fit.rds")

# Model summary
summary(a_fit)

# Extract posterior draws and summarise
a_draws <- as_draws_df(a_fit)
a_draws_sum <- summarise_draws(a_draws)

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
    mean_log = mean(draw),
    median_log = median(draw),
    upper_log = quantile(draw, 0.975),
    lower_log = quantile(draw, 0.025),
    mean_or = mean(exp(draw)),
    median_or = median(exp(draw)),
    upper_or = quantile(exp(draw), 0.975),
    lower_or = quantile(exp(draw), 0.025)
  ) %>% 
  filter(grepl("^b_", class)) %>% 
  mutate(class = gsub("b_class", "", class))
  
# Save summaries
write_csv(a_draws_cnvg, "output/obj1/sum_cnvg.csv")
write_csv(a_draws_split, "output/obj1/sum_res.csv")

# Predictive check on arm-level attrition counts
plot_rec_trial <- pp_check(a_fit, type = "intervals", ndraws = 803)

ggsave(
  "output/obj1/ppc_arm_level.jpg",
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
plot_odds <- a_draws_split %>% 
  filter(!class == "b_Intercept") %>% 
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
