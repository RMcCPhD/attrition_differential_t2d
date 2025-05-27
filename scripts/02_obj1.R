
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

summary(b_fit)

# Posterior checks -------------------------------------------------------------

pp_check(b_fit, ndraws = 803)
pp_check(b_fit, type = "intervals", ndraws = 803)

c_draws_df <- as_draws_df(b_fit) %>% 
  pivot_longer(everything(), names_to = "class", values_to = "draw") %>% 
  group_by(class) %>% 
  reframe(
    mean = mean(draw),
    upper = quantile(draw, 0.975),
    lower = quantile(draw, 0.025)
  ) %>% 
  filter(grepl("^b_", class)) %>% 
  mutate(class = gsub("b_class_short", "", class))

c_draws_df %>% 
  filter(!class == "b_Intercept") %>% 
  ggplot(
    aes(
      x = mean, xmin = lower, xmax = upper, 
      y = fct_rev(class)
    )
  ) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(n.breaks = 20) +
  labs(x = "Mean log-odds (95% credible intervals)", y = "Treatment class")
