
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
b_set_ref <- a_imp_df %>% 
  mutate(
    class = factor(class),
    class = relevel(class, ref = "placebo"),
    arm_compl = arm_n - arm_attr,
    across(arm_n:arm_compl, as.integer)
  )

# Create network
c_agg <- set_agd_arm(
  data = b_set_ref,
  study = trial_id,
  trt = atc,
  trt_ref = "placebo",
  r = arm_attr,
  n = arm_n
)

# Fit model
d_mdl <- nma(
  c_agg,
  trt_effects = "random",
  prior_intercept = normal(scale = 10),
  prior_trt = normal(scale = 10),
  seed = 123,
  chains = 4, 
  cores = 4,
  warmup = 1000,
  iter = 3000
)

# Summarise relative treatment effects
e_sum <- summary(d_mdl$stanfit, pars = "d")$summary %>% 
  as_tibble(rownames = "atc") %>% 
  mutate(atc = gsub("d\\[|\\]", "", atc)) %>% 
  left_join(b_set_ref %>% distinct(class, atc)) %>% 
  mutate(
    mean_or = exp(mean),
    lo_or = exp(`2.5%`),
    hi_or = exp(`97.5%`)
  )

write_csv(e_sum, "output/obj1/sum_res_multinma.csv")

# Plot
plot_odds <- e_sum %>% 
  ggplot(aes(x = mean_or, xmin = lo_or, xmax = hi_or, y = fct_rev(class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

plot_odds

ggsave(
  "output/obj1/plots/plot_odds_multinma.png",
  plot_odds,
  width = 8,
  height = 4,
  units = "in"
)

# Create comparison plot
a_imp_brm <- read_csv("output/obj1/sum_res.csv")

comp_df <- e_sum %>% 
  select(class, mean_or:hi_or) %>% 
  mutate(type = "multinma") %>% 
  full_join(
    a_imp_brm %>% 
      select(class, mean_or, lo_or = lower_or, hi_or = upper_or) %>% 
      mutate(type = "brms")
  ) %>% 
  arrange(class)

comp_plot <- comp_df %>% 
  ggplot(aes(x = mean_or, xmin = lo_or, xmax = hi_or, y = fct_rev(class), colour = type)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

comp_plot

ggsave(
  "output/obj1/plots/plot_odds_comp.png",
  comp_plot,
  width = 8,
  height = 4,
  units = "in"
)
