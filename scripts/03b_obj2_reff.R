
source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model
a_imp_mdl <- readRDS("output/obj2/inter_fit_n89.rds")

# Import ATC lookup
a_imp_atc <- read_csv("created_metadata/updated_class_names_codes.csv") %>% 
  mutate(
    new_class = case_match(
      new_class,
      "a_gluc" ~ "agluc",
      "biguanides" ~ "biguanide",
      .default = new_class
    )
  ) %>% 
  select(class = new_class, name = atc_short)

# Import IPD (for class-specific age10 means)
a_imp_ipd <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd) %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(readRDS("processed_data/res_n90.rds") %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Get class age10 means
a_ipd_mean <- a_imp_ipd %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  group_by(class) %>% 
  reframe(mean_age = mean(age10))

# class mean_age
# <chr>    <dbl>
# 1 dpp4      5.66
# 2 glp1      5.72
# 3 sglt2     5.57

# Get posterior data -----------------------------------------------------------

# Get tibble of posterior draws
b_pos <- as.data.frame(a_imp_mdl$stanfit) %>% as_tibble(rownames = "iter")

# Build a tidy table of draws for each newer antidiabetic
get_names <- c(dpp4 = "A10BH", glp1 = "A10BJ", sglt2 = "A10BK")

get_draws <- function(class_name, atc_name) {
  tibble(
    iter = b_pos$iter,
    class = class_name,
    d_main = b_pos[[sprintf("d[%s]", atc_name)]],
    b_age = b_pos[[sprintf("beta[x2:.trtclass%s]", class_name)]],
    b_sex = b_pos[[sprintf("beta[x1TRUE:.trtclass%s]", class_name)]]
  )
}

b_tidy_draws <- map2_dfr(names(get_names), unname(get_names), get_draws)

# Check slopes for age10 and sex
b_tidy_draws %>% 
  group_by(class) %>% 
  reframe(
    mean_age10 = mean(b_age),
    q2.5_age10 = quantile(b_age, 0.025),
    q97.5_age10 = quantile(b_age, 0.975),
    mean_sex = mean(b_sex),
    q2.5_sex = quantile(b_sex, 0.025),
    q97.5_sex = quantile(b_sex, 0.975)
  )

# Get relative effects ---------------------------------------------------------

# Prediction grid - age10 4-8, sex 0 and 1
c_pred_grid <- expand.grid(
  class = c("dpp4", "glp1", "sglt2"),
  age10 = seq(4, 8, by = 0.1),
  sex   = c(0L, 1L)
)

# Compute relative effects for log-odds and odds ratio
c_pred_reff <- b_tidy_draws %>% 
  inner_join(c_pred_grid) %>% 
  mutate(
    reff_log = d_main + b_age * age10 + b_sex * sex,
    reff_or  = exp(reff_log)
  )

# Check difference between average relative effect and main effect per class
c_pred_reff %>% 
  group_by(class) %>% 
  reframe(
    mean_main = mean(exp(d_main)),
    mean_reff = mean(reff_or)
  ) %>% 
  mutate(diff = mean_main - mean_reff)

# Summarise and plot -----------------------------------------------------------

# Summarised relative effects for odds ratio per class, age and sex
d_sum <- c_pred_reff %>% 
  group_by(class, age10, sex) %>% 
  reframe(
    or_mean  = mean(reff_or),
    or_med   = median(reff_or),
    or_lo    = quantile(reff_or, 0.025),
    or_hi    = quantile(reff_or, 0.975),
  )

# Plot
plot_reff_or <- d_sum %>% 
  mutate(
    sex = case_match(sex, 0L ~ "Female", 1L ~ "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>% 
  rename(Sex = sex) %>% 
  ggplot(aes(x = age10, y = or_mean, colour = Sex, fill = Sex)) +
  geom_ribbon(
    aes(ymin = or_lo, ymax = or_hi), 
    alpha = 0.05,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_line(aes(y = or_lo), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = or_hi), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ class, scales = "free") +
  theme_classic() +
  labs(x = "Mean age in decades", y = "Treatment effect on odds of attrition")

plot_reff_or

ggsave(
  "output/obj2/plots/plot_or_new.png",
  plot_reff_or,
  width = 10,
  height = 4,
  units = "in"
)
