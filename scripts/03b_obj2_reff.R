
source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model
a_imp_mdl <- readRDS("output/obj2/inter_fit_n90.rds")

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

# Import synthetic IPD (for class-specific age10 means)
a_imp_ipd <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd) %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(readRDS("processed_data/res_n90.rds") %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Check age-sex distributions - normal, similar between sexes
plot_agesex_dist <- a_imp_ipd %>%
  mutate(sex = if_else(sex == 0L, "Female", "Male")) %>% 
  ggplot(aes(x = age10, fill = as.factor(sex))) + 
  geom_density(colour = "black", alpha = 0.3) + 
  facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(n.breaks = 6) +
  theme_classic() +
  labs(x = "Age measured in decades", y = "Density", fill = "Sex")

ggsave(
  "output/plot_agesex_dist_ipd.png",
  plot_agesex_dist,
  width = 8,
  height = 8,
  units = "in"
)

# Get class age10 means
a_ipd_mean <- a_imp_ipd %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  group_by(class) %>% 
  reframe(
    mean_age = mean(age10),
    prop_f = sum(sex == 0L)/sum(sex)
  )

# class mean_age
# <chr>    <dbl>
# 1 dpp4      5.66
# 2 glp1      5.72
# 3 sglt2     5.57

# Examine interaction posterior estimates --------------------------------------

# Get summary data
b_sum_inter <- summary(a_imp_mdl) %>%
  as_tibble() %>%
  filter(grepl("beta", parameter)) %>%
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("beta", "", parameter)
  ) %>%
  separate(
    parameter,
    into = c("term", "class"),
    sep = "\\."
  ) %>%
  mutate(
    class = gsub("trtclass", "", class),
    term = gsub("\\:", "", term),
    term = case_match(
      term,
      "x1TRUE" ~ "Male",
      "x2" ~ "Age",
      .default = term
    ),
    or_mean = exp(mean),
    or_low = exp(`2.5%`),
    or_hi = exp(`97.5%`)
  ) %>% 
  rename(
    logor_mean = mean,
    logor_low = `2.5%`,
    logor_hi = `97.5%`
  ) %>% 
  mutate(across(logor_mean:or_hi, ~ round(., 3))) %>% 
  select(class, term, contains("logor"), starts_with("or"), Bulk_ESS:Rhat)

# Save
write_csv(b_sum_inter, "output/obj2/inter_sum.csv")

# Plot odds ratio
plot_or <- b_sum_inter %>%
  mutate(term = factor(term, levels = c("Age", "Male"))) %>% 
  filter(!is.na(class), class %in% c("dpp4", "glp1", "sglt2")) %>%
  ggplot(aes(x = or_mean, xmin = or_low, xmax = or_hi, y = fct_rev(term))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  facet_wrap(~class, scales = "free", ncol = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0.8, 1.25)) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = NULL)

plot_or

ggsave(
  "output/obj2/plots/plot_inter.png",
  plot_or,
  width = 4,
  height = 4,
  units = "in"
)

# Get posterior data -----------------------------------------------------------

# Get tibble of posterior draws
b_pos <- as.data.frame(a_imp_mdl$stanfit) %>% 
  select(-contains("delta")) %>% 
  as_tibble(rownames = "iter")

summary(b_pos)

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

# Check main effects (treatment)
b_tidy_draws %>% 
  select(class, d_main) %>% 
  group_by(class) %>% 
  reframe(
    mean_main = mean(d_main),
    q2.5_main = quantile(d_main, 0.025),
    q97.5_main = quantile(d_main, 0.975)
  )

# Check age10 and sex interactions
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
  iter = seq(1, 12000, by = 1),
  class = c("dpp4", "glp1", "sglt2"),
  age10 = seq(4, 8, by = 0.1),
  sex   = c(0L, 1L),
  stringsAsFactors = FALSE
)

# Compute relative effects for log-odds and odds ratio
c_pred_reff <- b_tidy_draws %>% 
  inner_join(c_pred_grid %>% mutate(iter = as.character(iter))) %>% 
  mutate(
    reff_log = d_main + b_age * age10 + b_sex * sex,
    reff_or  = exp(reff_log)
  )

# Unadjusted effects
c_unadj <- c("dpp4" = 0.724, "glp1" = 0.941, "sglt2" = 0.699) # All trials
# c_unadj <- c("dpp4" = 0.717, "glp1" = 0.857, "sglt2" = 0.665) # Sens, >=10 attr

# Check difference between average relative effect and unadjusted effects
c_pred_reff %>% 
  group_by(class) %>% 
  reframe(mean_reff = mean(reff_or)) %>% 
  left_join(
    c_unadj %>% 
      enframe() %>% 
      rename(class = name, mean_unadj = value)
  ) %>% 
  mutate(diff = abs(mean_reff - mean_unadj))

# Original version
# class mean_reff mean_unadj   diff
# <chr>     <dbl>      <dbl>  <dbl>
# 1 dpp4      0.713      0.743 0.0294
# 2 glp1      0.848      0.986 0.138 
# 3 sglt2     0.665      0.695 0.0305

# Latest version with fixes to agg data (e.g. A10BX became glp1) and using
# random effects for every model
# class mean_reff mean_unadj   diff
# <chr>     <dbl>      <dbl>  <dbl>
# 1 dpp4      0.709      0.724 0.0147
# 2 glp1      0.853      0.941 0.0883
# 3 sglt2     0.670      0.699 0.0294

# Summarise and plot -----------------------------------------------------------

# Summarised relative effects for odds ratio per class, age and sex
d_sum <- c_pred_reff %>% 
  group_by(class, age10, sex) %>% 
  reframe(
    or_mean  = mean(reff_or),
    or_med   = median(reff_or),
    or_lo    = quantile(reff_or, 0.025),
    or_hi    = quantile(reff_or, 0.975),
  ) %>% 
  mutate(
    across(or_mean:or_hi, ~ round(., 3)),
    class = case_match(
      class,
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2"
    ),
    class = factor(class, levels = c("DPP4", "GLP-1", "SGLT2"))
  )

# Save
write_csv(d_sum, "output/obj2/reff_sum.csv")

# Rug data
ipd_rug <- a_imp_ipd %>%
  transmute(
    class = dplyr::recode(class, dpp4 = "DPP4", glp1 = "GLP-1", sglt2 = "SGLT2"),
    sex,
    age_years = age10 * 10
  )

# Plot females
plot_reff_or_f <- d_sum %>% 
  filter(sex == 0L) %>% 
  mutate(age10 = age10 * 10) %>% 
  ggplot(aes(x = age10, y = or_mean, ymin = or_lo, ymax = or_hi)) +
  geom_line(colour = "steelblue") +
  geom_line(aes(y = or_lo), colour = "steelblue", linetype = "dashed") +
  geom_line(aes(y = or_hi), colour = "steelblue", linetype = "dashed") +
  geom_ribbon(fill = "steelblue", alpha = 0.1) +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(40, 80)) +
  scale_y_continuous(limits = c(0.4, 1.3)) +
  geom_rug(
    data = ipd_rug %>% dplyr::filter(sex == 0L),
    inherit.aes = FALSE,
    aes(x = age_years),
    sides = "b",
    alpha = 0.1,
    colour = "black"
  ) +
  facet_wrap(~ class, scales = "free") +
  theme_classic() +
  labs(
    x = "Mean female age", 
    y = "Relative effect for attrition (odds ratio)"
  )

plot_reff_or_f

ggsave(
  "output/obj2/plots/plot_or_female.png",
  plot_reff_or_f,
  width = 10,
  height = 4,
  units = "in"
)

# Plot males
plot_reff_or_m <- d_sum %>% 
  filter(sex == 1L) %>% 
  mutate(age10 = age10 * 10) %>% 
  ggplot(aes(x = age10, y = or_mean, ymin = or_lo, ymax = or_hi)) +
  geom_line(colour = "steelblue") +
  geom_line(aes(y = or_lo), colour = "steelblue", linetype = "dashed") +
  geom_line(aes(y = or_hi), colour = "steelblue", linetype = "dashed") +
  geom_ribbon(fill = "steelblue", alpha = 0.1) +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(40, 80)) +
  scale_y_continuous(limits = c(0.4, 1.2)) +
  geom_rug(
    data = ipd_rug %>% dplyr::filter(sex == 1L),
    inherit.aes = FALSE,
    aes(x = age_years),
    sides = "b",
    alpha = 0.1,
    colour = "black"
  ) +
  facet_wrap(~ class, scales = "free") +
  theme_classic() +
  labs(
    x = "Mean male age", 
    y = "Relative effect for attrition (odds ratio)"
  )

plot_reff_or_m

ggsave(
  "output/obj2/plots/plot_or_male.png",
  plot_reff_or_m,
  width = 10,
  height = 4,
  units = "in"
)
