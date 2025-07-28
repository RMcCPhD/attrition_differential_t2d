
# Note: Reference for sex is female (0/FALSE - female, 1/TRUE - male)
# Objective 2 model checks (convergence, ess, predictive check)
# Initial plotting of log-odds for treatment-covariate interations
# Prepare interactions posterior draws (i.e. join beta and delta per iter)
# Prepare distribution data for IPD trials

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model
a_imp_mdl <- readRDS("output/obj2/inter_fit.rds")

# Import age distributions for IPD trials
a_imp_df <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd)

# Import interaction results to get IPD trial IDs
a_imp_id <- readRDS("processed_data/res_n92.rds")

# Import atc metadata
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

summary(a_imp_mdl)

# Initial plots of interactions ------------------------------------------------

# Commented out after running once
# Model summary
# b_sum <- summary(a_imp_mdl)
# 
# # Extract interactions summary
# b_sum_inter <- b_sum %>%
#   as_tibble() %>%
#   filter(grepl("beta", parameter)) %>%
#   mutate(
#     parameter = gsub("\\[|\\]", "", parameter),
#     parameter = gsub("beta", "", parameter)
#   ) %>%
#   separate(
#     parameter,
#     into = c("term", "class"),
#     sep = "\\."
#   ) %>%
#   mutate(
#     class = gsub("trtclass", "", class),
#     term = gsub("\\:", "", term)
#   )
# 
# # Plot log-odds
# plot_log <- b_sum_inter %>%
#   mutate(term = case_match(term, "sexTRUE" ~ "sexMale", .default = term)) %>%
#   filter(!is.na(class)) %>%
#   ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = term)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   geom_linerange(position = position_dodge(width = 0.5)) +
#   geom_vline(xintercept = 0, colour = "red") +
#   geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
#   geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
#   facet_wrap(~class, scales = "free") +
#   theme_bw() +
#   scale_x_continuous(n.breaks = 6) +
#   labs(x = "Mean log-odds (95% credible intervals)", y = NULL)
# 
# plot_log
# 
# ggsave(
#   "output/obj2/plots/plot_log_odds.png",
#   plot_log,
#   width = 8,
#   height = 4,
#   units = "in"
# )

# Prepare interaction comparisons data -----------------------------------------

# Get posterior draws as tibble
c_pos <- as.data.frame(a_imp_mdl$stanfit) %>% as_tibble(rownames = "iter")

# Identify parameters
c_params <- c_pos %>% 
  pivot_longer(-iter) %>% 
  mutate(
    param = str_extract(name, "^[^\\[]+"),
    name = case_when(
      grepl("NCT", name) ~ gsub("mu\\[|\\]", "", name), 
      grepl("beta|^d", name) ~ gsub("(beta|d)\\[|\\]", "", name),
      TRUE ~ name
    ),
    name = gsub("\\.trtclass", "", name),
    name = gsub("sexTRUE", "sexMale", name)
  )

# Check contents (12k draws max - 3k across 4 chains after warmup)
# Beta: 192,000 rows (168k for age and sex inters, remaining ref inter and main)
# Delta: 84,000 rows (12k for 7 treatment classes)
# Mu: 1,104,000 rows (12k for each trial)
# Log-likelihood: 12,000 (one per iteration)
c_params %>% count(param)

# Nest posterior data
c_pivot <- c_params %>% group_by(param) %>% nest(.key = "iter")

# Pull out beta and delta
c_beta <- c_pivot %>% filter(param == "beta") %>% unnest(iter) %>% ungroup()
c_delta <- c_pivot %>% filter(param == "d") %>% unnest(iter) %>% ungroup()

# Remove reference interactions
# Extract treatment class name to join delta
c_beta_rm <- c_beta %>% 
  filter(name != "age10" & name != "sexMale") %>% 
  mutate(class = str_extract(name, "(?<=:)\\w+")) %>% 
  left_join(a_imp_atc) %>% 
  select(-param) %>% 
  rename(beta = value)

# Replace atc codes with treatment class in delta nest
c_delta_names <- c_delta %>% 
  left_join(a_imp_atc) %>% 
  select(iter, class, delta = value)

# Join beta and delta (no NAs)
c_join <- c_beta_rm %>% left_join(c_delta_names)
checkNA(c_join)

# Extract interactions and rename terms
c_inter <- c_join %>% 
  mutate(
    inter = str_extract(name, "^[^:]+"),
    iter = as.integer(iter)
  ) %>% 
  select(-name) %>% 
  pivot_wider(
    names_from = "inter", 
    values_from = c("beta"),
    names_vary = "slowest"
  ) %>% 
  rename(main = delta)

# Get summarised parameter values
c_inter_sum <- c_inter %>% 
  group_by(class) %>% 
  reframe(
    across(
      main:sexMale,
      list(
        mean = ~ mean(.),
        median = ~ median(.),
        q2.5 = ~ quantile(., 0.025),
        q97.5 = ~ quantile(., 0.975)
      ),
      .names = "{.col}_{.fn}"
    )
  )

c_inter_sum %>% mutate(across(contains("main"), ~ exp(.)))

# Get distribution variables for 92 IPD trials
# Join treatment class names
# Get age10
c_dist_df <- a_imp_df %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(a_imp_id %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  select(nct_id:age, class) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

a_imp_id %>% anti_join(c_dist_df %>% select(nct_id))

c_dist_df %>% 
  group_by(class, sex) %>% 
  reframe(n = n(), mean_age = mean(age10), sd_age = sd(age10))

# Get mean age per sex
# Centre on mean
c_mean_age10_dpp4 <- mean(c_dist_df$age10[c_dist_df$class == "dpp4"])
c_mean_age10_glp1 <- mean(c_dist_df$age10[c_dist_df$class == "glp1"])
c_mean_age10_sglt2 <- mean(c_dist_df$age10[c_dist_df$class == "sglt2"])

c_dist_sum <- c_dist_df %>% 
  mutate(sex = if_else(sex == 0L, 1L, 0L)) %>% 
  group_by(nct_id, sex, class) %>% 
  reframe(
    age10 = mean(age10),
    age10s = case_when(
      class == "dpp4" ~ age10 - c_mean_age10_dpp4,
      class == "glp1" ~ age10 - c_mean_age10_glp1,
      TRUE ~ age10 - c_mean_age10_sglt2
    )
  ) %>% 
  distinct()

# Check age-sex distribution
c_dist_sum %>%
  ggplot(aes(x = age10, fill = as.factor(sex))) + 
  geom_density(colour = "black", alpha = 0.3) + 
  facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(n.breaks = 6) +
  theme_classic() +
  labs(x = "Age measured in decades", y = "Density", fill = "Sex")

# Join summarised parameters
# Calculate treatment effect
c_join <- c_dist_sum %>% 
  left_join(c_inter_sum) %>% 
  mutate(
    effect = age10s * age10_mean + main_mean + sex * sexMale_mean,
    effect_lcri = age10s * age10_q2.5 + main_q2.5 + sex * sexMale_q2.5,
    effect_ucri = age10s * age10_q97.5 + main_q97.5 + sex * sexMale_q97.5,
    across(effect:effect_ucri, ~ exp(.), .names = "{.col}_or")
  ) %>% 
  select(nct_id:age10s, effect:effect_ucri_or)

c_join %>% group_by(class) %>% reframe(mean_eff = mean(exp(effect)))

# Test plots
tst_plot <- c_join %>% 
  mutate(
    sex = case_match(sex, 0L ~ "Female", 1L ~ "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>% 
  rename(Sex = sex) %>% 
  ggplot(aes(x = age10, y = effect, colour = Sex, fill = Sex)) +
  geom_ribbon(
    aes(ymin = effect_lcri, ymax = effect_ucri), 
    alpha = 0.05,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_line(aes(y = effect_lcri), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = effect_ucri), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ class, scales = "free") +
  scale_x_continuous(limits = c(4, 8)) +
  theme_classic() +
  labs(x = "Mean age in decades", y = "Treatment effect on log-odds of attrition")

tst_plot

tst_plot2 <- c_join %>% 
  mutate(
    sex = case_match(sex, 0L ~ "Female", 1L ~ "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>% 
  rename(Sex = sex) %>% 
  ggplot(aes(x = age10, y = effect_or, colour = Sex, fill = Sex)) +
  geom_ribbon(
    aes(ymin = effect_lcri_or, ymax = effect_ucri_or), 
    alpha = 0.05,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_line(aes(y = effect_lcri_or), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = effect_ucri_or), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ class, scales = "free") +
  scale_x_continuous(limits = c(4, 8)) +
  theme_classic() +
  labs(x = "Mean age in decades", y = "Treatment effect on odds of attrition")

tst_plot2

ggsave(
  "output/obj2/plots/plot_inter_effect_OR2.png",
  tst_plot,
  width = 10,
  height = 4,
  units = "in"
)




