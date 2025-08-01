
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
# Extract interactions summary
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
    term = gsub("\\:", "", term)
  )
# 
# # Plot log-odds
plot_log <- b_sum_inter %>%
  mutate(
    term = case_match(
      term, 
      "x1TRUE" ~ "sexMale", 
      "x2" ~ "age10",
      .default = term
    )
  ) %>%
  filter(!is.na(class), class %in% c("dpp4", "glp1", "sglt2")) %>%
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = term)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  facet_wrap(~class, scales = "free", ncol = 1) +
  theme_bw() +
  scale_x_continuous(n.breaks = 6) +
  labs(x = "Mean log-odds (95% credible intervals)", y = NULL)

plot_log

ggsave(
  "output/obj2/plots/plot_inter_rejig1.png",
  plot_log,
  width = 4,
  height = 4,
  units = "in"
)

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
    name = gsub("x1TRUE", "sexMale", name),
    name = gsub("x2", "age10", name)
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
# Join reference interactions for age10 and sex
c_join <- c_beta_rm %>% 
  left_join(c_delta_names) %>% 
  left_join(
    c_beta %>% 
      filter(name %in% c("age10", "sexMale")) %>% 
      pivot_wider(names_from = "name", values_from = "value") %>% 
      select(iter, age10, sexMale) %>% 
      rename(age10_ref = age10, sexMale_ref = sexMale)
  )

# NAs in insulin
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
    values_from = "beta",
    names_vary = "slowest"
  ) %>% 
  rename(main = delta)

# Get estimates for newer antidiabetics
# Nest
c_newer <- c_inter %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  rename(age10_est = age10, sexMale_est = sexMale) %>% 
  # group_by(class) %>% 
  # slice_sample(n = 2000) %>% 
  # ungroup() %>% 
  select(-c(iter, contains("_ref")))

# Estimate relative effect -----------------------------------------------------

# Prepare synthetic IPD for newer antidiabetics
# Calculate age10
# Reference is female
d_ipd <- a_imp_df %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(a_imp_id %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Check age-sex distributions - normal, similar between sexes
d_ipd %>%
  ggplot(aes(x = age10, fill = as.factor(sex))) + 
  geom_density(colour = "black", alpha = 0.3) + 
  facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(n.breaks = 6) +
  theme_classic() +
  labs(x = "Age measured in decades", y = "Density", fill = "Sex")

# Create prediction grid
e_pred_grid <- expand_grid(
  trt = c("dpp4", "glp1", "sglt2"),
  age10 = seq(4, 8, by = 0.1),
  sex = c(0L, 1L)
)

# Join prediction grid with posterior draws
e_join_grid <- c_newer %>% 
  crossing(e_pred_grid) %>% 
  mutate(
    reff_log = main + age10 * age10_est + sex * sexMale_est,
    reff_odds = exp(reff_log)
  )

# Check class distribution for relative effects
e_join_grid %>% 
  ggplot(aes(x = reff_log, fill = as.factor(sex))) +
  geom_density(colour = "black", alpha = 0.3) +
  facet_wrap(~age10, scales = "free") +
  theme_classic()

# Summarise relative effect
f_sum_reff <- e_join_grid %>% 
  group_by(class, age10, sex) %>% 
  reframe(
    across(
      contains("reff"),
      list(
        mean = ~ mean(.),
        median = ~ median(.),
        q2.5 = ~ quantile(., 0.025),
        q97.5 = ~ quantile(., 0.975)
      ),
      .names = "{.col}_{.fn}"
    )
  )

# Plot odds
tst_plot <- f_sum_reff %>% 
  mutate(
    sex = case_match(sex, 0L ~ "Female", 1L ~ "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>% 
  rename(Sex = sex) %>% 
  ggplot(aes(x = age10, y = reff_odds_mean, colour = Sex, fill = Sex)) +
  geom_ribbon(
    aes(ymin = reff_odds_q2.5, ymax = reff_odds_q97.5), 
    alpha = 0.05,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_line(aes(y = reff_odds_q2.5), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = reff_odds_q97.5), linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ class, scales = "free") +
  theme_classic() +
  labs(x = "Mean age in decades", y = "Treatment effect on odds of attrition")

ggsave(
  "output/obj2/plots/plot_or_rejig1.png",
  tst_plot,
  width = 10,
  height = 4,
  units = "in"
)


# Check difference between average relative effect and main effect per class
e_join_grid %>% 
  group_by(class) %>% 
  reframe(
    mean_main = mean(main),
    mean_reff = mean(reff_log)
  ) %>% 
  mutate(diff = mean_main - mean_reff)
