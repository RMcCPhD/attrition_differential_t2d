
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
a_imp_id <- readRDS("processed_data/res_n92.rds") %>% distinct(nct_id)

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

# Initial plots of interactions ------------------------------------------------

# Commented out after running once

# Model summary
#
# b_sum <- summary(a_imp_mdl)

# Extract interactions summary
#
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

# Plot log-odds
#
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
# Beta: 192,000 rows (168k for age and sex inters, remaining for ref inters)
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
    values_from = c("beta", "delta"),
    names_vary = "slowest"
  )

# Get distribution variables for 92 IPD trials
# Join treatment class names
# Get age10
c_dist_df <- a_imp_df %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  rename(name = trtcls5) %>% 
  left_join(a_imp_atc) %>% 
  inner_join(a_imp_id) %>% 
  select(nct_id:age, class) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Get mean age (no NAs)
c_mean <- c_dist_df %>% 
  group_by(nct_id, sex) %>% 
  reframe(mean_age10_sex = mean(age10))







