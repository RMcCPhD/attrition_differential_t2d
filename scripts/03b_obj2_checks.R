

# This script sets up the interactions network and fits the NMA
# Common interactions and fixed treatment effects assumed across trials

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import aggregate data and prepared ipd outputs
a_imp_res <- readRDS("processed_data/res_n90.rds")
a_imp_vcov <- readRDS("processed_data/vcov_n90.rds")

# Import atc lookup
a_imp_atc <- read_csv("created_metadata/atc.csv")
a_imp_atc %>% count(code_short, class_short)

# Add reference rows to regression dataset
b_add_ref <- a_imp_res %>% 
  bind_rows(
    a_imp_res %>% 
      group_by(nct_id) %>% 
      transmute(nct_id, class = ref, atc = ref, estimate = NA, ref) %>% 
      ungroup()
  ) %>% 
  distinct() %>% 
  arrange(nct_id)

# Fix atc where reference is class name
# Add covariate levels
b_add_covs <- b_add_ref %>% 
  mutate(
    atc = case_when(
      class == "insulin"   | grepl("insulin\\:", term) ~ "A10A",
      class == "biguanide" | grepl("biguanide\\:", term)  ~ "A10BA",
      class == "sulf"      | grepl("sulf\\:", term)  ~ "A10BB",
      class == "agluc"     | grepl("agluc\\:", term)  ~ "A10BF",
      class == "thia"      | grepl("thia\\:", term)  ~ "A10BG",
      class == "dpp4"      | grepl("dpp4\\:", term)  ~ "A10BH",
      class == "glp1"      | grepl("glp1\\:", term)  ~ "A10BJ",
      class == "sglt2"     | grepl("sglt2\\:", term)  ~ "A10BK",
      TRUE ~ atc
    ),
    sex = case_when(
      grepl("sex", term) ~ TRUE, 
      class == ref ~ FALSE,
      TRUE ~ NA
    ),
    age10 = case_when(
      grepl("age10", term) ~ 1, 
      class == ref ~ 0,
      TRUE ~ NA
    )
  )

# Inspect interaction estimates
# Most estimates cross zero
plot_inspect_res <- b_add_covs %>%
  mutate(
    term = if_else(term == "(Intercept)", paste0("aaIntercept_", ref), term),
    term = as.factor(term),
    term = fct_relevel(term, grep("aaIntercept_", levels(term), value = TRUE)[1])
  ) %>%
  arrange(nct_id, term) %>%
  select(nct_id, term, estimate, std.error) %>%
  na.omit() %>%
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error
  ) %>%
  ggplot(aes(x = fct_rev(term), y = estimate, ymin = lci, ymax = uci)) +
  geom_point() +
  geom_linerange() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~nct_id, scales = "free") +
  coord_flip()

# ggsave(
#   "output/obj2/plots/inspect_res.png",
#   plot_inspect_res,
#   width = 30,
#   height = 20,
#   units = "in"
# )

# Keep core variables and rename to mimic vignette example
b_ready <- b_add_covs %>% 
  rename(trt = atc, x1 = sex, x2 = age10) %>% 
  select(estimate, std.error, trt, x1, x2, nct_id, class, term, ref)

# Prepare vcov matrices
c_vcov_lst <- a_imp_vcov$cov_matrix
names(c_vcov_lst) <- unique(b_ready$nct_id)

# Remove a trial with very large variance for a treatment:sex interaction
# 89 trials remaining
b_ready_rm <- b_ready %>% filter(nct_id != "NCT02597049")
c_vcov_rm <- c_vcov_lst[names(c_vcov_lst) != "NCT02597049"]

# Get spread of trial-specific slopes for age10 and sex
# This is the sticking point - pooled estimates for dpp4 and glp1 are ~ -2.5
# Estimates nowhere near that size among individual regressions
b_ready_rm %>%
  filter(str_detect(term, ":age")) %>%
  group_by(class) %>%
  summarise(
    mean = mean(estimate, na.rm=TRUE),
    sd   = sd(estimate, na.rm=TRUE),
    p10  = quantile(estimate, .10, na.rm=TRUE),
    p50  = quantile(estimate, .50, na.rm=TRUE),
    p90  = quantile(estimate, .90, na.rm=TRUE)
  )

b_ready_rm %>%
  filter(str_detect(term, ":sex")) %>%
  group_by(class) %>%
  summarise(
    mean = mean(estimate, na.rm=TRUE),
    sd   = sd(estimate, na.rm=TRUE),
    p10  = quantile(estimate, .10, na.rm=TRUE),
    p50  = quantile(estimate, .50, na.rm=TRUE),
    p90  = quantile(estimate, .90, na.rm=TRUE)
  )

# parameter                           mean      sd   `2.5%`   `25%`   `50%`   `75%`  `97.5%` Bulk_ESS Tail_ESS  Rhat
# <chr>                              <dbl>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>    <dbl> <dbl>
# 1 beta[x2]                         0.0547  0.00754  0.0397   0.0496  0.0547  0.0598  0.0694      964.    2147.  1.00
# 2 beta[x1TRUE]                    -0.161   0.0174  -0.195   -0.173  -0.161  -0.150  -0.127     12045.    9963.  1.00
# 3 beta[x2:.trtclassagluc]         -1.10    0.413   -1.92    -1.37   -1.09   -0.819  -0.289     13047.    9881.  1.00
# 4 beta[x2:.trtclassbiguanide]     -0.887   0.0512  -0.986   -0.922  -0.887  -0.853  -0.786      3308.    5932.  1.00
# 5 beta[x2:.trtclassdpp4]          -2.23    0.0229  -2.28    -2.25   -2.23   -2.22   -2.19       1468.    3695.  1.00
# 6 beta[x2:.trtclassglp1]           0.298   0.0246   0.250    0.281   0.298   0.314   0.346      2649.    6715.  1.00
# 7 beta[x2:.trtclassinsulin]        0.756   0.0490   0.659    0.723   0.756   0.789   0.852      2036.    3759.  1.00
# 8 beta[x2:.trtclasssglt2]         -2.28    0.0279  -2.33    -2.30   -2.28   -2.26   -2.22       2380.    5219.  1.00
# 9 beta[x2:.trtclasssulf]          -3.04    0.0255  -3.09    -3.05   -3.04   -3.02   -2.99       1966.    4235.  1.00
# 10 beta[x2:.trtclassthia]           0.0922  0.0947  -0.0943   0.0279  0.0926  0.157   0.277     10766.    7945.  1.00
# 11 beta[x1TRUE:.trtclassagluc]     -2.29    0.644   -3.56    -2.72   -2.29   -1.86   -1.03      23545.    8546.  1.00
# 12 beta[x1TRUE:.trtclassbiguanide]  0.116   0.0576   0.00243  0.0775  0.116   0.155   0.229     23637.    9039.  1.00
# 13 beta[x1TRUE:.trtclassdpp4]       0.0990  0.0357   0.0285   0.0749  0.0985  0.123   0.168     18770.    9838.  1.00
# 14 beta[x1TRUE:.trtclassglp1]       0.304   0.0567   0.192    0.266   0.304   0.342   0.416     15146.   10253.  1.00
# 15 beta[x1TRUE:.trtclassinsulin]   -0.226   0.0434  -0.313   -0.255  -0.226  -0.197  -0.141     18228.    9917.  1.00
# 16 beta[x1TRUE:.trtclasssglt2]     -0.650   0.0537  -0.757   -0.685  -0.649  -0.614  -0.544     21641.    8775.  1.00
# 17 beta[x1TRUE:.trtclasssulf]      -0.0575  0.0261  -0.109   -0.0752 -0.0572 -0.0399 -0.00607   19518.    9502.  1.00
# 18 beta[x1TRUE:.trtclassthia]       0.00993 0.175   -0.334   -0.107   0.0110  0.127   0.358     23211.    8527.  1.00
# 19 d[A10A]                         -0.506   0.0457  -0.596   -0.537  -0.506  -0.475  -0.417      2043.    3723.  1.00
# 20 d[A10BA]                         0.968   0.0514   0.867    0.933   0.968   1.00    1.07       3748.    6016.  1.00
# 21 d[A10BB]                         3.41    0.0229   3.37     3.40    3.41    3.43    3.46       1742.    4047.  1.00
# 22 d[A10BF]                         0.755   0.332    0.105    0.532   0.752   0.975   1.41      13122.    9673.  1.00
# 23 d[A10BG]                        -0.810   0.0903  -0.986   -0.872  -0.809  -0.748  -0.633     10620.    8071.  1.00
# 24 d[A10BH]                         2.70    0.0212   2.65     2.68    2.70    2.71    2.74       1424.    3615.  1.00
# 25 d[A10BJ]                        -0.649   0.0221  -0.692   -0.664  -0.650  -0.634  -0.606      2428.    5288.  1.00
# 26 d[A10BK]                         1.46    0.0249   1.41     1.44    1.46    1.48    1.51       2243.    4993.  1.00

# Construct network
c_network <- set_agd_regression(
  data = b_ready_rm,
  study = nct_id,
  trt = trt,
  estimate = estimate,
  cov = c_vcov_rm,
  # regression = ~ .trt,
  regression = ~ .trt * (x2 + x1),
  trt_ref = "placebo",
  trt_class = class
)

plot(c_network)

# Fit NMA model
mdl <- nma(
  c_network,
  trt_effects = "fixed",
  link = "identity",
  likelihood = "normal",
  class_interactions = "common",
  # regression = ~ .trt,
  regression = ~ .trt * (x2 + x1),
  prior_intercept = normal(scale = 2),
  prior_trt = normal(scale = 2),
  prior_reg = normal(scale = 1),
  seed = 123,
  chains = 4, 
  cores = 4,
  warmup = 1000,
  iter = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# saveRDS(mdl, "output/obj2/inter_fit_trt_only.rds")
saveRDS(mdl, "output/obj2/inter_fit_n89.rds")


# Note: Reference for sex is female (0/FALSE - female, 1/TRUE - male)
# Objective 2 model checks (convergence, ess, predictive check)
# Initial plotting of log-odds for treatment-covariate interations
# Prepare interactions posterior draws (i.e. join beta and delta per iter)
# Prepare distribution data for IPD trials

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model
a_imp_mdl <- readRDS("output/obj2/inter_fit_n89.rds")
# a_imp_mdl <- readRDS("output/obj2/inter_fit_trt_only.rds")

# Import age distributions for IPD trials
a_imp_df <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd)

# Import interaction results to get IPD trial IDs
a_imp_id <- readRDS("processed_data/res_n90.rds")

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

# Posterior summaries
summary(a_imp_mdl) %>% 
  as_tibble() %>% 
  filter(grepl("beta\\[|d\\[", parameter)) %>% 
  print(n = 26)

# Plot unadjusted treatment effects (model with ~ .trt)
# plot_unadj_trt <- summary(a_imp_mdl) %>%
#   as_tibble() %>% 
#   filter(grepl("d\\[", parameter)) %>% 
#   mutate(
#     parameter = gsub("\\[|\\]", "", parameter),
#     parameter = gsub("d", "", parameter)
#   ) %>%
#   rename(name = parameter) %>% 
#   left_join(a_imp_atc) %>% 
#   filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
#   ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = class)) +
#   geom_point(size = 0.5) +
#   geom_linerange() +
#   geom_vline(xintercept = 0, colour = "red") +
#   theme_bw() +
#   scale_x_continuous(n.breaks = 6) +
#   labs(x = "Mean log-odds (95% credible intervals)", y = NULL)
# 
# ggsave(
#   "output/checks/plot_unadj_nma.png",
#   plot_unadj_trt,
#   width = 4,
#   height = 2,
#   units = "in"
# )

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

# Plot log-odds
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

# ggsave(
#   "output/obj2/plots/plot_inter.png",
#   plot_log,
#   width = 4,
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
    name = gsub("x1TRUE", "sexMale", name),
    name = gsub("x2", "age10", name)
  )

# Check contents (12k draws max - 3k across 4 chains after warmup)
# Beta: 192,000 rows (168k for age and sex inters, remaining ref inter and main)
# Delta: 84,000 rows (12k for 7 treatment classes)
# Mu: 1,104,000 rows (12k for each trial) - 672k in reduced version (56 trials)
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
c_newer <- c_inter %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  rename(age10_est = age10, sexMale_est = sexMale) %>% 
  select(-c(contains("_ref")))

# Estimate relative effect -----------------------------------------------------

# Prepare synthetic IPD for newer antidiabetics
# Get age10
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

# Get mean age per treatment class
d_ipd %>% group_by(class) %>% reframe(mean_age = mean(age10))

# Create prediction grid
e_pred_grid <- expand_grid(
  trt = c("dpp4", "glp1", "sglt2"),
  age10 = seq(4, 8, by = 0.1),
  sex = c(0L, 1L)
  ) %>% 
  mutate(
    age10c = case_when(
      trt == "dpp4" ~ age10 - 5.66,
      trt == "glp1" ~ age10 - 5.72,
      trt == "sglt2" ~ age10 - 5.57
    )
  )

# Join prediction grid with posterior draws
e_join_grid <- c_newer %>% 
  # group_by(class) %>% 
  # slice_sample(n = 2000) %>% 
  # ungroup() %>% 
  crossing(e_pred_grid) %>%
  mutate(
    reff_log = main + age10 * age10_est + sex * sexMale_est,
    reff_odds = exp(reff_log),
  )

# Plot relative effect for estimates calculated per draw
tst_plot_unscaled <- e_join_grid %>% 
  mutate(
    sex = case_match(sex, 0L ~ "Female", 1L ~ "Male"),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>% 
  rename(Sex = sex) %>% 
  ggplot(aes(x = age10, y = reff_odds, colour = Sex, fill = Sex, group = iter)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~class, scales = "free") +
  theme_classic() +
  labs(x = "Mean age in decades", y = "Treatment effect on odds of attrition")

tst_plot_unscaled

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

# Plot odds (not scaled)
tst_plot_unscaled <- f_sum_reff %>% 
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

tst_plot_unscaled

ggsave(
  "output/obj2/plots/plot_or.png",
  tst_plot_unscaled,
  width = 10,
  height = 4,
  units = "in"
)

# Check difference between average relative effect and main effect per class
e_join_grid %>% 
  group_by(class) %>% 
  reframe(
    mean_main = mean(exp(main)),
    mean_reff = mean(reff_odds)
  ) %>% 
  mutate(diff = mean_main - mean_reff)
