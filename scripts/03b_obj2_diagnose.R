
summary(a_imp_mdl) %>% 
  as_tibble() %>% 
  filter(grepl("d\\[", parameter)) %>% 
  mutate(parameter = gsub("d\\[|\\]", "", parameter)) %>% 
  left_join(a_imp_atc %>% rename(parameter = name)) %>% 
  select(class, mean, `2.5%`, `97.5%`) %>% 
  ggplot(aes(x = fct_rev(class), y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_point() +
  geom_linerange() +
  geom_hline(yintercept = 0, colour = "red") +
  theme_bw() +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "Treatment class", y = "Mean log-odds (95% credible intervals)") +
  coord_flip()

# Script to explore why estimates for relative effect are strange
# Mean difference between relative effects and main effects are substantial

# Mean difference
# class mean_main mean_reff  diff
# <chr>     <dbl>     <dbl> <dbl>
# 1 dpp4     15.0   0.000243  15.0 
# 2 glp1      0.546 3.53      -2.98
# 3 sglt2     4.35  0.0000426  4.35

# Pooled estimates from brm model of only IPD agg trials
# class mean_log median_log upper_log lower_log mean_or median_or upper_or lower_or
# <chr>    <dbl>      <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>    <dbl>
# 1 dpp4    -0.300     -0.302   -0.167     -0.425   0.742     0.740    0.846    0.653
# 2 glp1    -0.168     -0.167   -0.0183    -0.324   0.848     0.846    0.982    0.723
# 3 sglt2   -0.438     -0.440   -0.300     -0.571   0.647     0.644    0.741    0.565

# First check: NMA of main effects to compare to obj1

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import analysis ready data (res, vcov)
a_imp_res <- readRDS("output/obj2/res_ready_n89.rds")
a_imp_vcov <- readRDS("output/obj2/vcov_ready_n89.rds")

# NMA of main effects from the regression models -------------------------------

# Isolate to main effects
b_main_res <- a_imp_res %>% 
  filter(!grepl(""))

# Construct network
c_network <- set_agd_regression(
  data = a_imp_res,
  study = nct_id,
  trt = trt,
  estimate = estimate,
  cov = a_imp_vcov,
  regression = ~ .trt,
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
  regression = ~ .trt * (x2 + x1),
  prior_intercept = normal(scale = 10),
  prior_trt = normal(scale = 10),
  prior_reg = normal(scale = 10),
  seed = 123,
  chains = 4, 
  cores = 4,
  warmup = 1000,
  iter = 4000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)