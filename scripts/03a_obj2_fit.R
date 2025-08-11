
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

saveRDS(mdl, "output/obj2/inter_fit_n89.rds")
