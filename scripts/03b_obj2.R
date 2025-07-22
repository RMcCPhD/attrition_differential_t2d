
source("scripts/00_packages.R")
source("scripts/00_config.R")

# Import aggregate data and prepared ipd outputs
a_imp_ipd_res <- readRDS("processed_data/res_n92.rds")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Prepare interaction results
b_prep_res <- a_imp_ipd_res %>% 
  select(-spec, -ref, -std.error) %>% 
  bind_rows(
    a_imp_ipd_res %>% 
      group_by(nct_id) %>% 
      transmute(
        nct_id,
        class_short = "placebo",
        atc_short = "placebo",
        estimate = NA
      )
  ) %>% 
  mutate(
    int_terms = class_short,
    class_short = case_when(
      grepl("age10|sex|\\:", class_short) ~ NA,
      class_short == "(Intercept)" ~ NA,
      TRUE ~ class_short
    ),
    age10 = NA,
    sex = NA
  ) %>% 
  arrange(nct_id) %>% 
  select(estimate, class_short, sex, age10, everything()) %>% 
  distinct()

# Construct nma dataset - identify interpretations
# Reference treatment and reference covariate levels
# Intercept, treatment effect and main covariate effects
# Interactions of age and sex on treatment effect 
b_add_constr <- b_prep_res %>% 
  mutate(
    sex = if_else(class_short == "placebo", FALSE, sex),
    age10 = if_else(class_short == "placebo", 0, age10),
    sex = case_when(
      int_terms == "sex" | grepl("\\:sex", int_terms) ~ TRUE,
      TRUE ~ sex
    ),
    age10 = case_when(
      int_terms == "age10" | grepl("\\:age10", int_terms) ~ 1, 
      TRUE ~ age10
    )
  ) %>% 
  rename(trt = atc_short) %>% 
  select(estimate, trt, sex, age10, nct_id, class_short, int_terms)

# Add treatment and class_short to interaction terms
b_add_trt <- b_add_constr %>% 
  group_by(nct_id) %>% 
  mutate(
    trt = case_when(
      is.na(trt) & grepl("biguanide\\:", int_terms) ~ "A10BA",
      is.na(trt) & grepl("sulf\\:", int_terms) ~ "A10BB",
      is.na(trt) & grepl("agluc\\:", int_terms) ~ "A10BF",
      is.na(trt) & grepl("thia\\:", int_terms) ~ "A10BG",
      is.na(trt) & grepl("dpp4\\:", int_terms) ~ "A10BH",
      is.na(trt) & grepl("glp1\\:", int_terms) ~ "A10BJ",
      is.na(trt) & grepl("sglt2\\:", int_terms) ~ "A10BK",
      TRUE ~ trt
    ),
    class_short = case_when(
      is.na(class_short) & grepl("biguanide\\:", int_terms) ~ "biguanide",
      is.na(class_short) & grepl("sulf\\:", int_terms) ~ "sulf",
      is.na(class_short) & grepl("agluc\\:", int_terms) ~ "agluc",
      is.na(class_short) & grepl("thia\\:", int_terms) ~ "thia",
      is.na(class_short) & grepl("dpp4\\:", int_terms) ~ "dpp4",
      is.na(class_short) & grepl("glp1\\:", int_terms) ~ "glp1",
      is.na(class_short) & grepl("sglt2\\:", int_terms) ~ "sglt2",
      TRUE ~ class_short
    )
  )

# Prepare vcov
b_vcov <- a_imp_vcov %>% filter(modeltype == "glm_int")
b_vcov_ready <- list(b_vcov$cov_matrix)
b_vcov_ready <- b_vcov_ready[[1]]
names(b_vcov_ready) <- unique(b_vcov$nct_id)

# Keep non-NA cells
b_vcov_prep <- map2(b_vcov_ready, names(b_vcov_ready), function(mat, id) {
  
  est_terms <- b_add_constr %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(int_terms)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

c_network <- set_agd_regression(
  data = b_add_constr,
  study = nct_id,
  trt = trt,
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt * (age10 + sex),
  trt_ref = "placebo",
  trt_class = class_short
)

plot(c_network)

# With default priors
mdl <- nma(
  c_network,
  trt_effects = "fixed",
  link = "identity",
  class_interactions = "common",
  likelihood = "normal",
  regression = ~ .trt * (age10 + sex),
  prior_intercept = normal(scale = 1),
  prior_trt = normal(scale = 1),
  prior_reg = normal(scale = 1),
  seed = 123,
  chains = 4, 
  cores = 4,
  warmup = 1000,
  iter = 4000,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

saveRDS(mdl, "inter_fit.rds")