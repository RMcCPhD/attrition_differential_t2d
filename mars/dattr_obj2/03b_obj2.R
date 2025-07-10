
source("00_config.R")
source("00_packages.R")

# Import aggregate data and prepared ipd outputs
a_imp_ipd_res <- readRDS("data/res_n92.rds")
a_imp_vcov <- readRDS("data/vcov_n92.rds")

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
b_add_covs <- b_prep_res %>% 
  mutate(
    across(sex:age10, ~ if_else(class_short == "placebo", 0, .)),
    sex = case_when(
      int_terms == "sex" | grepl("\\:sex", int_terms) ~ 1,
      TRUE ~ sex
    ),
    age10 = case_when(
      int_terms == "age10" | grepl("\\:age10", int_terms) ~ 1, 
      TRUE ~ age10
    )
  )

# Prepare vcov
b_vcov <- a_imp_vcov %>% filter(modeltype == "glm_int")
b_vcov_ready <- list(b_vcov$cov_matrix)
b_vcov_ready <- b_vcov_ready[[1]]
names(b_vcov_ready) <- unique(b_vcov$nct_id)

# Keep non-NA cells
b_vcov_prep <- map2(b_vcov_ready, names(b_vcov_ready), function(mat, id) {
  
  est_terms <- b_add_covs %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(int_terms)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

c_network <- set_agd_regression(
  data = b_add_covs,
  study = nct_id,
  trt = class_short,
  trt_ref = "placebo",
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt * (age10 + sex),
  trt_class = atc_short
)

plot(c_network)

# With default priors
mdl <- nma(
  c_network,
  trt_effects = "fixed",
  link = "identity",
  regression = ~ .trt * (age10 + sex),
  class_interactions = "common",
  chains = 4, 
  cores = 4,
  control = list(max_treedepth = 15)
)

saveRDS(mdl, "inter_fit.rds")
