
# This script is run within a high performance computing environment
# The slurm file used to run the script has the same file name
# Sets up the interactions networks and fits the network meta-analysis
# Common interactions and fixed treatment effects assumed across trials

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import aggregate data and prepared ipd outputs
a_imp_ipd_res <- readRDS("processed_data/res_n92.rds")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Import atc lookup
a_imp_atc <- read_csv("created_metadata/atc.csv")
a_imp_atc %>% count(code_short, class_short)

# Prepare regression dataset ---------------------------------------------------

# Add reference rows to regression dataset
b_add_plc <- a_imp_ipd_res %>% 
  bind_rows(
    a_imp_ipd_res %>% 
      group_by(nct_id) %>% 
      transmute(nct_id, class = ref, atc = ref, estimate = NA, ref) %>% 
      ungroup()
  ) %>% 
  distinct() %>% 
  arrange(nct_id) %>% 
  select(-c(spec, std.error))

# Create a terms variable
# Remove interaction terms from class for covariate effects
# Trim class to name for interaction terms
# Add empty age10 and sex
b_sep_terms <- b_add_plc %>% 
  mutate(
    terms = class,
    class = case_when(
      grepl("^sex|^age10", class) ~ NA,
      grepl("\\:", class) ~ gsub("\\:age10|\\:sex", "", class),
      class == "(Intercept)" ~ NA,
      TRUE ~ class
    )
  )

# Add covariate levels
# Complete treatment codes, corrects reference where not placebo
b_add_levels <- b_sep_terms %>% 
  mutate(
    atc = case_when(
      class == "insulin" ~ "A10A",
      class == "biguanide" ~ "A10BA",
      class == "sulf" ~ "A10BB",
      class == "agluc" ~ "A10BF",
      class == "thia" ~ "A10BG",
      class == "dpp4" ~ "A10BH",
      class == "glp1" ~ "A10BJ",
      class == "sglt2" ~ "A10BK",
      TRUE ~ atc
    ),
    sex = case_when(
      grepl("sex", terms) ~ TRUE, 
      class == ref ~ FALSE,
      TRUE ~ NA
    ),
    age10 = case_when(
      grepl("age10", terms) ~ 1, 
      class == ref ~ 0,
      TRUE ~ NA
    )
  )

# Keep core variables and rename to mimic vignette example
b_ready <- b_add_levels %>% 
  rename(trt = atc, x1 = sex, x2 = age10) %>% 
  select(estimate, trt, x1, x2, nct_id, class, terms, ref)

dpp4 <- b_ready %>% group_by(nct_id) %>% filter(any(class == "dpp4")) %>% ungroup()
glp1 <- b_ready %>% group_by(nct_id) %>% filter(any(class == "glp1")) %>% ungroup()
sglt2 <- b_ready %>% group_by(nct_id) %>% filter(any(class == "sglt2")) %>% ungroup()

# Test: Remove placebo rows where it was not the reference
b_test <- b_ready %>% 
  group_by(nct_id) %>% 
  mutate(terms = if_else(ref != "placebo" & terms == "placebo", NA, terms))

# Prepare variance-covariance matrices -----------------------------------------

b_vcov <- a_imp_vcov %>% filter(modeltype == "glm_int")
b_vcov_ready <- list(b_vcov$cov_matrix)
b_vcov_ready <- b_vcov_ready[[1]]
names(b_vcov_ready) <- unique(b_vcov$nct_id)

# Keep non-NA cells
b_vcov_prep <- map2(b_vcov_ready, names(b_vcov_ready), function(mat, id) {
  
  est_terms <- b_ready %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(terms)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

c_network <- set_agd_regression(
  data = b_ready,
  study = nct_id,
  trt = trt,
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt * (x2 + x1),
  trt_ref = "placebo",
  trt_class = class
)

plot(c_network)

# With default priors
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

saveRDS(mdl, "inter_fit.rds")