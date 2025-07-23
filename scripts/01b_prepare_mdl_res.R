
# This script prepares the model results exported from Vivli
# Model point estimates for interaction models are extracted
# VCOV constructed from correlation data

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Import IPD model exports (coefficients, coefficient correlation)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_ipd_cor <- read_csv("vivli/coef_cor_n92.csv")

# Prepare ipd aggregate data
b_prep_ipd_res <- a_imp_ipd_res %>% 
  filter(modeltype == "glm" & spec == "int" & std.error < 11) %>% 
  select(nct_id, spec:std.error) %>% 
  rename(class = term) %>% 
  left_join(a_imp_df %>% select(class, atc) %>% distinct())

# Save
saveRDS(b_prep_ipd_res, "processed_data/res_n92.rds")

# Construct VCOV matrices-------------------------------------------------------

# Get logistic interaction model set
# Pivot and nest correlations by trial
# Add variable names as a list per nest
b_vars_models <- a_imp_ipd_cor %>% 
  filter(modeltype == "glm_int") %>% 
  inner_join(b_prep_ipd_res %>% select(nct_id, row = class)) %>% 
  select(nct_id, modeltype, row, col) %>% 
  pivot_longer(row:col, values_to = "variable") %>% 
  distinct(nct_id, modeltype, variable) %>% 
  group_by(nct_id, modeltype) %>% 
  reframe(vars = list(unique(variable)))

# Add triangles and diagonals for each trial to create symmetrical matrices
b_sym <- a_imp_ipd_cor %>% 
  rename(var1 = row, var2 = col) %>% 
  bind_rows(a_imp_ipd_cor %>% rename(var1 = col, var2 = row)) %>% 
  group_by(nct_id, modeltype) %>% 
  group_split() %>% 
  map_dfr(function(data_part) {
    key <- data_part %>% select(nct_id, modeltype) %>% slice(1)
    vars <- unique(c(data_part$var1, data_part$var2))
    
    # Get diagonals
    diags <- expand.grid(var1 = vars, var2 = vars) %>% 
      filter(var1 == var2) %>% 
      mutate(value = 1) %>% 
      bind_cols(key[rep(1, nrow(.)), ])
    
    # Bind diagonals and triangles
    bind_rows(data_part, diags)
  })

# Construct correlation matrices
b_matrices <- b_sym %>% 
  group_by(nct_id, modeltype) %>% 
  nest() %>% 
  mutate(
    cor_matrix = map(data, ~ {
      vars <- union(.x$var1, .x$var2)
      cor_tbl <- xtabs(value ~ var1 + var2, data = .x)[vars, vars]
      matrix(
        as.numeric(cor_tbl),
        nrow = length(vars),
        ncol = length(vars),
        dimnames = list(vars, vars)
      )
    })
  ) %>% 
  select(-data)

# Convert to diagonals and triangles to variance and covariance
b_vcov <- b_matrices %>% 
  mutate(
    cov_matrix = pmap(list(nct_id, modeltype, cor_matrix), function(id, type, cor_mat) {
      se_df <- a_imp_ipd_res %>% filter(nct_id == id, modeltype == modeltype)
      se_vec <- setNames(se_df$std.error, se_df$term)
      se_order <- se_vec[rownames(cor_mat)]
      cov_mat <- outer(se_order, se_order) * cor_mat
      cov_mat
    })
  )

# Save
saveRDS(b_vcov, "processed_data/vcov_n92.rds")
