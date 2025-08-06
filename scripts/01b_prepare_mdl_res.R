
# This script prepares the model results exported from Vivli
# Model point estimates for interaction models are extracted
# VCOV constructed from correlation data
# Two trials removed because of extreme standard errors

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Import IPD model exports (coefficients, coefficient correlation)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_ipd_cor <- read_csv("vivli/coef_cor_n92.csv")

# Get logistic interaction model sets
b_log_res <- a_imp_ipd_res %>% filter(modeltype == "glm" & spec == "int")
b_log_cor <- a_imp_ipd_cor %>% filter(modeltype == "glm_int")

# Check distributions
hist(b_log_res$estimate) # Some estimates greater than +/- 10
hist(b_log_res$std.error) # Some se very large
hist(b_log_cor$value) # All within +/- 1

# Remove any trials where standard errors are blown up
# 90 trials remaining
b_res_rm <- b_log_res %>% 
  group_by(nct_id) %>% 
  filter(!any(std.error > 11)) %>% 
  ungroup()

# Tidy regression data
# Add class variable
# Join atc lookup to get class codes
c_res_tidy <- b_res_rm %>% 
  select(nct_id, term:std.error) %>% 
  mutate(
    class = term,
    class = case_when(
      class %in% c("age10", "sex", "(Intercept)") ~ NA,
      grepl("\\:sex|\\:age10", class) ~ gsub("\\:sex|\\:age10", "", class),
      TRUE ~ class
    )
  ) %>%
  left_join(a_imp_df %>% select(class, atc) %>% distinct()) %>% 
  select(nct_id:ref, class:atc, everything()) %>% 
  arrange(nct_id)

# Save tidied regression data
saveRDS(c_res_tidy, "processed_data/res_n90.rds")

# Get correlations for 90 trials
c_cor <- b_log_cor %>% 
  inner_join(c_res_tidy %>% select(nct_id) %>% distinct()) %>% 
  select(-modeltype)

# List unique covariate names per trial
c_cor_vars <- c_cor %>% 
  pivot_longer(row:col, values_to = "variable") %>% 
  distinct(nct_id, variable) %>% 
  group_by(nct_id) %>% 
  reframe(vars = list(unique(variable)))

# Bind correlation and rename components for upper and lower triangles
# Create diagonals (i.e. cor = 1)
c_cor_parts <- c_cor %>% 
  bind_rows(c_cor %>% rename(col = row, row = col)) %>% 
  group_by(nct_id) %>% 
  group_split() %>% 
  map_dfr(function(data_part) {
    key <- data_part %>% select(nct_id) %>% slice(1)
    vars <- unique(c(data_part$row, data_part$col))
    
    # Get diagonals
    diags <- expand.grid(row = vars, col = vars) %>% 
      filter(row == col) %>% 
      mutate(value = 1) %>% 
      bind_cols(key[rep(1, nrow(.)), ])
    
    # Bind diagonals and triangles
    bind_rows(data_part, diags)
  })

# Construct correlation matrices
c_cor_construct <- c_cor_parts %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  mutate(
    cor_matrix = map(data, ~ {
      vars <- union(.x$row, .x$col)
      cor_tbl <- xtabs(value ~ row + col, data = .x)[vars, vars]
      matrix(
        as.numeric(cor_tbl),
        nrow = length(vars),
        ncol = length(vars),
        dimnames = list(vars, vars)
      )
    })
  ) %>% 
  select(-data)

# Get VCOV matrices
c_vcov <- c_cor_construct %>% 
  mutate(
    cov_matrix = pmap(list(nct_id, cor_matrix), function(id, cor_mat) {
      se_df <- a_imp_ipd_res %>% filter(nct_id == id)
      se_vec <- setNames(se_df$std.error, se_df$term)
      se_order <- se_vec[rownames(cor_mat)]
      cov_mat <- outer(se_order, se_order) * cor_mat
      cov_mat
    })
  )

# Save
saveRDS(c_vcov, "processed_data/vcov_n90.rds")
