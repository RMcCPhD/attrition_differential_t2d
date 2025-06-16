
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Import IPD model exports (coefficients, coefficient correlation)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_ipd_cor <- read_csv("vivli/coef_cor_n92.csv")

# Construct v-cov matrix -------------------------------------------------------

# Trial and model-specific variable sets
b_vars_models <- a_imp_ipd_cor %>% 
  select(nct_id, modeltype, row, col) %>% 
  pivot_longer(row:col, values_to = "variable") %>% 
  distinct(nct_id, modeltype, variable) %>% 
  group_by(nct_id, modeltype) %>% 
  reframe(vars = list(unique(variable)), .groups = "drop")

# Build symmetrical and diagonal entries per model and trial
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

# Convert to v-cov matrices
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
