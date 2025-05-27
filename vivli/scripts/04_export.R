
source("Scripts/00_config.R")
library(tidyverse)

# Import model outputs and aggregate data
a_imp_res <- readRDS("Outputs/tidy_mdl_res.rds")
a_imp_diag <- readRDS("Outputs/tidy_mdl_diag.rds")
a_imp_vcov <- readRDS("Outputs/tidy_mdl_vcov.rds")
a_imp_agg <- read_csv("Processed_data/agg_n92.csv")

# Check sponsor frequency
a_imp_agg %>% count(sponsor, nct_id) %>% count(sponsor) %>% mutate(pcnt = n/92 * 100)

# Tidy results
b_res <- a_imp_res %>% mutate(across(estimate:uci, ~ round(., 3)))

# Tidy diagnostics
b_diag <- a_imp_diag %>% 
  ungroup() %>% 
  mutate(across(statistic.log:df.residual, ~ round(., 3)))

# Tidy long-format correlation
b_cor <- a_imp_vcov %>% 
  mutate(
    value = round(value, 3),
    across(row:col, ~ gsub("ref_class", "", .)),
    across(row:col, ~ gsub("_[0-9]", "", .))
  )

# Save
write_csv(a_imp_agg, "Export/agg_n92.csv")
write_csv(b_res, "Export/res_n92.csv")
write_csv(b_diag, "Export/diag_n92.csv")
write_csv(b_cor, "Export/coef_cor_n92.csv")
