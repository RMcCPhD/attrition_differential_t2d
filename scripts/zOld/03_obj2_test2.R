
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n386_with_bl.rds") %>% 
  mutate(
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Import IPD model exports (coefficients, variance-covariance)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Age- and sex-treatment interactions ------------------------------------------

