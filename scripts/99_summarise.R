
# Change file number when additional scripts added for obj2 etc
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )