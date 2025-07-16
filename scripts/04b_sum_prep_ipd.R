
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepped summary data
a_imp_sum <- read_csv("output/sum/sum_data.csv")

# Import arm_id lookup
a_imp_lkp <- read_csv("created_metadata/arm_id_trt.csv")

# Isolate ipd trials
b_ipd <- a_imp_sum %>% filter(source == "ipd")
