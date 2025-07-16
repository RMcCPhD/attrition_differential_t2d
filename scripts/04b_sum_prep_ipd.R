
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepped summary data
a_imp_sum <- read_csv("output/sum/sum_data.csv")

# Import arm_id lookup
a_imp_lkp <- read_csv("created_metadata/arm_id_trt.csv")

# Isolate ipd trials
b_ipd <- a_imp_sum %>% filter(source == "ipd")

# Import baseline characteristics from nma ipd trials
b_files_bl <- list.files("data/nma_agesex", pattern = "age_dist", recursive = TRUE)
b_files_path <- paste0("data/nma_agesex/", b_files_bl)
b_files_lst <- lapply(b_files_path, read.csv)
b_bl <- bind_rows(b_files_lst) %>% as_tibble() %>% select(nct_id:sd)

# Get dattr trials
b_bl_dattr <- b_ipd %>%
  left_join(a_imp_lkp) %>% 
  select(trial_id:class, arm_id, n, n_male) %>% 
  left_join(b_bl %>% rename(trial_id = nct_id, arm_id = arm_id_unq))

# Isolate dattr trials with missing baseline
# Mismatch in arm_id_unq and nma ipd trials
c_na <- b_bl_dattr %>% filter(is.na(sex)) %>% select(trial_id) %>% distinct()
c_bl <- b_bl %>% filter(nct_id %in% c_na$trial_id)
