
source("scripts/00_packages.R")

a_imp_agg <- read_csv("data/trials_n391_dattr.csv")
# a_imp_ipd <- read_csv("data/ipd_n92_dattr.csv")
a_imp_atc <- read_csv("created_metadata/atc.csv")

# Add treatment type
# Set A10Bs in 3 trials as "oad" for now
b_trt <- a_imp_agg %>% 
  select(trial_id, source, atc_new, n:arm_n_attr) %>% 
  mutate(
    atc_short = case_when(
      atc_new == "placebo" ~ "placebo",
      atc_new != "placebo" ~ substr(atc_new, 1, 5)
    )
  ) %>% 
  select(-atc_new) %>% 
  left_join(
    a_imp_atc %>% 
      select(class_short, atc_short = code_short) %>% 
      distinct()
  ) %>% 
  mutate(
    class_short = case_when(
      is.na(class_short) & atc_short == "placebo" ~ "placebo", 
      is.na(class_short) & atc_short == "A10B" ~ "oad",
      TRUE ~ class_short
    ),
    trttype = case_when(
      class_short %in% c("dpp4", "glp1", "sglt2") ~ "exp",
      class_short == "placebo" ~ "plc",
      TRUE ~ "ac"
    )
  )

b_trt %>% count(trial_id, source) %>% count(source)

# Prepare IPD