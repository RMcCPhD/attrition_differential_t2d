
source("Scripts/00_config.R")
library(tidyverse)

# Import data
a_imp_df <- readRDS("V:/From 8697/work_folder/attrition_attrisk/Processed_data/n95.rds")

# Classify treatment type based on treatment class
# Remove patients missing age and sex (n=83)
# Remove trials comparing one antidiabetic
# 92 trials
b_df <- a_imp_df %>% 
  select(sponsor:event, age10:sex) %>% 
  mutate(
    trttype = case_when(
      trt_class == "placebo" ~ "plc",
      trt_class %in% c("dpp4", "glp1", "sglt2") ~ "exp",
      TRUE ~ "ac"
    )
  ) %>% 
  filter(!is.na(age10), !is.na(sex)) %>% 
  group_by(nct_id) %>% 
  filter(length(unique(trt_class)) != 1) %>% 
  ungroup()

b_df %>% count(nct_id)
b_df %>% count(arm_lvl, trt_class) %>% print(n = 32)
b_df %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save full prepared dataset
saveRDS(b_df, "Processed_data/n92.rds")

# Get aggregated attrition counts per arm in trial
c_agg <- b_df %>% 
  group_by(sponsor, nct_id, trt_class, arm_lvl) %>% 
  reframe(
    n = n(),
    attr = sum(event),
    compl = n - attr,
    min_dur = min(weeks),
    max_dur = max(weeks),
    med_dur = median(weeks)
  ) %>% 
  mutate(
    trttype = case_when(
      trt_class == "placebo" ~ "plc",
      trt_class %in% c("dpp4", "glp1", "sglt2") ~ "exp",
      TRUE ~ "ac"
    ),
    across(contains("dur"), ~ round(., 2))
  ) %>% 
  select(sponsor:trt_class, contains("dur"), arm_lvl, trttype, everything())

c_agg %>% count(nct_id, sponsor) %>% count(sponsor) %>% mutate(pcnt = n/92 * 100)
c_agg %>% count(nct_id)
c_agg %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save aggregated dataset
write_csv(c_agg, "Processed_data/agg_n92.csv")













