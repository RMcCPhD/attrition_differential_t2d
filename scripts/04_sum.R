
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import agg data (387 trials - 295 ctgov, 92 ipd)
a_imp_387 <- readRDS("processed_data/tidy_agg_n391.rds")
a_imp_92 <- read_csv("vivli/agg_n92.csv")

# Import characteristics table from nma_agesex_public
# Only for dattr trials
a_imp_char <- read_csv("data/base_dsp.csv") %>% 
  inner_join(
    a_imp_387 %>% 
      select(trial_id) %>% 
      distinct()
  )

# One trial missing - NCT04170998
a_imp_387 %>% 
  select(trial_id) %>% 
  distinct() %>% 
  anti_join(a_imp_char %>% select(trial_id) %>% distinct())

# Get relevant variables
b_char_var <- a_imp_char %>% 
  filter(!(variable %in% c("ethnicity", "race"))) %>% 
  select(trial_id, arm_id, variable, measure_type, first_format:second)

# Inspect
b_char_var %>% count(variable, measure_type)
b_char_var %>% count(variable, first_format)
b_char_var %>% count(variable, second_format)

# Handle variables separately
b_age <- b_char_var %>% 
  filter(variable == "age") %>% 
  select(trial_id, arm_id, first, second) %>% 
  rename(
    mean_age = first,
    sd_age = second
  )

b_sex <- b_char_var %>% 
  filter(variable == "male") %>% 
  select(trial_id, arm_id, first, second) %>% 
  rename(
    n_male = first,
    pcnt_male = second
  )

b_n <- b_char_var %>% 
  filter(variable == "n") %>% 
  select(trial_id, arm_id, n = first)

# Join together
b_joined <- b_n %>% 
  left_join(b_age) %>% 
  left_join(b_sex)

b_joined %>% reframe(across(everything(), ~ sum(is.na(.))))

# Label arms -------------------------------------------------------------------

# Import nma_agesex_public working agg data to get arms
# Has all of them
c_imp_nma_df <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$agg) %>% 
  select(nct_id, arm_lvl, arm_id_unq) %>% 
  rename(trial_id = nct_id) %>% 
  inner_join(a_imp_387 %>% select(trial_id) %>% distinct())

# Get ipd arms
c_imp_ipd <- read_csv("data/arm_lvl_info_as_analysed.csv") %>% 
  rename(trial_id = nct_id) %>% 
  inner_join(b_joined %>% select(trial_id) %>% distinct())

# Join together
c_arms <- c_imp_nma_df %>% 
  full_join(c_imp_ipd) %>% 
  select(-c(reference_arm:trtcls4)) %>% 
  rename(arm_id = arm_id_unq)

c_join <- b_joined %>% 
  left_join(c_arms) %>% 
  mutate()

# Get drug names for missing arm_lvl
c_fix_nas <- c_join %>%
  filter(is.na(arm_lvl)) %>% 
  left_join(
    read_csv("data/arm_data_all_cleaned.csv") %>% 
      mutate(
        arm_id_unq = case_when(
          !is.na(arm_id_subgroup) ~ arm_id_subgroup,
          TRUE ~ arm_id_unq
        )
      ) %>% 
      select(
        trial_id, arm_id = arm_id_unq, arm_id2 = arm_id, arm_label, drug_name
      )
  )

# Remove nas then join fixes
# Remove total arm levels
# Add missing arm_lvl
c_fixed_nas <- c_join %>% 
  filter(!is.na(arm_lvl)) %>% 
  full_join(c_fix_nas) %>% 
  arrange(trial_id) %>% 
  filter(!(!is.na(arm_label) & is.na(drug_name))) %>% 
  mutate(
    arm_lvl = case_when(
      is.na(arm_lvl) & grepl("tin$", drug_name) ~ "A10BH",
      is.na(arm_lvl) & grepl("zin$", drug_name) ~ "A10BK",
      is.na(arm_lvl) & grepl("insulin", drug_name) ~ "A10A",
      is.na(arm_lvl) & grepl("tide$", drug_name) ~ "A10BJ",
      is.na(arm_lvl) & grepl("metformin", drug_name) ~ "A10BA",
      is.na(arm_lvl) & grepl("azone$", drug_name) ~ "A10BG",
      is.na(arm_lvl) & drug_name == "placebo" ~ "placebo",
      TRUE ~ arm_lvl
    ),
    arm_lvl = case_when(
      arm_lvl == "placebo" ~ arm_lvl,
      TRUE ~ substr(arm_lvl, 1, 5)
    )
  )
