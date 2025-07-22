
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

# One trial missing - NCT04170998 (not missing, not removing in case I forgot the context)
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
  ) %>% 
  select(-c(arm_id2:drug_name))

c_fixed_nas %>% reframe(across(everything(), ~ sum(is.na(.))))

# Add treatment classes
d_imp_atc <- read_csv("created_metadata/atc.csv")

d_class <- c_fixed_nas %>% 
  rename(atc = arm_lvl) %>% 
  mutate(
    class = atc,
    class = case_match(
      class,
      "A10A" ~ "insulin",
      "A10B" ~ "oad",
      "A10BA" ~ "biguanides",
      "A10BB" ~ "sulf",
      "A10BF" ~ "a_gluc",
      "A10BG" ~ "thia",
      "A10BH" ~ "dpp4",
      "A10BJ" ~ "glp1",
      "A10BK" ~ "sglt2",
      "placebo" ~ "placebo"
    )
  )

# Save arm lookup (used here and attrisk)
write_csv(
  d_class %>% select(trial_id:arm_id, atc, class),
  "created_metadata/arm_id_trt.csv"
)

# Add events -------------------------------------------------------------------

# Join aggregate and ipd events
# Tidy atc codes
e_events <- a_imp_387 %>% 
  filter(source == "agg") %>%
  select(trial_id, atc_short, arm_n, arm_attr) %>% 
  rename(atc = atc_short) %>% 
  full_join(
    a_imp_92 %>% 
      select(nct_id, arm_lvl, n, attr) %>% 
      rename(
        trial_id = nct_id,
        atc = arm_lvl,
        arm_n = n,
        arm_attr = attr
      )
  ) %>% 
  mutate(
    atc = case_when(
      atc == "placebo" ~ atc,
      TRUE ~ substr(atc, 1, 5)
    )
  )

# Collapse treatment class characteristics for joining events
# Compute weighted mean age and total number of males per arm
e_clp <- d_class %>% 
  group_by(trial_id, atc, class) %>% 
  mutate(
    wgt_mean_age = mean_age * n,
    n_male = case_when(
      !is.na(n_male) ~ n_male,
      TRUE ~ n * (pcnt_male / 100)
    )
  ) %>% 
  reframe(
    n = sum(n),
    mean_a %>% ge = sum(wgt_mean_age) / n,
    n_male = round(sum(n_male))
  )

# Add pooled standard deviation
e_clp_add_sd <- e_clp %>% 
  left_join(
    d_class %>% 
      group_by(trial_id, atc, class) %>% 
      reframe(
        num = sum((n - 1) * sd_age ^ 2),
        denom = sum(n) - n(),
        sd_age = sqrt(num / denom)
      )
  ) %>% 
  select(trial_id:mean_age, sd_age, everything())

# Join events
# Remove arms with no attrition counts
f_add_events <- e_clp_add_sd %>% 
  left_join(
    e_events %>% 
      group_by(trial_id, atc) %>% 
      reframe(
        arm_n = sum(arm_n),
        arm_attr = sum(arm_attr)
      )
  ) %>% 
  select(-c(num, denom, arm_n)) %>% 
  rename(attr = arm_attr)

# 17 missing mean age
# 19 missing sd age
# 15 missing male count
# Look for these as manual inputs or use na.rm if not obtainable
f_add_events %>% reframe(across(everything(), ~ sum(is.na(.))))

# Distinguish source
g_add_source <- f_add_events %>% 
  mutate(
    source = case_when(
      trial_id %in% a_imp_92$nct_id ~ "ipd",
      TRUE ~ "agg"
    )
  )

# Save 
write_csv(g_add_source, "output/sum/sum_data.csv")
