
# This script prepares data for the analysis objectives
# AACT and IPD trial data are prepared separately and joined together

source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import datasets
# Aggregate data for 391 trials from AACT
# Aggregate data for 92 IPD trials from Vivli
# Metadata from ATC for antidiabetic classes
a_imp_aact <- read_csv("data/trials_n391_dattr.csv")
a_imp_ipd <- read_csv("vivli/agg_n92.csv")
a_imp_atc <- read_csv("created_metadata/atc.csv")

# Prepare atc lookup
b_atc_prep <- a_imp_atc %>% 
  select(class_short, atc_short = code_short) %>% 
  distinct()

# Prepare aggregate data from AACT trials --------------------------------------

# Remove trials labelled as being in IPD - added in from Vivli set
# Add medication class names by joining atc lookup
b_aact_prep <- a_imp_aact %>% 
  filter(source != "ipd") %>% 
  mutate(
    atc_short = case_when(
      atc_new == "placebo" ~ "placebo",
      atc_new != "placebo" ~ substr(atc_new, 1, 5)
    )
  ) %>% 
  left_join(b_atc_prep)

# Create metadata for class name and atc code combinations
b_class_meta <- b_aact_prep %>% 
  count(class_short, atc_short) %>% 
  mutate(
    new_class = case_when(
      atc_short == "A10BX" ~ "other",
      !is.na(class_short) ~ class_short,
      is.na(class_short) & atc_short == "placebo" ~ atc_short,
      is.na(class_short) & atc_short == "A10B" ~ "remove",
      TRUE ~ class_short
    )
  ) %>% 
  select(-n)

# write_csv(b_class_meta, "created_metadata/updated_class_names_codes.csv")

# Join metadata for class name and atc code combinations
# Remove three trials where comparators were general A10Bs, non-informative
# Add treatment type (experimental, active, placebo)
# Select core variables
b_aact_tidy <- b_aact_prep %>% 
  left_join(b_class_meta) %>% 
  group_by(trial_id) %>% 
  filter(!any(new_class == "remove")) %>% 
  ungroup() %>% 
  mutate(
    trttype = case_when(
      class_short %in% c("dpp4", "glp1", "sglt2") ~ "exp",
      class_short == "placebo" ~ "plc",
      TRUE ~ "ac"
    )
  ) %>% 
  select(
    trial_id:source, 
    class = new_class, atc = atc_short, arm_n, arm_attr = arm_n_attr
  )

# Check NAs (none)
checkNA(b_aact_tidy)

# Collapse trial counts by treatment arm
b_aact_clps <- b_aact_tidy %>% 
  group_by(trial_id, source, class, atc) %>% 
  reframe(across(arm_n:arm_attr, ~ sum(.)))

# Check data
b_aact_clps %>% count(trial_id, source) %>% count(source)
b_aact_clps %>% count(class)

# Prepare aggregate data from IPD trials ---------------------------------------

# Rename variables
# Add source (to distinguish when joined with aact trials)
# Extract atc short code
# Tidy class names (e.g. agluc -> a_gluc)
# Select core variables
c_ipd_prep <- a_imp_ipd %>% 
  mutate(source = "ipd") %>% 
  rename(
    trial_id = nct_id,
    atc = arm_lvl,
    class = trt_class,
    arm_n = n,
    arm_attr = attr
  ) %>% 
  mutate(
    atc = gsub("_d[0-9.]+$|_dcmplx", "", atc),
    atc = if_else(atc == "A10BH__07", "A10BH07", atc),
    atc = case_when(
      class == "insulin" ~ substr(atc, 1, 4),
      class == "placebo" ~ class,
      TRUE ~ substr(atc, 1, 5)
    ),
    class = case_when(
      class == "biguanide" ~ "biguanides",
      class == "agluc" ~ "a_gluc",
      TRUE ~ class
    )
  ) %>% 
  select(-c(min_dur:med_dur))

# Check NAs (none)
checkNA(c_ipd_prep)

# Collapse counts by treatment arm
c_ipd_clps <- c_ipd_prep %>% 
  group_by(trial_id, source, class, atc) %>% 
  reframe(across(arm_n:arm_attr, ~ sum(.)))

# Check data
c_ipd_clps %>% count(trial_id, source) %>% count(source)
c_ipd_clps %>% count(class)

# Join datasets ----------------------------------------------------------------

# Join datasets, no NAs introduced
d_join <- b_aact_clps %>% full_join(c_ipd_clps)

# No NAs, 388 trials (296 aact, 92 ipd), 9 treatment classes and placebo
checkNA(d_join)
d_join %>% count(trial_id, source) %>% count(source)
d_join %>% count(class)

# Check unique treatment classes per trial
# 23 trials have one unique treatment class after collapsing
d_join %>% 
  group_by(trial_id) %>% 
  reframe(n_unique = length(unique(class))) %>% 
  count(n_unique)

# Inspect treatment classes where n=1
# Comparisons of different dpp4, glp1 or sglt2; one is insulin
d_join %>% 
  group_by(trial_id) %>% 
  filter(length(unique(class)) == 1) %>% 
  ungroup() %>% 
  count(class)

# Save
saveRDS(d_join, "processed_data/tidied_agg.rds")
