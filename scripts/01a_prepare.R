
source("scripts/00_config.R")
source("scripts/00_packages.R")

a_imp_agg <- read_csv("data/trials_n391_dattr.csv")
a_imp_ipd <- read_csv("vivli/agg_n92.csv")
a_imp_atc <- read_csv("created_metadata/atc.csv")

# Prepare ctgov trials
# Add treatment type
# Set A10Bs in 3 trials as "oad" for now
b_agg <- a_imp_agg %>% 
  filter(source != "ipd") %>% 
  select(trial_id, source, atc_new, n:arm_n_attr) %>% 
  rename(arm_attr = arm_n_attr) %>% 
  mutate(
    atc_short = case_when(
      atc_new == "placebo" ~ "placebo",
      atc_new != "placebo" ~ substr(atc_new, 1, 5)
    )
  ) %>% 
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
  ) %>% 
  select(trial_id, source, n:arm_attr, atc = atc_new, everything())

# Collapse treatment arms (ignoring dose) for ctgov trials
b_agg_clp <- b_agg %>% 
  group_by(trial_id, source, atc, atc_short, class_short, trttype) %>% 
  reframe(
    arm_n = sum(arm_n),
    arm_attr = sum(arm_attr)
  )

# 299 trials
b_agg_clp %>% count(trial_id, source) %>% count(source)
b_agg_clp %>% count(class_short)

# Prepare IPD
b_ipd <- a_imp_ipd %>% 
  mutate(source = "ipd") %>% 
  select(trial_id = nct_id, source, arm_lvl:attr) %>% 
  rename(
    atc = arm_lvl,
    class_short = trt_class,
    arm_n = n,
    arm_attr = attr
  ) %>% 
  mutate(
    atc = gsub("_d[0-9.]+$|_dcmplx", "", atc),
    atc = if_else(atc == "A10BH__07", "A10BH07", atc),
    atc_short = case_when(
      class_short == "insulin" ~ substr(atc, 1, 4),
      class_short == "placebo" ~ class_short,
      TRUE ~ substr(atc, 1, 5)
    ),
    class_short = case_when(
      class_short == "biguanide" ~ "biguanides",
      class_short == "agluc" ~ "a_gluc",
      TRUE ~ class_short
    )
  )

# Collapse treatment arms for ipd trials
b_ipd_clp <- b_ipd %>% 
  group_by(trial_id, source, atc, atc_short, class_short, trttype) %>% 
  reframe(
    arm_n = sum(arm_n),
    arm_attr = sum(arm_attr)
  )

# 92 trials
b_ipd_clp %>% count(trial_id, source) %>% count(source)
b_ipd_clp %>% count(class_short)

# Join datasets
# Calculate n
c_join <- b_agg_clp %>% 
  full_join(b_ipd_clp) %>% 
  group_by(trial_id) %>% 
  mutate(n = sum(arm_n)) %>% 
  ungroup()

# 92 trials, no NAs
c_join %>% count(trial_id, source) %>% count(source)
c_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Check unique treatment classes per trial
c_join %>% 
  group_by(trial_id, atc) %>% 
  reframe(n = n()) %>% 
  group_by(trial_id) %>% 
  reframe(n_atc = length(atc)) %>% 
  View()

# Remove trials with one arm (no comparator)
d_rem_1arm <- c_join %>% group_by(trial_id) %>% filter(n() != 1)

# Check unique treatment classes per trial
d_rem_1arm %>% 
  group_by(trial_id, atc, class_short, trttype) %>% 
  reframe(n = n()) %>% 
  View()

d_rem_1arm %>% 
  group_by(trial_id, trttype) %>% 
  reframe(n = n()) %>% 
  View()

# Four trials where active were multiple OADs
z_trials_oad <- c("NCT03730662", "NCT03882970", "NCT04093752")

d_rem_1arm %>% 
  group_by(trial_id, atc) %>% 
  reframe(n = n()) %>% 
  group_by(trial_id) %>% 
  reframe(n_atc = length(atc)) %>% 
  View()

# Conditionally factorise treatment class (if placebo or active comparator)
e_fmt_class <- d_rem_1arm %>% 
  group_by(trial_id, atc, class_short, trttype) %>% 
  reframe() %>% 
  group_by(trial_id) %>% 
  mutate(
    has_plc = any(trttype == "plc"),
    has_ac = any(trttype == "ac"),
    ref_type = case_when(has_plc ~ "plc", !has_plc & has_ac ~ "ac"),
    ref = if_else(trttype == ref_type, 1L, 0L),
    ref = case_when(
      # Cases where there was no placebo and multiple active comparators
      trial_id == "NCT00676338" & class_short != "biguanides" ~ 0L,
      trial_id %in% c("NCT01147627", z_trials_oad) & class_short != "insulin" ~ 0L,
      TRUE ~ ref
    )
  ) %>% 
  ungroup()

# Conditionally format treatment class (where there were only experimental arms)
# DPP4 as the reference class in these cases
e_fmt_dpp4 <- e_fmt_class %>% 
  group_by(trial_id) %>% 
  filter(is.na(ref)) %>% 
  mutate(
    has_dpp4 = any(class_short == "dpp4"),
    has_sglt2 = any(class_short == "sglt2"),
    ref_type = case_when(
      sum(class_short == "dpp4") == 1 & has_dpp4 ~ "dpp4",
      sum(class_short == "sglt2") == 1 & !has_dpp4 & has_sglt2 ~ "sglt2"
    ),
    ref = if_else(class_short == ref_type, 1L, 0L)
  ) %>% 
  ungroup()

# Deal with remaining 19 trials with experimental arms (same classes compared)
# Manually assign by checking against ctgov then import
e_fmt_rest <- e_fmt_dpp4 %>% 
  group_by(trial_id) %>% 
  filter(is.na(ref)) %>% 
  mutate(
    ref_type = case_when(atc == min(atc) ~ atc),
    ref = if_else(atc == ref_type, 1L, 0L),
    ref = if_else(is.na(ref), 0L, ref)
  ) %>% 
  fill(ref_type, .direction = "downup") %>% 
  ungroup()

# Combine conditional references
e_fmt_all <- e_fmt_class %>% 
  filter(!is.na(ref)) %>% 
  full_join(e_fmt_dpp4 %>% filter(!is.na(ref))) %>% 
  full_join(e_fmt_rest %>% filter(!is.na(ref))) %>% 
  select(-c(contains("has")))

# Join references
e_fmt_join <- d_rem_1arm %>% left_join(e_fmt_all)

# Numerically order treatment arms
# 0 = reference, 1... = other arms
e_order_refs <- e_fmt_join %>% 
  group_by(trial_id, atc, trttype, class_short, ref_type, ref) %>% 
  reframe() %>% 
  group_by(trial_id) %>% 
  mutate(ordering = row_number() - 1) %>% 
  arrange(trial_id, ordering) %>% 
  mutate(ref_class = paste0(class_short, "_", ordering)) %>% 
  ungroup()

# Join ordered treatment arms/classes with data
e_order_join <- e_fmt_join %>% 
  left_join(e_order_refs) %>% 
  select(trial_id:trttype, ref_class, arm_n:arm_attr) %>% 
  ungroup()

# 387 trials (295 ctgov, 92 ipd)
e_order_join %>% count(trial_id, source) %>% count(source)

# No NAs
e_order_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save
saveRDS(e_order_join, "processed_data/tidy_agg_n391.rds")
