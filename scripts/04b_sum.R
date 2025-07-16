
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepped summary data
a_imp_sum <- read_csv("output/sum/sum_data.csv")

# Overall ----------------------------------------------------------------------

# Total number of participants
sum(a_imp_sum$n)

# Number of trials
a_imp_sum %>% count(trial_id, source) %>% count(source)

# Proportion female
a_imp_sum %>% 
  reframe(
    n_female = sum(n) - sum(n_male, na.rm = TRUE),
    pcnt_female = 1 - sum(n_male, na.rm = TRUE) / sum(n)
  )

# Weighted mean age
a_imp_sum %>%
  mutate(wgt_mean_age = mean_age * n) %>% 
  reframe(mean_age = sum(wgt_mean_age, na.rm = TRUE) / sum(n))

# Pooled sd age
a_imp_sum %>% 
  filter(!is.na(sd_age)) %>% 
  reframe(
    num = sum((n - 1) * sd_age ^ 2),
    denom = sum(n) - n(),
    sd_age = sqrt(num / denom)
  )

# Attrition
a_imp_sum %>% 
  reframe(
    sum_attr = sum(attr),
    pct_attr = sum_attr / sum(n)
  )

# Sources ----------------------------------------------------------------------

# Total number of participants
a_imp_sum %>% group_by(source) %>% reframe(n = sum(n))

# Proportion female
a_imp_sum %>% 
  group_by(source) %>% 
  reframe(
    n_female = sum(n) - sum(n_male, na.rm = TRUE),
    pcnt_female = 1 - sum(n_male, na.rm = TRUE) / sum(n)
  )

# Weighted mean age
a_imp_sum %>% 
  group_by(source) %>% 
  mutate(wgt_mean_age = mean_age * n) %>% 
  reframe(mean_age = sum(wgt_mean_age, na.rm = TRUE) / sum(n))

# Pooled sd age
a_imp_sum %>% 
  filter(!is.na(sd_age)) %>% 
  group_by(source) %>% 
  reframe(
    num = sum((n - 1) * sd_age ^ 2),
    denom = sum(n) - n(),
    sd_age = sqrt(num / denom)
  )

# Attrition
a_imp_sum %>% 
  group_by(source) %>% 
  reframe(
    sum_attr = sum(attr),
    pct_attr = sum_attr / sum(n)
  )

# Number of comparisons (1 is where a dpp4 was compared to another, is 2)
a_imp_sum %>% 
  group_by(trial_id, source) %>% 
  reframe(
    n_comps = length(unique(class))
  ) %>% 
  count(source, n_comps) %>% 
  mutate(n_comps = if_else(n_comps == 1L, 2L, n_comps)) %>% 
  group_by(source, n_comps) %>% 
  reframe(n = sum(n))

# Number of participants per drug class
a_imp_sum %>% 
  group_by(trial_id, source, class) %>% 
  reframe(n = n()) %>% 
  count(source, class) %>% 
  arrange(source, desc(n))

# Total duration
a_imp_dur <- read_csv("data/trial_duration.csv")

a_dur <- a_imp_sum %>% 
  select(trial_id, source) %>% 
  distinct() %>% 
  left_join(a_imp_dur %>% mutate(duration = round(duration)))

a_dur %>% 
  group_by(source) %>% 
  reframe(
    med_dur = median(duration),
    min_dur = min(duration),
    max_dur = max(duration)
  )

# Source and drug class --------------------------------------------------------

# Identify main class
b_class <- a_imp_sum %>% 
  group_by(trial_id) %>% 
  mutate(
    main_class = case_when(
      class == "glp1" ~ "glp1",
      class == "dpp4" ~ "dpp4",
      class == "sglt2" ~ "sglt2"
    )
  ) %>% 
  fill(main_class, .direction = "downup") %>% 
  ungroup()

# Deal with trials where >1 class was tagged as main
b_class_mt1 <- b_class %>% 
  group_by(trial_id) %>% 
  filter(length(unique(main_class)) > 1) %>% 
  ungroup()

# write_csv(b_class_mt1, "scratch_data/main_class_mt1.csv")

# Join fixes for main treatment class (checked against ctgov for 28 trials)
b_imp_fix <- read_csv("scratch_data/fixed_main_class_mt1.csv")

b_class_fixed <- b_class %>% 
  filter(!(trial_id %in% b_imp_fix$trial_id)) %>% 
  full_join(b_imp_fix) %>% 
  arrange(trial_id)

# Total number of participants
b_class_fixed %>% group_by(source, main_class) %>% reframe(n = sum(n))

# Number of trials
b_class_fixed %>% count(trial_id, source, main_class) %>% count(source, main_class)

# Proportion female
b_class_fixed %>% 
  group_by(source, main_class) %>% 
  reframe(
    n_female = sum(n) - sum(n_male, na.rm = TRUE),
    pcnt_female = 1 - sum(n_male, na.rm = TRUE) / sum(n)
  )

# Weighted mean age
b_class_fixed %>% 
  group_by(source, main_class) %>% 
  mutate(wgt_mean_age = mean_age * n) %>% 
  reframe(mean_age = sum(wgt_mean_age, na.rm = TRUE) / sum(n))

# Pooled sd age
b_class_fixed %>% 
  filter(!is.na(sd_age)) %>% 
  group_by(source, main_class) %>% 
  reframe(
    num = sum((n - 1) * sd_age ^ 2),
    denom = sum(n) - n(),
    sd_age = sqrt(num / denom)
  )

# Attrition
b_class_fixed %>% 
  group_by(source, main_class) %>% 
  reframe(
    sum_attr = sum(attr),
    pct_attr = sum_attr / sum(n)
  )

# Number of comparisons (1 is where a dpp4 was compared to another, is 2)
b_class_fixed %>% 
  group_by(trial_id, source, main_class) %>% 
  reframe(n_comps = length(unique(class))) %>% 
  count(source, main_class, n_comps) %>% 
  mutate(n_comps = if_else(n_comps == 1L, 2L, n_comps)) %>% 
  group_by(source, n_comps, main_class) %>% 
  reframe(n = sum(n)) %>% 
  arrange(source, main_class)

# Duration
a_dur %>% 
  left_join(b_class_fixed) %>% 
  group_by(source, main_class) %>% 
  reframe(
    med_dur = median(duration),
    min_dur = min(duration),
    max_dur = max(duration)
  )

# Drug classes across trials
b_class_fixed %>% 
  group_by(trial_id, source, main_class, class) %>% 
  reframe(n = n()) %>% 
  count(source, main_class, class) %>% 
  arrange(source, main_class, desc(n)) %>% 
  print(n = 42)

# Attrition --------------------------------------------------------------------


