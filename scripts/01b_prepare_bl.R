
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import dattr data and baseline chars
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds")
a_imp_bl <- read_csv("data/base_dsp.csv")

# Extract age and male_prop from baseline chars
b_bl_ext <- a_imp_bl %>% filter(variable %in% c("age", "male", "n"))

b_n <- b_bl_ext %>% 
  filter(variable == "n") %>% 
  select(trial_id:arm_id, sample = first)

b_male <- b_bl_ext %>% 
  filter(variable == "male") %>% 
  select(trial_id:arm_id, first, second) %>% 
  rename(
    male_n = first,
    male_pcnt = second
  ) %>% 
  left_join(b_n) %>% 
  mutate(
    male_pcnt = case_when(
      is.na(male_pcnt) & !is.na(male_n) ~ male_n / sample,
      TRUE ~ male_pcnt / 100
    )
  ) %>% 
  select(trial_id, arm_id, male_prop = male_pcnt)

b_age <- b_bl_ext %>% 
  filter(variable == "age") %>% 
  mutate(
    second = case_when(
      second_format == "median" ~ NA,
      TRUE ~ second
    )
  ) %>% 
  rename(
    age_mean = first,
    age_sd = second
  ) %>% 
  select(trial_id:arm_id, age_mean, age_sd)

# Join baseline characteristics
# 380/387 trials present
b_join_bl <- b_n %>% 
  left_join(b_male) %>% 
  left_join(b_age) %>% 
  inner_join(a_imp_df %>% select(trial_id) %>% distinct())

b_join_bl %>% count(trial_id)
b_join_bl %>% reframe(across(everything(), ~ sum(is.na(.))))

# Add atc class
# Collapse by treatment group
# Weighted average for male_prop and age_mean
# Pooled standard deviation for age_sd
b_add_class <- read_csv("created_metadata/arm_id_trt.csv") %>%
  left_join(b_join_bl) %>% 
  rename(class_short = class) %>% 
  group_by(trial_id, class_short) %>% 
  reframe(
    arm_n = sum(sample),
    male_prop = sum(male_prop * sample) / sum(sample),
    age_mean = sum(age_mean * sample) / sum(sample),
    age_sd = sqrt(
      sum(
        (sample - 1) * age_sd^2 
        + sample * (age_mean - sum(age_mean * sample) / sum(sample))^2
      ) / (sum(sample) - 1)
    )
  )

b_add_class %>% count(trial_id)
b_add_class %>% reframe(across(everything(), ~ sum(is.na(.))))

# Join with dattr
c_join <- a_imp_df %>% left_join(b_add_class %>% select(-arm_n))
c_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Impute small amount of missingness (results not posted)
library(mice)
d_imp <- mice(c_join, m = 20, method = "pmm", seed = 123)
d_values <- complete(d_imp, 1) %>% as_tibble()
d_values %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save
saveRDS(d_imp, "processed_data/tidy_agg_n386_with_bl.rds")


