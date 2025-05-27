
# Replicating summaries from JAMA paper
# N participants
# Attrition - n (%)
# Age - mean (sd), 5-95th percentile
# Sex - n (%)
# N trials
# N comparisons within trials
# Trial duration - median, 5-95th percentile
# N participants per drug class incl. placebo

source("Scripts/00_config.R")
library(tidyverse)

# Data
# Descale age
# Get higher-level grouping for summary statistics
a_imp_df <- readRDS("Processed_data/n92_updated.rds") %>% 
  mutate(
    age = age10 * 10,
    sum_grp = case_when(
      trttype == "exp" ~ trt_class,
      trttype == "plc" ~ "placebo",
      trttype == "ac" ~ "ac"
    )
  ) %>% 
  group_by(nct_id) %>% 
  mutate(
    n_arms = length(unique(trt_class)),
    duration = round(max(weeks))
  ) %>% 
  ungroup()

# Inspect
hist(a_imp_df$age10)
table(a_imp_df$sex)
hist(a_imp_df$weeks)

# Overall summaries
sum_chars <- a_imp_df %>% 
  reframe(
    n = n(),
    n_attr = sum(event),
    pcnt_attr = round(n_attr/n * 100, 2),
    mean_age = mean(age),
    sd_age = sd(age),
    lci_age = mean_age + 1.96 * (sd_age/sqrt(n())),
    uci_age = mean_age - 1.96 * (sd_age/sqrt(n())),
    n_male = sum(sex == 1L),
    pcnt_male = round(n_male/n() * 100, 2),
    n_female = sum(sex == 0L),
    pcnt_female = round(n_female/n() * 100, 2),
    n_trials = length(unique(nct_id)),
    median_dur = median(duration),
    min_dur = min(duration),
    max_dur = max(duration)
  )

sum_arms <- a_imp_df %>% 
  count(nct_id, n_arms) %>% 
  count(n_arms) %>% 
  mutate(pcnt_arms = round(n/92 * 100, 2))

sum_classes <- a_imp_df %>% 
  count(nct_id, trt_class) %>% 
  count(trt_class) %>% 
  mutate(pcnt = round(n/92 * 100, 2))

# Save
write_csv(sum_chars, "Outputs/sum_overall_chars.csv")
write_csv(sum_arms, "Outputs/sum_overall_arms.csv")
write_csv(sum_classes, "Outputs/sum_overall_classes.csv")

# Summarise by treatment class
sum_chars_class <- a_imp_df %>%
  group_by(trt_class) %>% 
  reframe(
    n = n(),
    n_attr = sum(event),
    pcnt_attr = round(n_attr/n * 100, 2),
    mean_age = mean(age),
    sd_age = sd(age),
    lci_age = mean_age + 1.96 * (sd_age/sqrt(n())),
    uci_age = mean_age - 1.96 * (sd_age/sqrt(n())),
    n_male = sum(sex == 1L),
    pcnt_male = round(n_male/n() * 100, 2),
    n_female = sum(sex == 0L),
    pcnt_female = round(n_female/n() * 100, 2),
    n_trials = length(unique(nct_id)),
    median_dur = median(duration),
    min_dur = min(duration),
    max_dur = max(duration)
  )

sum_arms_class <- a_imp_df %>%
  count(nct_id, trt_class, n_arms) %>% 
  count(trt_class, n_arms) %>% 
  mutate(pcnt = round(n/92 * 100, 2))

sum_classes_class <- a_imp_df %>% 
  group_by(nct_id) %>% 
  filter(any(trt_class %in% "dpp4")) %>% 
  ungroup() %>% 
  count(nct_id, trt_class) %>% 
  count(trt_class) %>% 
  mutate(
    pcnt = round(n/n[trt_class == "dpp4"] * 100, 2),
    class = "dpp4"
  ) %>% 
  full_join(
    a_imp_df %>% 
      group_by(nct_id) %>% 
      filter(any(trt_class %in% "glp1")) %>% 
      ungroup() %>% 
      count(nct_id, trt_class) %>% 
      count(trt_class) %>% 
      mutate(
        pcnt = round(n/n[trt_class == "glp1"] * 100, 2),
        class = "glp1"
      )
  ) %>% 
  full_join(
    a_imp_df %>% 
      group_by(nct_id) %>% 
      filter(any(trt_class %in% "sglt2")) %>% 
      ungroup() %>% 
      count(nct_id, trt_class) %>% 
      count(trt_class) %>% 
      mutate(
        pcnt = round(n/n[trt_class == "sglt2"] * 100, 2),
        class = "sglt2"
      )
  )

# Save
write_csv(sum_chars_class, "Outputs/sum_class_chars.csv")
write_csv(sum_arms_class, "Outputs/sum_class_arms.csv")
write_csv(sum_classes_class, "Outputs/sum_class_classes.csv")