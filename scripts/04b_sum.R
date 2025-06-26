
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepped summary data
a_imp_sum <- read_csv("output/sum/sum_data.csv")

# Overall ----------------------------------------------------------------------

# Total number of participants
sum(a_imp_sum$n)

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

# Overview
# 211,591 participants
# 97,998 (46.3%) female
# Weighted mean, pooled SD: 55.7 (10.2)
# 14.1% (29.751) attrition

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

# Overview
# Agg: 155,463, IPD: 56,128 participants
# Agg: 71,620 (46.1%), IPD: 26,378 (47%) female
# Agg: 55.6 (10.3), IPD: 56 (9.7) age
# Agg: 12.7% (19,688), IPD: 17.9% (10)

