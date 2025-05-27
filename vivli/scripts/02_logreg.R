
source("Scripts/00_config.R")
library(tidyverse)
library(survival)
library(reshape2)

# Data
# 0 - female, 1 - male
a_imp_df <- readRDS("Processed_data/n92.rds")

# Identify reference class
# If placebo is present it is the reference
# Otherwise ac is the reference
# Accounts for trials where both are present
# Set as dpp4 for 6 trials comparing exp to exp
b_id_refs <- a_imp_df %>%
  group_by(nct_id, trt_class, trttype, arm_lvl) %>% 
  reframe() %>% 
  group_by(nct_id) %>% 
  mutate(
    has_plc = any(trttype == "plc"),
    has_ac = any(trttype == "ac"),
    ref_type = case_when(has_plc ~ "plc", !has_plc & has_ac ~ "ac"),
    ref = if_else(trttype == ref_type, 1L, 0L),
    ref = case_when(
      is.na(ref) & trt_class == "dpp4" ~ 1L,
      is.na(ref) & trt_class != "dpp4" ~ 0L,
      TRUE ~ ref
    )
  ) %>% 
  ungroup() %>% 
  select(nct_id, trt_class, trttype, arm_lvl, ref)

# Review trials with only experimental arms
# Identify active comparator
# 6 trials
# DPP-4 always the comparator or listed second in trial descriptions
# Taking as the active comparator
z_trials_exp <- a_imp_df %>% 
  left_join(b_id_refs) %>% 
  group_by(nct_id) %>% 
  filter(any(ref == 1L & trt_class == "dpp4")) %>% 
  ungroup() %>% 
  distinct(nct_id)

# Join datasets
# Assign active for dpp4 in 6 trials comparing novel antidiabetics
b_join <- a_imp_df %>% 
  left_join(b_id_refs) %>% 
  mutate(
    trttype = case_when(
      nct_id %in% z_trials_exp$nct_id & trt_class == "dpp4" ~ "ac",
      TRUE ~ trttype
    )
  )

b_join %>% count(nct_id) # 92 trials

# dpp4 as ac in 6 trials comparing novel antidiabetics
b_join %>% 
  filter(nct_id %in% z_trials_exp$nct_id) %>% 
  group_by(nct_id) %>% 
  count(trttype, trt_class)

# Save for summary statistics
saveRDS(b_join, "processed_data/n92_updated.rds")

# Numerically order treatment class for each trial
# 0 - reference, 1... = other arms
b_format_refs <- b_join %>% 
  group_by(nct_id, trt_class, ref) %>% 
  reframe() %>%
  arrange(nct_id, desc(ref)) %>% 
  group_by(nct_id) %>% 
  mutate(
    n = n(),
    ref_class = case_when(
      ref == 1L ~ paste0(trt_class, "_0"),
      ref == 0L ~ paste0(trt_class, "_", row_number() - 1)
    )
  ) %>% 
  ungroup()

# Join ordered references
# Format reference class and get reference level
b_format_df <- b_join %>%
  left_join(b_format_refs %>% select(-n)) %>% 
  select(sponsor:trt_class, ref, ref_class, everything())

# All have a reference class (i.e. _0)
b_format_df %>% 
  group_by(nct_id) %>% 
  reframe(has_ref = any(grepl("_0$", ref_class))) %>% 
  count(has_ref)

# Nest
# Relevel factors based on reference class for each nest
b_nest <- b_format_df %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  mutate(
    data = map(data, function(d) {
      ref <- d$ref_class[grepl("_0$", d$ref_class)][1]
      if (!is.na(ref)) {
        d$ref_class <- fct_relevel(d$ref_class, ref)
      } else {
        warning("No reference class found")
      }
      d$ref_class <- factor(d$ref_class, levels = c(ref, setdiff(unique(d$ref_class), ref)))
      d
    })
  )

# Fit models
# Treatment class
b_nest$cox_unadj <- map(b_nest$data,  ~ coxph(Surv(weeks, event) ~ ref_class, data = .x))
b_nest$glm_unadj <- map(b_nest$data, ~ glm(event ~ ref_class, data = .x, family = "binomial"))

# Main effects (age and sex)
b_nest$cox_adj <- map(b_nest$data,  ~ coxph(Surv(weeks, event) ~ sex + age10, data = .x))
b_nest$glm_adj <- map(b_nest$data, ~ glm(event ~ sex + age10, data = .x, family = "binomial"))

# Age and sex-treatment interactions
b_nest$cox_int <- map(b_nest$data,  ~ coxph(Surv(weeks, event) ~ ref_class * (sex + age10), data = .x))
b_nest$glm_int <- map(b_nest$data, ~ glm(event ~ ref_class * (sex + age10), data = .x, family = "binomial"))

# Gather results
c_gather_df <- b_nest %>% gather("modeltype", "model", cox_unadj:glm_int)
c_gather_df$res <- map(c_gather_df$model, ~ broom::tidy(.x))
c_gather_df$diag <- map(c_gather_df$model, ~ broom::glance(.x))
c_gather_df$vcov <- map(c_gather_df$model, ~ vcov(.x))

# Tidy model dataset
# Add reference class
# Filter out intercepts and NA estimates (where class was not present)
# Remove extreme estimates (high standard error >1000)
d_tidy_res <- c_gather_df %>% 
  select(nct_id, modeltype, res) %>% 
  unnest(res) %>% 
  left_join(
    b_id_refs %>% 
      filter(ref == 1) %>% 
      select(nct_id, ref = trt_class) %>% 
      distinct()
  ) %>% 
  filter(!is.na(estimate)) %>% 
  separate(modeltype, into = c("modeltype", "spec"), sep = "_", extra = "merge") %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
    term = gsub("ref_class", "", term),
    term = gsub("_[0-9]", "", term)
  ) %>% 
  # group_by(nct_id, modeltype, main_inter) %>%
  # filter(!any(std.error > 1000)) %>%
  # ungroup() %>% 
  rename(covariate = term) %>% 
  select(nct_id:covariate, ref, everything())

# Investigate trials returning a coefficient for the reference
# 65 trials
d_tidy_res %>% filter(covariate == ref) %>% distinct(nct_id)
# None after group-level factor ordering
d_tidy_res %>% filter(covariate == ref) %>% distinct(nct_id)

# Tidy diagnostics
d_tidy_diag <- c_gather_df %>% 
  select(nct_id, modeltype, diag) %>% 
  unnest(diag)

# Convert vcov to correlation matrices
# Keep upper triangle (removing lower and diagonals), rest as NA
d_vcov <- c_gather_df %>% 
  select(nct_id, modeltype, vcov) %>% 
  mutate(
    vcov = map(vcov, function(v) {
      cor <- cov2cor(v)
      cor[!upper.tri(cor)] <- NA_real_
      cor
    })
  )

# Convert to long-format data
# Long-format and filter out NAs (lower triangle, diagonals)
d_tidy_vcov <- d_vcov %>% 
  mutate(
    vcov_long = map(vcov, function(v) {
      df <- as.data.frame(v) %>% 
        as_tibble(rownames = "row") %>% 
        gather("col", "value", -row) %>% 
        filter(!is.na(value))
      df
    })
  ) %>% 
  select(nct_id, modeltype, vcov_long) %>% 
  unnest(vcov_long) %>% 
  ungroup()

d_tidy_vcov %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save tidy results
saveRDS(d_tidy_res, "Outputs/tidy_mdl_res.rds")
saveRDS(d_tidy_diag, "Outputs/tidy_mdl_diag.rds")
saveRDS(d_tidy_vcov, "Outputs/tidy_mdl_vcov.rds")

# Nest for plotting
d_tidy_nest <- d_tidy_res %>% group_by(spec) %>% nest()

# Plots
e_plot_trt <- d_tidy_nest$data[[1]] %>% 
  ggplot(
    aes(
      x = nct_id,
      y = estimate,
      ymin = lci,
      ymax = uci,
      colour = modeltype
    )
  ) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~covariate, scales = "free", ncol = 3) +
  theme_bw() +
  coord_flip() +
  labs(x = "ClinicalTrials.gov ID", y = "Log estimate")

ggsave(
  "Outputs/plot_treatment.png",
  e_plot_trt,
  width = 15,
  height = 10,
  units = "in"
)

# Main effects (age, sex)
e_plot_main <- d_tidy_nest$data[[2]] %>% 
  ggplot(
    aes(
      x = nct_id,
      y = estimate,
      ymin = lci,
      ymax = uci,
      colour = modeltype
    )
  ) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~covariate, scales = "free", ncol = 4) +
  theme_bw() +
  coord_flip() +
  labs(x = "ClinicalTrials.gov ID", y = "Log estimate")

ggsave(
  "Outputs/plot_age_sex.png",
  e_plot_main,
  width = 15,
  height = 10,
  units = "in"
)

# Interactions
e_plot_int <- d_tidy_nest$data[[3]] %>% 
  group_by(nct_id, modeltype, covariate) %>%
  filter(!any(std.error > 1000)) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = nct_id,
      y = estimate,
      ymin = lci,
      ymax = uci,
      colour = modeltype
    )
  ) +
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, colour = "red") +
  facet_wrap(~covariate, scales = "free", ncol = 4) +
  theme_bw() +
  labs(x = "Log estimate", y = "ClinicalTrials.gov ID") +
  coord_flip()

ggsave(
  "Outputs/plot_interactions.png",
  e_plot_int,
  width = 15,
  height = 10,
  units = "in"
)
