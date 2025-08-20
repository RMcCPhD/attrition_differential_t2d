
source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import tidied data
a_imp_df <- readRDS("processed_data/tidied_agg.rds")

# Set placebo as global reference treatment
b_set_ref <- a_imp_df %>% 
  mutate(
    class = factor(class),
    class = relevel(class, ref = "placebo"),
    arm_compl = arm_n - arm_attr,
    across(arm_n:arm_compl, as.integer)
  )

# Conditional reference for trial-level regression models
# Identify placebo or ac references
b_ref <- a_imp_df %>% 
  mutate(
    trttype = case_when(
      !class %in% c("dpp4", "glp1", "sglt2") & class == "placebo" ~ "plc",
      !class %in% c("dpp4", "glp1", "sglt2", "placebo") ~ "ac"
    ),
  ) %>% 
  group_by(trial_id) %>% 
  mutate(
    has_plc = any(trttype == "plc"),
    has_ac = any(trttype == "ac"),
    across(has_plc:has_ac, ~ if_else(is.na(.), FALSE, .)),
    ref_type = case_when(
      has_plc & !has_ac | has_plc & has_ac ~ "plc", 
      has_ac & !has_plc ~ "ac"
    ),
    ref = if_else(trttype == ref_type, 1L, NA)
  ) %>% 
  ungroup()

# Add remaining references among 44 trials (e.g. dpp4 in dpp4 vs sglt2)
b_ref_rest <- b_ref %>% 
  group_by(trial_id) %>% 
  filter(sum(is.na(ref)) == length(trial_id)) %>% 
  mutate(
    ref_type = case_when(
      is.na(ref) & class == "dpp4" & !any(class == "sglt2") ~ "dpp4",
      is.na(ref) & class == "sglt2" & !any(class == "dpp4") ~ "sglt2",
      is.na(ref) & class == "dpp4" & !any(class == "glp1") ~ "dpp4",
      TRUE ~ ref_type
    ),
    ref = if_else(class == ref_type, 1L, ref)
  )

# Join together
b_ref_all <- b_ref %>% 
  filter(!trial_id %in% b_ref_rest$trial_id) %>% 
  full_join(b_ref_rest) %>% 
  arrange(trial_id) %>% 
  mutate(ref = if_else(is.na(ref), 0L, ref))

# Check for any missing references (none)
b_ref_all %>% 
  group_by(trial_id) %>% 
  filter(sum(is.na(ref)) == length(trial_id))

# Order references
b_ref_order <- b_ref_all %>% 
  group_by(trial_id, class, ref) %>% 
  reframe() %>%
  arrange(trial_id, desc(ref)) %>% 
  group_by(trial_id) %>% 
  mutate(
    n = n(),
    ref_class = case_when(
      ref == 1L ~ paste0(class, "_0"),
      ref == 0L ~ paste0(class, "_", row_number() - 1)
    )
  ) %>% 
  ungroup()

# Join ordered references
# Format reference class and get reference level
# Manually fix reference levels for two trials with two acs
b_ref_ordered <- b_ref_all %>%
  left_join(b_ref_order %>% select(-n)) %>% 
  select(-c(trttype:ref)) %>% 
  mutate(
    ref_class = case_when(
      trial_id %in% c("NCT00676338", "NCT01147627") & grepl("thia", ref_class) ~ "thia_1",
      TRUE ~ ref_class
    )
  )

# All have a reference class (i.e. _0)
b_ref_ordered %>% 
  group_by(trial_id) %>% 
  reframe(trial_id = any(grepl("_0$", ref_class))) %>% 
  count(trial_id)

# Nest, assign factor reference levels, fit glm models
c_nest <- b_ref_ordered %>% 
  group_by(trial_id) %>% 
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
    }),
    reg = map(data, function(df) {
      glm(cbind(arm_attr, arm_n - arm_attr) ~ ref_class, data = df, family = binomial())
    })
  )

# Meta-analyse with multinma ---------------------------------------------------

# Get tidy results and VCOV
d_tidy <- c_nest %>% 
  select(trial_id, reg) %>% 
  mutate(
    res = map(reg, function(mdl){broom::tidy(mdl)}),
    vcov = map(reg, function(mdl){vcov(mdl)})
  ) %>% 
  arrange(trial_id)

# Prepare regression df
d_tidy_df <- d_tidy %>% 
  select(trial_id, res) %>% 
  unnest(res) %>% 
  mutate(across(estimate:p.value, ~ round(., 3))) %>% 
  left_join(
    c_nest %>% 
      select(trial_id, data) %>% 
      unnest(data) %>% 
      select(trial_id, ref_class) %>% 
      filter(grepl("_0$", ref_class)) %>% 
      mutate(ref_class = gsub("_0$", "", ref_class)) %>% 
      distinct()
  ) %>% 
  mutate(term = gsub("ref_class|_[0-9]", "", term))

# Add reference rows to regression df
d_add_ref <- d_tidy_df %>% 
  select(trial_id, term, estimate, ref_class) %>% 
  bind_rows(
    d_tidy_df %>% 
      group_by(trial_id) %>% 
      transmute(trial_id, term = ref_class, estimate = NA, ref_class) %>% 
      ungroup()
  ) %>% 
  distinct() %>% 
  arrange(trial_id) %>% 
  mutate(
    trt = case_when(
      term == "(Intercept)" ~ NA,
      TRUE ~ term
    ),
    trt_class = trt
  )

# Prepare vcov
d_vcov <- d_tidy$vcov
names(d_vcov) <- d_tidy$trial_id

# Construct network
e_network <- set_agd_regression(
  data = d_add_ref,
  study = trial_id,
  trt = trt,
  estimate = estimate,
  cov = d_vcov,
  regression = ~ .trt,
  trt_ref = "placebo",
  trt_class = trt_class
)

plot(e_network)

# Fit NMA model
mdl <- nma(
  e_network,
  trt_effects = "random",
  link = "identity",
  likelihood = "normal",
  class_interactions = "common",
  regression = ~ .trt,
  prior_intercept = normal(scale = 10),
  prior_trt = normal(scale = 10),
  prior_reg = normal(scale = 10),
  seed = 123,
  chains = 4, 
  cores = 4,
  warmup = 1000,
  iter = 3000,
  control = list(adapt_delta = 0.95)
)

# Summarise treatment estimates
sum_pos <- as.data.frame(mdl$stanfit) %>% 
  select(-contains("delta")) %>% 
  as_tibble(rownames = "iter") %>% 
  select(iter, starts_with("d")) %>% 
  pivot_longer(-iter) %>% 
  group_by(name) %>% 
  reframe(
    mean_main = mean(value),
    q2.5_main = quantile(value, 0.025),
    q97.5_main = quantile(value, 0.975)
  ) %>% 
  mutate(
    name = gsub("d\\[|\\]", "", name),
    or = exp(mean_main),
    q2.5_or = exp(q2.5_main),
    q97.5_or = exp(q97.5_main)
  )

# Save data for comparison
write_csv(sum_pos, "output/obj1/sum_res_agg_log.csv")

# Plot
plot_odds <- sum_pos %>% 
  ggplot(aes(x = or, xmin = q2.5_or, xmax = q97.5_or, y = fct_rev(name))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 1, colour = "red") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean odds ratio (95% credible intervals)", y = "Treatment class")

plot_odds

ggsave(
  "output/obj1/plots/plot_odds_triallvl_models_ma.png",
  plot_odds,
  width = 8,
  height = 4,
  units = "in"
)
