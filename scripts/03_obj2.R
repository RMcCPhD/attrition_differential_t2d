
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Import IPD model exports (coefficients, variance-covariance)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Get sample sizes
b_samples <- a_imp_df %>%
  select(trial_id, class_short, arm_n) %>% 
  rename(nct_id = trial_id) %>% 
  filter(nct_id %in% a_imp_ipd_res$nct_id)

# Get logistic reg outputs
b_log <- a_imp_ipd_res %>% 
  filter(modeltype == "glm") %>% 
  filter(std.error < 11)

# Test for unadjusted models
b_unadj <- b_log %>% filter(spec == "unadj" & term != "(Intercept)")

b_unadj_vcov <- a_imp_vcov %>% filter(modeltype == "glm_unadj")

# Add placebo ref
b_add_ref <- b_unadj %>% 
  bind_rows(
    b_unadj %>% 
      group_by(nct_id) %>% 
      filter(n() == 1) %>% 
      ungroup() %>% 
      transmute(
        nct_id,
        modeltype,
        spec,
        term = "placebo",
        ref = NA,
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        lci = NA,
        uci = NA
      )
  ) %>% 
  arrange(nct_id) %>% 
  distinct() %>% 
  select(nct_id, term, estimate, std.error, ref) %>% 
  rename(se = std.error)

# Add intercept (placebo est/se) where >1 arm
# Add sample size
b_int_targets <- b_add_ref %>% 
  group_by(nct_id) %>% 
  filter(sum(is.na(estimate)) == 0) %>% 
  ungroup()

b_int_plc <- b_log %>% 
  filter(spec == "unadj") %>% 
  inner_join(b_int_targets %>% select(nct_id) %>% distinct()) %>% 
  mutate(term = if_else(grepl("Inter", term), "placebo", term)) %>% 
  filter(term == "placebo") %>% 
  select(nct_id, term, estimate, std.error, ref) %>% 
  rename(se = std.error)

b_all_plc <- b_add_ref %>% 
  full_join(b_int_plc) %>% 
  arrange(nct_id) %>% 
  mutate(estimate = if_else(term == "placebo", NA, estimate)) %>% 
  select(-ref)

b_all_plc %>% group_by(nct_id) %>% filter(n() > 2) %>% View()
  
# Test for treatment-only
b_network <- set_agd_contrast(
  data = b_all_plc,
  study = nct_id,
  trt = term,
  y = estimate,
  se = se,
  trt_ref = "placebo"
)

# Plot network
plot(b_network, level = "treatment")

# Meta-regression
b_mdl <- nma(
  b_network,
  trt_effects = "random",
  prior_trt = normal(0, 100),
  prior_het = half_normal(5)
)

# Quick plot of estimates
b_plot <- as.data.frame(summary(b_mdl))[1:7,] %>% 
  mutate(
    parameter = gsub("d\\[", "", parameter),
    parameter = gsub("\\]", "", parameter)
  ) %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(parameter))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean log-odds (95% credible intervals)", y = "Treatment class")

b_plot

ggsave(
  "output/obj2/plot_nma_log_odds.png",
  b_plot,
  width = 8,
  height = 4,
  units = "in"
)
