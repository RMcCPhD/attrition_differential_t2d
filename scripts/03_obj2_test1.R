
# Testing setup for multinma
# Unadjusted treatment
# Combination of aggregated data and ipd outputs
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import prepared aggregate data
a_imp_df <- readRDS("processed_data/tidy_agg_n391.rds") %>% 
  mutate(
    class_short = case_match(
      class_short,
      "a_gluc" ~ "agluc",
      "biguanides" ~ "biguanide",
      .default = class_short
    ),
    atc = factor(atc, levels = c("placebo", setdiff(unique(atc), "placebo"))),
    class_short = factor(class_short),
    class_short = relevel(class_short, ref = "placebo")
  )

# Import IPD model exports (coefficients, variance-covariance)
a_imp_ipd_res <- read_csv("vivli/res_n92.csv")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Unadjusted treatment (test) --------------------------------------------------

# Setup aggregate data
b_agg_df <- a_imp_df %>% 
  filter(source != "ipd") %>% 
  select(trial_id, class_short, atc_short, ref_class:arm_attr) %>% 
  mutate(
    atc_short = if_else(class_short == "oad", "A10BX", atc_short)
  )

# Setup ipd data
b_ipd_filter <- a_imp_ipd_res %>% 
  filter(
    spec == "unadj",
    modeltype == "glm"
  ) %>% 
  rename(class_short = term) %>% 
  select(-c(std.error:uci)) %>% 
  left_join(
    b_agg_df %>% 
      select(class_short, atc_short) %>% 
      distinct()
  )

# Add placebo rows with NA estimate for reference level
b_ipd_df <- b_ipd_filter %>%
  filter(!grepl("Inter", class_short)) %>% 
  bind_rows(
    b_ipd_filter %>% 
      group_by(nct_id) %>% 
      transmute(
        nct_id,
        modeltype,
        spec,
        class_short = "placebo",
        atc_short = "placebo",
        ref,
        estimate = NA
      )
  ) %>% 
  arrange(nct_id) %>% 
  distinct()

# Prepare vcov, setting NA for intercepts
b_vcov <- a_imp_vcov %>% 
  filter(modeltype == "glm_unadj") %>% 
  mutate(
    cov2 = map(cov_matrix, function(mat) {
      
      rownames(mat) <- gsub("Intercept", "placebo", rownames(mat))
      rownames(mat) <- gsub("\\(|\\)", "", rownames(mat))
      colnames(mat) <- gsub("Intercept", "placebo", colnames(mat))
      colnames(mat) <- gsub("\\(|\\)", "", colnames(mat))
      
      placebo_cell <- which(rownames(mat) == "placebo")
      
      if (length(placebo_cell) > 0) {
        mat[placebo_cell, ] <- NA
        mat[, placebo_cell] <- NA
      }
      
      return(mat)
      
    })
  )

b_vcov_ready <- list(b_vcov$cov2)
b_vcov_ready <- b_vcov_ready[[1]]
names(b_vcov_ready) <- unique(b_vcov$nct_id)

# Keep non-NA cells
b_vcov_prep <- map2(b_vcov_ready, names(b_vcov_ready), function(mat, id) {
  
  est_terms <- b_ipd_df %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(class_short)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

# Setup networks
c_agg <- set_agd_arm(
  data = b_agg_df,
  study = trial_id,
  trt = class_short,
  trt_ref = "placebo",
  r = arm_attr,
  n = arm_n,
  trt_class = atc_short
)
  
plot(c_agg)
  
c_ipd <- set_agd_regression(
  data = b_ipd_df,
  study = nct_id,
  trt = class_short,
  trt_ref = "placebo",
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt,
  trt_class = atc_short
)

plot(c_ipd)

# Combine
c_cmbn <- combine_network(c_agg, c_ipd)
plot(c_cmbn)
summary(c_cmbn)

# Add integration
d_mdl <- nma(
  c_cmbn,
  trt_effects = "random",
  regression = ~ .trt,
  class_interactions = "common",
  prior_intercept = normal(0, 5),
  prior_trt = normal(0, 2.5),
  prior_het = half_normal(1),
  prior_reg = normal(0, 2.5),
  chains = 4, 
  cores = 4
)

plot_prior_posterior(d_mdl, prior = c("trt"))
dic(d_mdl)

# Extract fixed effects
e_ext <- summary(d_mdl$stanfit, pars = "d")$summary

# Plot
plot_unadj <- as.data.frame(e_ext) %>% 
  rownames_to_column(var = "trt_class") %>% 
  mutate(
    trt_class = gsub("d\\[|\\]", "", trt_class)
  ) %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(trt_class))) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean log-odds (95% credible intervals)", y = "Treatment class")

ggsave(
  "output/obj2/test_unadj_trt_log_odds.png",
  plot_unadj,
  width = 8,
  height = 4,
  units = "in"
)
