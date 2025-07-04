
source("00_config.R")
source("00_packages.R")

# Import aggregate data and prepared ipd outputs
a_imp_agg <- readRDS("data/agg_n387.rds")
a_imp_ipd_res <- readRDS("data/res_n92.rds")
a_imp_vcov <- readRDS("data/vcov_n92.rds")

# Remove intercepts
# Add placebo (ref) rows to ipd results
# Add placebo interactions with age10 and sex
# Get interaction regression results
b_add_plc <- a_imp_ipd_res %>% 
  filter(!grepl("Inter", class_short)) %>% 
  bind_rows(
    a_imp_ipd_res %>% 
      group_by(nct_id, spec) %>% 
      transmute(
        nct_id,
        spec,
        class_short = "placebo",
        atc_short = "placebo",
        ref,
        estimate = NA
      )
  ) %>% 
  rename(trial_id = nct_id) %>% 
  arrange(trial_id, spec) %>% 
  distinct() %>% 
  filter(
    !(spec == "adj" & class_short == "placebo"),
    spec == "int"
  )

# Add sex and age_mean
b_add_bl <- b_add_plc %>% 
  left_join(
    a_imp_agg %>% 
      select(trial_id, atc_short, sex, age_mean) %>% 
      filter(trial_id %in% b_add_plc$trial_id)
  ) %>% 
  rename(age10 = age_mean) %>% 
  distinct() %>% 
  mutate(sex = 1 - sex)

# Additional row where there were two dpp4 arms, adds age_mean
b_add_plc %>% 
  group_by(trial_id) %>% 
  reframe(n = n()) %>% 
  left_join(b_add_bl %>% group_by(trial_id) %>% reframe(n2 = n())) %>% 
  mutate(flag = if_else(n != n2, 1, 0)) %>% 
  filter(flag == 1)

# Remove NCT02597049 int, extreme covariance values but reasonable estimates
b_rm_n1 <- b_add_bl %>% filter(trial_id != "NCT02597049")

# Prepare variables for network
b_prep <- b_rm_n1 %>% 
  mutate(
    int_terms = class_short,
    class_short = case_when(
      grepl("age10|sex|\\:", class_short) ~ NA,
      TRUE ~ class_short
    )
  )

# Prepare vcov, setting NA for intercepts
b_vcov <- a_imp_vcov %>% 
  filter(modeltype == "glm_int") %>% 
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
  
  est_terms <- b_add_bl %>% 
    filter(trial_id == id & !is.na(estimate)) %>% 
    pull(class_short)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

# Remove NCT02597049
b_vcov_prep <- b_vcov_prep[-90]

# Get age-sex correlation
b_bl_cor <- read_csv("vivli/coef_cor_n92.csv") %>% 
  rename(trial_id = nct_id) %>% 
  inner_join(b_prep %>% select(trial_id) %>% distinct()) %>% 
  filter(
    modeltype == "glm_int", 
    row %in% c("sex", "age10")
  )

mean(b_bl_cor$value) # 0.007


b_cor_mat <- b_bl_cor %>%
  group_by(trial_id) %>%
  nest() %>%
  mutate(
    cor_matrix = map(data, ~{
      df <- .x
      vars <- union(df$row, df$col)
      mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
                    dimnames = list(vars, vars))
      for(i in seq_len(nrow(df))) {
        mat[df$row[i], df$col[i]] <- df$value[i]
      }
      diag(mat) <- ifelse(is.na(diag(mat)), 1, diag(mat))
      mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
      mat
    })
  )

# Setup networks ---------------------------------------------------------------

# Aggregate arm-level data - counts, binomially distributed
c_agg <- set_agd_arm(
  data = a_imp_agg %>% filter(source != "ipd"),
  study = trial_id,
  trt = class_short,
  trt_ref = "placebo",
  r = arm_attr,
  n = arm_n,
  trt_class = atc_short
)

plot(c_agg)

# IPD coefficients and variance - log-odds, multivariate normal
c_ipd <- set_agd_regression(
  data = b_prep,
  study = trial_id,
  trt = class_short,
  trt_ref = "placebo",
  estimate = estimate,
  cov = b_vcov_prep,
  regression = ~ .trt * (age10 + sex),
  trt_class = atc_short
)

plot(c_ipd)

# Combine
c_cmbn <- combine_network(c_agg, c_ipd)
plot(c_cmbn)
summary(c_cmbn)

# Simulate pseudo IPD from aggregate age and sex to get correlation
a_imp_agg %>% 
  filter(source != "ipd") %>% 
  reframe(
    mean = mean(age_mean), # 5.73
    sd = sd(age_mean) # 0.43
  )

mean(a_imp_agg$sex) # 0.45

c_pseudo <- data.frame(
  age10 = rnorm(10000, mean = 5.73, sd = 0.43),
  sex = rbinom(10000, size = 1, prob = 0.45)
)

c_pseudo_cor <- cor(c_pseudo)

# Add integration
c_integ <- add_integration(
  c_cmbn,
  age10 = distr(qnorm, mean = age_mean, sd = age_sd),
  sex = distr(qbern, prob = sex),
  cor = c_pseudo_cor,
  n_int = 64
)

# Fit model --------------------------------------------------------------------

d_mdl <- nma(
  c_cmbn,
  trt_effects = "random",
  regression = ~ .trt * (age10 + sex),
  class_interactions = "common",
  prior_intercept = normal(0, 5),
  prior_trt = normal(0, 2.5),
  prior_het = half_normal(1),
  prior_reg = normal(0, 2.5),
  chains = 4, 
  cores = 4,
  control = list(max_treedepth = 15)
)

plot_prior_posterior(d_mdl, prior = c("trt"))
dic(d_mdl)

# Extract fixed effects
e_ext <- summary(d_mdl$stanfit, pars = "d")$summary

# Plot
plot_int <- as.data.frame(e_ext) %>% 
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
  "log_odds.png",
  plot_int,
  width = 8,
  height = 4,
  units = "in"
)
