
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import aggregate data and prepared ipd outputs
a_imp_agg <- readRDS("processed_data/agg_n391.rds")
a_imp_ipd_res <- readRDS("processed_data/res_n92.rds")
a_imp_vcov <- readRDS("processed_data/vcov_n92.rds")

# Remove intercepts
# Add placebo (ref) rows to ipd results
# Add placebo interactions with age10 and sex
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
  # bind_rows(
  #   a_imp_ipd_res %>% 
  #     filter(spec == "int") %>% 
  #     transmute(
  #       nct_id,
  #       spec,
  #       class_short = "placebo:age10",
  #       atc_short = "placebo",
  #       ref,
  #       estimate = NA
  #     )
  # ) %>% 
  # bind_rows(
  #   a_imp_ipd_res %>% 
  #     filter(spec == "int") %>% 
  #     transmute(
  #       nct_id,
  #       spec,
  #       class_short = "placebo:sex",
  #       atc_short = "placebo",
  #       ref,
  #       estimate = NA
  #     )
  # ) %>% 
  arrange(nct_id, spec) %>% 
  distinct() %>% 
  filter(!(spec == "adj" & class_short == "placebo"))

# Get interaction ipd results
b_ipd_int <- b_add_plc %>% filter(spec == "int")

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
  
  est_terms <- b_ipd_int %>% 
    filter(nct_id == id & !is.na(estimate)) %>% 
    pull(class_short)
  
  matched_terms <- est_terms[est_terms %in% rownames(mat)]
  
  out <- mat[matched_terms, matched_terms, drop = FALSE]
  return(out)
  
})

# Setup networks ---------------------------------------------------------------

# Aggregate arm-level data - counts, binomially distributed
c_agg <- set_agd_arm(
  data = a_imp_agg,
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
  data = b_ipd_int,
  study = nct_id,
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
