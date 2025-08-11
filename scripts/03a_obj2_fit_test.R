
library(dplyr)
library(stringr)
library(forcats)
library(purrr)
library(tidyr)

# --- Load ---------------------------------------------------------------------
a_imp_res   <- readRDS("processed_data/res_n90.rds")
a_imp_vcov  <- readRDS("processed_data/vcov_n90.rds")
a_imp_atc   <- readr::read_csv("created_metadata/atc.csv")

# --- Add reference rows, fix ATC, set covariate coding ------------------------
b_add_ref <- a_imp_res %>%
  bind_rows(
    a_imp_res %>%
      group_by(nct_id) %>%
      transmute(nct_id, class = ref, atc = ref, estimate = NA_real_, ref) %>%
      ungroup()
  ) %>%
  distinct() %>%
  arrange(nct_id)

# map class/regex → ATC
map_atc <- c(
  insulin   = "A10A",
  biguanide = "A10BA",
  sulf      = "A10BB",
  agluc     = "A10BF",
  thia      = "A10BG",
  dpp4      = "A10BH",
  glp1      = "A10BJ",
  sglt2     = "A10BK",
  placebo   = "placebo"
)

b_add_covs <- b_add_ref %>%
  mutate(
    atc = case_when(
      class %in% names(map_atc) ~ map_atc[class],
      str_detect(term, "insulin:")   ~ "A10A",
      str_detect(term, "biguanide:") ~ "A10BA",
      str_detect(term, "sulf:")      ~ "A10BB",
      str_detect(term, "agluc:")     ~ "A10BF",
      str_detect(term, "thia:")      ~ "A10BG",
      str_detect(term, "dpp4:")      ~ "A10BH",
      str_detect(term, "glp1:")      ~ "A10BJ",
      str_detect(term, "sglt2:")     ~ "A10BK",
      TRUE ~ atc
    ),
    # covariate indicators: NA means "not part of this coefficient"
    sex  = case_when(str_detect(term, "sex")   ~ TRUE,  class == ref ~ FALSE, TRUE ~ NA),
    age10= case_when(str_detect(term, "age10") ~ 1,     class == ref ~ 0,     TRUE ~ NA)
  )

# --- Canonical term key (this is what we align the covariance to) -------------
canon_term <- function(atc, term) {
  # study intercept stays "(Intercept)"
  if (term == "(Intercept)") return("(Intercept)")
  # prognostic main effects
  if (term == "age10") return("age10")
  if (str_detect(term, "^sex")) return("sex")
  # treatment main effects are ATC codes
  if (!is.na(atc) && !str_detect(term, ":")) return(atc)
  # interactions: use ATC:age10 or ATC:sex
  if (str_detect(term, ":age10")) return(paste0(atc, ":age10"))
  if (str_detect(term, ":sex"))   return(paste0(atc, ":sex"))
  # fallback: original
  term
}

b_ready <- b_add_covs %>%
  mutate(
    trt = atc,
    x1  = sex,
    x2  = age10,
    key = map2_chr(atc, term, canon_term)
  ) %>%
  select(nct_id, estimate, std.error, trt, x1, x2, class, term, ref, key)

# --- Remove the outlier trial (as you had) ------------------------------------
drop_id <- "NCT02597049"
b_ready_rm <- filter(b_ready, nct_id != drop_id)

# bring in the covariance list and name by nct_id (ensure exact names)
V_list <- a_imp_vcov$cov_matrix
names(V_list) <- unique(b_ready$nct_id)
V_list <- V_list[names(V_list) != drop_id]

# --- Per-study ordering, covariance reindexing, and validation ----------------
reindex_covariance <- function(df_study, V) {
  # vector of non-NA coefficients in THIS study, in the row order we will pass to multinma
  df_nonNA <- df_study %>% filter(!is.na(estimate))
  keys     <- df_nonNA$key
  
  # If V has dimnames, try to reorder by names
  if (!is.null(rownames(V)) && !is.null(colnames(V))) {
    # Are names a permutation match?
    if (setequal(keys, rownames(V)) && setequal(keys, colnames(V))) {
      V <- V[keys, keys, drop = FALSE]
    } else {
      warning("Covariance names do not match keys in study ", df_study$nct_id[1],
              " — attempting positional alignment; please verify upstream.")
      # fall through to positional subset below
    }
  }
  # If no names or mismatch, subset by position assuming V is already in the same order
  if (nrow(V) != nrow(df_nonNA)) {
    stop("Dimension mismatch in study ", df_study$nct_id[1],
         ": V has ", nrow(V), " rows, estimates have ", nrow(df_nonNA))
  }
  
  # Diagonal ≈ std.error^2 check
  se <- df_nonNA$std.error
  if (all(is.finite(se))) {
    d <- sqrt(diag(V))
    if (max(abs(d - se), na.rm = TRUE) > 5e-2) {
      warning("SE vs sqrt(diag(V)) differs in study ", df_study$nct_id[1],
              " (max abs diff = ", signif(max(abs(d - se), na.rm = TRUE), 3), ")")
    }
  }
  list(df = df_study, V = V, keys = keys)
}

# apply per study
split_df <- split(b_ready_rm, b_ready_rm$nct_id)
aligned <- map(names(split_df), function(id) {
  if (is.null(V_list[[id]])) stop("Missing covariance for study ", id)
  reindex_covariance(split_df[[id]], V_list[[id]])
})

# rebuild aligned objects
b_ready_aligned <- bind_rows(map(aligned, "df"))
V_aligned <- set_names(map(aligned, "V"), names(split_df))

# --- Hard validation: x2 steps & interaction rows -----------------------------
# All age interactions should have x2 == 1 (per decade)
age_inter_steps <- b_ready_aligned %>%
  filter(str_detect(term, ":age10", ignore.case = TRUE)) %>%
  count(class, x2, name = "n")
print(age_inter_steps)
if (any(is.na(age_inter_steps$x2) | age_inter_steps$x2 != 1)) {
  warning("Found age interaction rows with x2 != 1 (per decade) — fix before NMA.")
}

# Reference rows (estimate NA) should have x1==FALSE, x2==0
ref_ok <- b_ready_aligned %>%
  filter(is.na(estimate)) %>%
  summarise(all_ok = all(x1 == FALSE & x2 == 0, na.rm = TRUE)) %>% pull(all_ok)
if (!ref_ok) warning("Some reference rows do not have x1=FALSE and x2=0.")