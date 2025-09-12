
# Prepare supplementary tables
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Supplementary figure 1 - IPD age-sex distribution ----------------------------

# Import ATC lookup
a_imp_atc <- read_csv("created_metadata/updated_class_names_codes.csv") %>% 
  mutate(
    new_class = case_match(
      new_class,
      "a_gluc" ~ "agluc",
      "biguanides" ~ "biguanide",
      .default = new_class
    )
  ) %>% 
  select(class = new_class, name = atc_short)

# Import synthetic IPD (for class-specific age10 means)
a_imp_ipd <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd) %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(readRDS("processed_data/res_n90.rds") %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Check age-sex distributions - normal, similar between sexes
sf1_plot <- a_imp_ipd %>%
  mutate(
    class = case_match(
      class,
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2"
    ),
    class = factor(class, levels = c("DPP4", "GLP-1", "SGLT2"))
  ) %>% 
  mutate(sex = if_else(sex == 0L, "Female", "Male")) %>% 
  ggplot(aes(x = age10, fill = as.factor(sex))) + 
  geom_density(colour = "black", alpha = 0.3) + 
  facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(n.breaks = 6) +
  theme_classic() +
  labs(x = "Age measured in decades", y = "Density", fill = "Sex")

sf1_plot

ggsave(
  "output/suppl/supp_fig1.png",
  sf1_plot,
  width = 8,
  height = 4,
  units = "in"
)

# Supplementary table 1 - main treatment effects for attrition -----------------

# Import obj1 model
a_imp_fit <- readRDS("output/obj1/brm_fit.rds")

# Extract posterior draws and summarise
st1_draws <- as_draws_df(a_imp_fit)
st1_draws_sum <- summarise_draws(st1_draws)

st1_tidy <- st1_draws_sum %>% 
  filter(grepl("^b_", variable) & !grepl("Intercept", variable)) %>% 
  mutate(
    class = case_match(
      variable,
      "b_classa_gluc" ~ "Alpha-glucosidase",
      "b_classbiguanides" ~ "Biguanides",
      "b_classdpp4" ~ "DPP4",
      "b_classglp1" ~ "GLP-1",
      "b_classinsulin" ~ "Insulin",
      "b_classsglt2" ~ "SGLT2",
      "b_classsulf" ~ "Sulfonylurea",
      "b_classthia" ~ "Thiazolidinedione"
    ),
    class = factor(
      class,
      levels = c(
        "DPP4", "GLP-1", "SGLT2",
        "Alpha-glucosidase", "Biguanides", "Insulin", "Sulfonylurea",
        "Thiazolidinedione"
      )
    )
  ) %>% 
  arrange(class) %>% 
  select(class, everything()) %>% 
  select(-variable, -median, -mad, -ess_tail) %>% 
  mutate(across(mean:ess_bulk, ~ round(., 3)))

st1_tidy
write_csv(st1_tidy, "output/suppl/supp_tbl1.csv")

# Supplementary figure 2 - posterior predictive check --------------------------

# Predictive check on arm-level attrition counts
sf1_ppc <- pp_check(a_imp_fit, type = "intervals")

ggsave(
  "output/suppl/supp_fig2.png",
  sf1_ppc,
  width = 10,
  height = 4,
  units = "in"
)

# Supplementary table 2 - interactions for attrition ---------------------------

a_imp_fit2 <- readRDS("output/obj2/inter_fit_n90.rds")

st2_sum <- summary(a_imp_fit2) %>%
  as_tibble() %>%
  filter(grepl("^beta|^d", parameter)) %>% 
  filter(!grepl("delta", parameter)) %>% 
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("beta|d\\[", "", parameter)
  ) %>%
  separate(
    parameter,
    into = c("inter", "term"),
    sep = "\\."
  )

st2_tidy <- st2_sum %>% 
  filter(
    !grepl("agluc|biguan|insul|sulf|thia", term),
    !grepl("A10A|A10BA|A10BB|A10BF|A10BG", inter)
  ) %>% 
  mutate(
    inter = gsub("\\:|d", "", inter),
    inter = case_match(
      inter, 
      "x2" ~ "Age", 
      "x1TRUE" ~ "Male",
      "A10BH" ~ "DPP4",
      "A10BJ" ~ "GLP-1",
      "A10BK" ~ "SGLT2",
    ),
    class = gsub("trtclass", "", term),
    class = case_match(
      class,
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2"
    ),
    class = if_else(inter %in% c("DPP4", "GLP-1", "SGLT2"), inter, class),
    class = if_else(is.na(class), "Adj", class),
    class = factor(class, levels = c("Adj", "DPP4", "GLP-1", "SGLT2")),
    across(mean:Rhat, ~ round(., 3))
  ) %>% 
  select(-Tail_ESS) %>%
  select(class, inter, mean:`2.5%`, `97.5%`, Rhat, Bulk_ESS)

st2_tidy
write_csv(st2_tidy, "output/suppl/tbl2.csv")

# Supplementary table 3 - model heterogeneity ----------------------------------

# Import models
a_imp_inter <- readRDS("output/obj2/inter_fit_n90.rds")

a_inter_tau <- summary(a_imp_inter) %>% 
  as_tibble() %>%
  filter(parameter == "tau")

# Supplementary figure 3 - Cox vs logreg ---------------------------------------

# Import tidied model results
a_imp_res <- read_csv("vivli/res_n92.csv") %>% 
  group_by(spec) %>% 
  nest()

# Generate plots
e_plot_main <- a_imp_res$data[[1]] %>% 
  filter(term != "(Intercept)") %>% 
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
  facet_wrap(~term, scales = "free", ncol = 4) +
  theme_classic() +
  coord_flip() +
  labs(x = NULL, y = "Log estimate")

e_plot_main

ggsave(
  "output/suppl/supp_fig4a.png",
  e_plot_main,
  width = 12,
  height = 10,
  units = "in"
)

e_plot_cov <- a_imp_res$data[[2]] %>% 
  filter(term != "(Intercept)") %>% 
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
  facet_wrap(~term, scales = "free", ncol = 4) +
  theme_classic() +
  coord_flip() +
  labs(x = NULL, y = "Log estimate")

e_plot_cov

ggsave(
  "output/suppl/supp_fig4b.png",
  e_plot_cov,
  width = 12,
  height = 10,
  units = "in"
)

e_plot_inter <- a_imp_res$data[[3]] %>% 
  filter(!std.error > 11 & term != "(Intercept)") %>%
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
  facet_wrap(~term, scales = "free", ncol = 4) +
  theme_classic() +
  theme(
    axis.text.y = element_blank()
  ) +
  coord_flip() +
  labs(x = NULL, y = "Log estimate")

e_plot_inter

ggsave(
  "output/suppl/supp_fig4c.png",
  e_plot_inter,
  width = 12,
  height = 10,
  units = "in"
)
