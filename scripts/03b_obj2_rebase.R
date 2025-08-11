
source("scripts/00_packages.R")
source("scripts/00_config.R")
source("scripts/00_functions.R")

# Import fitted model, get tibble of posterior draws
a_imp_mdl <- readRDS("output/obj2/inter_fit_n89.rds")
a_pos <- as.data.frame(a_imp_mdl$stanfit) %>% as_tibble(rownames = "iter")

# Import synthetic IPD
a_imp_ipd <- bind_rows(readRDS("data/agg_ipd_hba1c.Rds")$ipd)

# Import atc metadata
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

# Get mean age per class from synthetic IPD
c_ipd <- a_imp_ipd %>% 
  select(nct_id, sex, age, trtcls5) %>% 
  inner_join(readRDS("processed_data/res_n90.rds") %>% distinct(nct_id)) %>% 
  rename(name = trtcls5) %>% 
  mutate(name = if_else(name == "place", "placebo", name)) %>% 
  left_join(a_imp_atc) %>% 
  filter(class %in% c("dpp4", "glp1", "sglt2")) %>% 
  mutate(age = age/10) %>% 
  rename(age10 = age)

# Get mean age per treatment class
c_ipd_age10 <- c_ipd %>% 
  group_by(class) %>% 
  reframe(mean_age = mean(age10)) %>% 
  pivot_wider(names_from = "class", values_from = "mean_age")

# Check age-sex distributions - normal, similar between sexes
c_ipd %>%
  ggplot(aes(x = age10, fill = as.factor(sex))) + 
  geom_density(colour = "black", alpha = 0.3) + 
  facet_wrap(~class, scales = "free_y") +
  scale_x_continuous(n.breaks = 6) +
  theme_classic() +
  labs(x = "Age measured in decades", y = "Density", fill = "Sex")

# Function to estimate treatment effect at class-specific mean age
get_eff_at_mean <- function(trt, trt_name) {
  d   <- a_pos[[sprintf("d[%s]", trt)]]
  g_a <- a_pos[[sprintf("beta[x2:.trtclass%s]", trt_name)]]
  g_s <- a_pos[[sprintf("beta[x1TRUE:.trtclass%s]", trt_name)]]
  d + g_a * c_ipd_age10[[trt_name]] + g_s * 0
}

summary(get_eff_at_mean("A10BK", "sglt2"))
