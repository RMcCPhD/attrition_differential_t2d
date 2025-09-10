
source("scripts/00_config.R")
source("scripts/00_packages.R")

# Import aggregate data for IPD trials
# Get rate
# Add experimental
a_imp_agg <- readRDS("processed_data/tidied_agg.rds") %>% 
  filter(source == "ipd") %>% 
  mutate(attr_rate = arm_attr/arm_n) %>% 
  group_by(trial_id) %>% 
  mutate(
    exp = case_when(
      class %in% c("dpp4", "glp1", "sglt2") ~ class
    )
  ) %>%
  fill(exp, .direction = "downup") %>% 
  ungroup()

# Manually assign primary class for 23 comparing two newer antidiabs
a_imp_agg %>% 
  group_by(trial_id, exp) %>% 
  reframe(n = n()) %>% 
  group_by(trial_id) %>% 
  filter(length(unique(exp)) != 1) %>% 
  distinct(trial_id) %>%
  inner_join(a_imp_agg) %>% 
  write_csv("scratch_data/manual_assign_exp.csv")

# Join to rest
b_imp_fix <- read_csv("scratch_data/manual_assigned_exp.csv")

b_agg <- a_imp_agg %>% 
  filter(!trial_id %in% b_imp_fix$trial_id) %>% 
  full_join(b_imp_fix) %>% 
  arrange(trial_id) %>% 
  mutate(
    trttype = case_when(
      class %in% c("dpp4", "glp1", "sglt2") ~ class,
      class == "placebo" ~ class,
      TRUE ~ "classic"
    )
  )

# Plot
plot_rates <- b_agg %>% 
  mutate(
    trttype = case_match(
      trttype,
      "classic" ~ "Older ADs",
      "placebo" ~ "Placebo",
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2",
      .default = trttype
    ),,
    trttype = factor(
      trttype,
      levels = c("DPP4", "GLP-1", "SGLT2", "Older ADs", "Placebo")
    ),
    exp = case_match(
      exp,
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2"
    )
  ) %>% 
  ggplot(aes(x = attr_rate, y = trial_id, colour = trttype)) +
  geom_line(aes(group = trial_id), linewidth = 0.4, colour = "black") +
  geom_point() +
  scale_color_manual(
    values = c(
      "Older ADs" = "forestgreen",
      "DPP4" = "steelblue",
      "GLP-1" = "gold1",
      "SGLT2" = "darkorange1",
      "Placebo" = "orangered1"
    )
  ) +
  facet_wrap(~exp, scales = "free") +
  labs(x = "Attrition proportion", y = NULL, colour = "Treatment type") +
  theme_classic()

plot_rates

ggsave(
  "output/sum/plot_rates.png",
  plot_rates,
  width = 12,
  height = 8,
  units = "in"
)
