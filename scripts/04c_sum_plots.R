
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

# Join to rest, add comparator id
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
    ),
    attr_prop = round(attr_rate * 100, 2),
    comp = case_when(
      exp != trttype ~ trttype,
      TRUE ~ NA
    )
  ) %>% 
  group_by(trial_id) %>% 
  fill(comp, .direction = "downup") %>% 
  mutate(
    comp = case_when(
      length(unique(comp)) == 2 & "placebo" %in% comp ~ "placebo",
      TRUE ~ comp
    )
  ) %>% 
  ungroup()

# Arrange by total attrition per trial, compute total attrition
# Setup factors for plot
b_agg_index <- b_agg %>% 
  group_by(trial_id) %>% 
  mutate(attr_total = round(sum(arm_attr)/sum(arm_n) * 100, 2)) %>% 
  ungroup() %>% 
  arrange(attr_total) %>% 
  mutate(
    order_id = match(trial_id, unique(trial_id)),
    trttype = case_match(
      trttype,
      "classic" ~ "Other ADs",
      "placebo" ~ "Placebo",
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2",
      .default = trttype
    ),
    trttype = factor(
      trttype,
      levels = c("DPP4", "GLP-1", "SGLT2", "Other ADs", "Placebo")
    ),
    exp = case_match(
      exp,
      "dpp4" ~ "DPP4",
      "glp1" ~ "GLP-1",
      "sglt2" ~ "SGLT2"
    ),
    comp = fct_relevel(factor(comp), "placebo", "classic", "dpp4")
  ) %>% 
  arrange(comp, order_id, trial_id) %>% 
  mutate(trial_id = factor(trial_id, levels = unique(trial_id)))

# Plot
plot_rates <- b_agg_index %>% 
  ggplot(aes(x = attr_prop, y = fct_rev(trial_id), colour = trttype)) +
  geom_line(aes(group = trial_id), linewidth = 0.4, colour = "black") +
  geom_point() +
  scale_color_manual(
    values = c(
      "Other ADs" = "forestgreen",
      "DPP4" = "steelblue",
      "GLP-1" = "gold1",
      "SGLT2" = "darkorange1",
      "Placebo" = "orangered1"
    )
  ) +
  facet_wrap(~exp, scales = "free") +
  labs(x = "Total attrition (%)", y = NULL, colour = "Treatment type") +
  scale_x_continuous(limits = c(0, 60), n.breaks = 6) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 15, colour = "white")
  )

plot_rates

ggsave(
  "output/sum/plot_rates.png",
  plot_rates,
  width = 14,
  height = 8,
  units = "in"
)

