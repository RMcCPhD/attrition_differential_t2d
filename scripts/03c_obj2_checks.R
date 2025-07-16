
source("scripts/00_config.R")
source("scripts/00_packages.R")

a_imp_mdl <- readRDS("processed_data/inter_fit.rds")

# Model summary
b_sum <- summary(a_imp_mdl)

# Extract interactions summary
b_sum_inter <- b_sum %>% 
  as_tibble() %>% 
  filter(grepl("beta", parameter)) %>% 
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("beta", "", parameter)
  ) %>% 
  separate(
    parameter,
    into = c("term", "class"),
    sep = "\\."
  ) %>% 
  mutate(
    class = gsub("trtclass", "", class),
    term = gsub("\\:", "", term),
    mean_or = exp(mean)
  )

# Treatment-covariate interation effect estimate -------------------------------

# Plot log-odds
plot_log <- b_sum_inter %>% 
  # filter(!is.na(class)) %>% 
  # ggplot(aes(x = mean, y = term)) +
  # geom_point(position = position_dodge(width = 0.5)) +
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = term)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  facet_wrap(~class, scales = "free") +
  theme_bw() +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Mean log-odds (95% credible intervals)", y = NULL)

plot_log

ggsave(
  "output/obj1/plot_log_odds.png",
  plot_log,
  width = 8,
  height = 4,
  units = "in"
)



