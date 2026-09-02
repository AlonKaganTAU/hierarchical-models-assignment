rm(list = ls())

#### SETUP ####
library(dplyr)
library(marginaleffects)
library(ggplot2)
library(patchwork)

FITS_DIR <- "analysis/fits"
OUT_PNG  <- "analysis/figures/simple_slopes.png"
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

results <- readRDS(file.path(FITS_DIR, "results_summary.rds"))
fit2 <- results$fit2
df   <- read.csv("data/trial_level.csv")

# K_c is constant within a subject (repeated once per trial), so its SD
# must be computed over the 35 subjects, not over all 5,025 trial rows.
k_subject <- unique(df[c("subject", "K_c")])
k_sd <- sd(k_subject$K_c)

#### PANEL A: model-implied predictions ####
# Build the panel from the prediction data frame, so each K level gets its
# own line and colour.
preds_A <- plot_predictions(
  fit2,
  condition = list("reward_oneback", K_c = c(-k_sd, 0, k_sd)),
  re.form = NA,
  draw = FALSE
) |>
  mutate(K_level = factor(K_c, levels = c(-k_sd, 0, k_sd),
                           labels = c("Low K (-1 SD)", "Mean K", "High K (+1 SD)")),
         reward_lab = factor(reward_oneback, levels = c(0, 1),
                              labels = c("No reward", "Reward")))

plt_A <- ggplot(preds_A, aes(reward_lab, estimate, color = K_level, group = K_level)) +
  geom_line(position = position_dodge(width = 0.15)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), shape = 18, size = 0.8,
                   position = position_dodge(width = 0.15)) +
  scale_color_manual(values = c("#E69F00", "#666666", "#0072B2")) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(x = "Previous trial", y = "P(stay)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_blank())

#### PANEL B: each subject's own observed WSLS effect against their K ####
# Each subject's own WSLS effect (a within-subject contrast: P(stay)
# after reward minus after no reward) is computed with group_by()/
# summarise() first, then plotted with a linear reference line and its
# confidence band.
subject_wsls_effects <- df |>
  group_by(subject, K) |>
  summarise(
    wsls_effect = mean(stay[reward_oneback == 1]) - mean(stay[reward_oneback == 0]),
    .groups = "drop"
  )

plt_B <- ggplot(subject_wsls_effects, aes(K, wsls_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", color = "#0072B2", fill = "#0072B2") +
  geom_point(color = "grey30") +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Working memory capacity (K)", y = "WSLS effect on P(stay)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

#### COMBINE ####
plt <- (plt_A + plt_B) + plot_annotation(tag_levels = "A")

ggsave(OUT_PNG, plt, width = 9, height = 4, dpi = 300, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
