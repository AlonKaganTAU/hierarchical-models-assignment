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

# K_c is constant within a participant (repeated once per trial), so its SD
# must be computed over the 35 subjects, not over all 5,025 trial rows.
k_participant <- unique(df[c("participant", "K_c")])
k_sd <- sd(k_participant$K_c)

#### PANEL A: model-implied predictions ####
# Let plot_predictions() draw the figure directly and add layers on top, as
# in every course example (week 3's own +-1 SD moderator plot, week 6's GLMM
# plot_predictions() usage), instead of extracting a data frame by hand.
plt_A <- plot_predictions(
  fit2,
  condition = list("reward_oneback", K_c = c(-k_sd, k_sd)),
  re.form = NA
) +
  scale_color_manual(
    "K",
    values = c("#E69F00", "#0072B2"),
    labels = c("-1 SD K", "+1 SD K"),
    aesthetics = c("fill", "color")
  ) +
  labs(x = "Previous-trial reward", y = "P(stay)", title = "Model-implied WSLS effect") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11))

#### PANEL B: each participant's own observed WSLS effect against their K ####
# Mirrors GLMM.R's own precedent for visualizing a Level-2 moderator against a
# subject-level outcome: ggplot(SFON_data, aes(weberFr, Attend)) +
# stat_summary(aes(group = ID), geom = "point"). Our outcome is a within-
# subject contrast (P(stay) after reward minus after no reward) rather than a
# single proportion, so it is computed with group_by()/summarise() first and
# then plotted the same way, with a plain OLS reference line (geom_smooth(
# method = "lm", se = FALSE)) matching the course's BLUP-scatter convention
# (Ex3solution.R) rather than a model-based band.
participant_wsls_effects <- df |>
  group_by(participant, K) |>
  summarise(
    wsls_effect = mean(stay[reward_oneback == 1]) - mean(stay[reward_oneback == 0]),
    .groups = "drop"
  )

plt_B <- ggplot(participant_wsls_effects, aes(K, wsls_effect)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Working memory capacity (K)", y = "Observed WSLS effect",
       title = "Observed, by participant") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

#### COMBINE ####
# plot_layout(guides = "collect") merges the shared legend, matching the
# course's own patchwork usage (week 7).
plt <- plt_A + plt_B + plot_layout(guides = "collect")

ggsave(OUT_PNG, plt, width = 8, height = 3.5, dpi = 300, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
