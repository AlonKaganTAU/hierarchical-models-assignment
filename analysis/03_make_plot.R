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
# Get predictions() as a data frame (as in week 4) and build the plot by hand,
# so the three K levels are guaranteed their own line and colour.
preds_A <- plot_predictions(
  fit2,
  condition = list("reward_oneback", K_c = c(-k_sd, 0, k_sd)),
  re.form = NA,
  draw = FALSE
) |>
  mutate(K_level = factor(K_c, levels = c(-k_sd, 0, k_sd),
                           labels = c("-1 SD K", "Mean K", "+1 SD K")))

plt_A <- ggplot(preds_A, aes(factor(reward_oneback), estimate, colour = K_level, group = K_level)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_colour_manual("K", values = c("#E69F00", "#666666", "#0072B2")) +
  labs(x = "Previous-trial reward", y = "P(stay)", title = "Model-implied WSLS effect") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11))

#### PANEL B: each subject's own observed WSLS effect against their K ####
# Each subject's own WSLS effect (a within-subject contrast: P(stay)
# after reward minus after no reward) is computed with group_by()/
# summarise() first, then plotted with a plain linear reference line rather
# than a model-based band.
subject_wsls_effects <- df |>
  group_by(subject, K) |>
  summarise(
    wsls_effect = mean(stay[reward_oneback == 1]) - mean(stay[reward_oneback == 0]),
    .groups = "drop"
  )

plt_B <- ggplot(subject_wsls_effects, aes(K, wsls_effect)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Working memory capacity (K)", y = "Observed WSLS effect",
       title = "Observed, by subject") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

#### COMBINE ####
# plot_layout(guides = "collect") merges the shared legend; the theme()
# after "&" applies to the combined figure so the merged legend keeps
# Panel A's bottom placement instead of defaulting to the right.
plt <- (plt_A + plt_B + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")

ggsave(OUT_PNG, plt, width = 8, height = 3.5, dpi = 300, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
