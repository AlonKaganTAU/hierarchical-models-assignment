rm(list = ls())

#### SETUP ####
library(dplyr)
library(marginaleffects)
library(ggplot2)

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

# Let plot_predictions() draw the figure directly and add layers on top, as
# in every course example (week 3's own +-1 SD moderator plot, week 6's GLMM
# plot_predictions() usage), instead of extracting a data frame by hand.
plt <- plot_predictions(
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
  labs(x = "Previous-trial reward", y = "P(stay)",
       title = "Reward effect by working memory capacity") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11))

ggsave(OUT_PNG, plt, width = 3.75, height = 3, dpi = 96, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
