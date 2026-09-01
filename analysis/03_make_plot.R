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

k_sd <- sd(df$K)

# Population-average predictions (re.form = NA) at -1 SD/+1 SD of K, matching
# the course's plot_predictions() workflow for GLMMs.
p <- plot_predictions(
  fit2,
  condition = list(
    reward_oneback = c(0, 1),
    K_c = c(-k_sd, k_sd)
  ),
  re.form = NA,
  draw = FALSE
) |>
  mutate(K_level = factor(K_c, levels = c(-k_sd, k_sd), labels = c("-1 SD K", "+1 SD K")))

plt <- ggplot(p, aes(x = factor(reward_oneback), y = estimate, group = K_level, colour = K_level)) +
  geom_line(linewidth = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.6) +
  scale_colour_manual(values = c("-1 SD K" = "#E69F00", "+1 SD K" = "#0072B2"), name = NULL) +
  labs(x = "Previous-trial reward", y = "P(stay)",
       title = "Reward effect by working memory capacity") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11))

ggsave(OUT_PNG, plt, width = 3.75, height = 3, dpi = 96, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
