rm(list = ls())

#### SETUP ####
library(dplyr)
library(brms)
library(ggplot2)

FITS_DIR <- "analysis/fits"
OUT_PNG  <- "analysis/figures/simple_slopes.png"
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

fit2 <- readRDS(file.path(FITS_DIR, "fit2_final.rds"))
df   <- read.csv("data/trial_level.csv")

k_sd <- sd(df$K)

newdata <- expand.grid(
  reward = c(0, 1),
  K_c    = c(-k_sd, k_sd)
) |>
  mutate(K_level = factor(K_c, levels = c(-k_sd, k_sd), labels = c("-1 SD K", "+1 SD K")))

preds <- fitted(fit2, newdata = newdata, re_formula = NA, summary = TRUE, probs = c(0.025, 0.975))
plot_df <- bind_cols(newdata, as.data.frame(preds))

p <- ggplot(plot_df, aes(x = factor(reward), y = Estimate, group = K_level, colour = K_level)) +
  geom_line(linewidth = 1) +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), size = 0.6) +
  scale_colour_manual(values = c("-1 SD K" = "#E69F00", "+1 SD K" = "#0072B2"), name = NULL) +
  labs(x = "Previous-trial reward", y = "P(stay)",
       title = "Reward effect by working memory capacity") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11))

ggsave(OUT_PNG, p, width = 3.75, height = 3, dpi = 96, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
