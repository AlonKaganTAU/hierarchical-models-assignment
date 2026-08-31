rm(list = ls())

#### SETUP ####
library(dplyr)
library(readr)
library(brms)
library(tidybayes)
library(ggplot2)
library(patchwork)

FITS_DIR <- "analysis/fits"
OUT_PNG  <- "analysis/figures/simple_slopes.png"
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

fit2 <- readRDS(file.path(FITS_DIR, "fit2_final.rds"))
df   <- read_csv("data/trial_level.csv", show_col_types = FALSE)

# The SD must be taken over the 35 PARTICIPANTS. Taking sd() over the 5,025
# rows would weight each participant by their number of trials and give the
# wrong "+/- 1 SD" values.
K_by_subject <- df |> distinct(participant, K)
k_sd <- sd(K_by_subject$K)

#### PANEL A - reward effect at low / average / high WM capacity ####

# Population-level predictions (re_formula = NA), i.e. for a participant whose
# random effects are 0. Because the logit link is non-linear these are the
# MEDIAN participant's probabilities, not the average across participants; the
# population-average effect is reported separately in the text (see 02).
grid_a <- expand.grid(reward = c(0, 1), K_c = c(-k_sd, 0, k_sd)) |>
  mutate(K_level = factor(K_c,
                          levels = c(-k_sd, 0, k_sd),
                          labels = c("Low K (−1 SD)", "Mean K", "High K (+1 SD)")))

preds_a <- fitted(fit2, newdata = grid_a, re_formula = NA, probs = c(0.025, 0.975)) |>
  as_tibble() |>
  bind_cols(grid_a)

p_a <- ggplot(preds_a, aes(factor(reward), Estimate, colour = K_level, group = K_level)) +
  geom_line(linewidth = 0.8) +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5), size = 0.45,
                  position = position_dodge(width = 0.08)) +
  scale_x_discrete(labels = c("No reward", "Reward")) +
  scale_y_continuous(limits = c(0, 0.8), labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("#E69F00", "#666666", "#0072B2"), name = NULL) +
  labs(x = "Previous trial", y = "P(stay)", tag = "A") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.key.height = unit(8, "pt"))

#### PANEL B - the moderation itself ####

# Each participant's observed reward effect against their WM capacity, with the
# model-implied effect overlaid. This is the cross-level interaction plotted
# directly: the slope of the line is gamma_11 (on the probability scale).
observed <- df |>
  group_by(participant, K, reward) |>
  summarise(p_stay = mean(stay), .groups = "drop") |>
  tidyr::pivot_wider(names_from = reward, values_from = p_stay, names_prefix = "r") |>
  mutate(reward_effect = r1 - r0)

grid_b <- tibble(K = seq(min(df$K), max(df$K), length.out = 60)) |>
  mutate(K_c = K - mean(K_by_subject$K))

preds_b <- grid_b |>
  tidyr::expand_grid(reward = c(0, 1)) |>
  add_epred_draws(fit2, re_formula = NA) |>
  ungroup() |>
  select(K, reward, .draw, .epred) |>
  tidyr::pivot_wider(names_from = reward, values_from = .epred, names_prefix = "r") |>
  mutate(reward_effect = r1 - r0) |>
  group_by(K) |>
  summarise(Estimate = median(reward_effect),
            Q2.5     = quantile(reward_effect, 0.025),
            Q97.5    = quantile(reward_effect, 0.975))

p_b <- ggplot(preds_b, aes(K, Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "#0072B2", alpha = 0.15) +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  geom_point(data = observed, aes(K, reward_effect), alpha = 0.6, size = 1.4,
             inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Working memory capacity (K)",
       y = "Reward effect on P(stay)", tag = "B") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

ggsave(OUT_PNG, p_a + p_b, width = 6.3, height = 2.9, dpi = 300, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
