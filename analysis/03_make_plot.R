rm(list = ls())

#### SETUP ####
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(patchwork)

FITS_DIR <- "analysis/fits"
OUT_PNG  <- "analysis/figures/reward_effect_by_K.png"
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)

fit3 <- readRDS(file.path(FITS_DIR, "fit3_final.rds"))
df   <- read.csv("data/trial_level.csv")

subjects <- df |>
  summarise(K = first(K), K_c = first(K_c), .by = participant)

REWARD_COLOURS <- c("0" = "#B8752B", "1" = "#1F6FB4")

#### POPULATION-LEVEL PREDICTIONS ACROSS THE OBSERVED RANGE OF K ####
grid <- expand_grid(reward_oneback = c(0, 1),
                    K_c = seq(min(subjects$K_c), max(subjects$K_c), length.out = 60))

predicted <- fitted(fit3, newdata = grid, re_formula = NA, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  bind_cols(grid) |>
  mutate(K = K_c + mean(subjects$K - subjects$K_c),
         reward_oneback = factor(reward_oneback))

# Observed stay rates per participant, for context behind the model lines.
observed <- df |>
  summarise(stay = mean(stay), .by = c(participant, K, reward_oneback)) |>
  mutate(reward_oneback = factor(reward_oneback))

panel_predictions <- ggplot(predicted, aes(K, Estimate, colour = reward_oneback, fill = reward_oneback)) +
  geom_point(aes(y = stay), data = observed, alpha = 0.45, size = 1.1) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = REWARD_COLOURS, name = "Previous trial",
                      labels = c("no reward", "reward")) +
  scale_fill_manual(values = REWARD_COLOURS, guide = "none") +
  labs(x = "Working memory capacity (Cowan's K)", y = "P(stay)", tag = "A") +
  coord_cartesian(ylim = c(0, 1))

#### THE REWARD EFFECT ITSELF, AS A FUNCTION OF K ####
# Difference in predicted P(stay) between rewarded and unrewarded previous
# trials, computed draw by draw so that the interval is the interval of the
# effect and not of the two predictions separately.
draws_reward <- posterior_epred(fit3, newdata = filter(grid, reward_oneback == 1), re_formula = NA) -
                posterior_epred(fit3, newdata = filter(grid, reward_oneback == 0), re_formula = NA)

reward_effect <- filter(grid, reward_oneback == 1) |>
  mutate(K        = K_c + mean(subjects$K - subjects$K_c),
         estimate = apply(draws_reward, 2, median),
         lower    = apply(draws_reward, 2, quantile, 0.025),
         upper    = apply(draws_reward, 2, quantile, 0.975))

observed_effect <- observed |>
  pivot_wider(names_from = reward_oneback, values_from = stay, names_prefix = "stay_") |>
  mutate(effect = stay_1 - stay_0)

panel_effect <- ggplot(reward_effect, aes(K, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(y = effect), data = observed_effect, alpha = 0.45, size = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, fill = "#1F6FB4") +
  geom_line(linewidth = 0.8, colour = "#1F6FB4") +
  labs(x = "Working memory capacity (Cowan's K)",
       y = "Reward effect on P(stay)", tag = "B")

#### COMBINE AND SAVE ####
figure <- (panel_predictions | panel_effect) &
  theme_bw(base_size = 8) &
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(t = -4),
        plot.tag = element_text(face = "bold", size = 9))

ggsave(OUT_PNG, figure, width = 5.4, height = 1.95, dpi = 300, bg = "white")
cat(sprintf("Saved -> %s\n", OUT_PNG))
