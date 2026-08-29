rm(list = ls())

#### SETUP ####
library(dplyr)
library(brms)
library(loo)
library(posterior)
library(bayesplot)

DATA_CSV  <- "data/trial_level.csv"
FITS_DIR  <- "analysis/fits"
LOG_FILE  <- "analysis/model_log.txt"
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

CHAINS  <- 4
ITER    <- 3000
WARMUP  <- 1000
CORES   <- 4

log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)

df <- read.csv(DATA_CSV) |>
  mutate(participant = factor(participant))

cat(sprintf("Loaded %d trials from %d subjects.\n\n", nrow(df), n_distinct(df$participant)))

# A model has converged if there are no post-warmup divergent transitions and
# every Rhat is below 1.01; anything else falls back to the simpler
# random-effects structure, and both attempts are kept for reporting.
model_has_converged <- function(fit) {
  divergences <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  rhat_ok     <- all(rhat(fit) < 1.01, na.rm = TRUE)
  divergences == 0 && rhat_ok
}

#### MODEL 0 — unconditional (intercept-only), for the ICC ####
cat("#### Fitting Model 0: stay ~ 1 + (1 | participant) ####\n")
fit0 <- brm(stay ~ 1 + (1 | participant), data = df, family = bernoulli(),
            chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
            refresh = 0, backend = "rstan")
print(summary(fit0))
saveRDS(fit0, file.path(FITS_DIR, "fit0.rds"))

tau00_draws <- as_draws_df(fit0)$sd_participant__Intercept^2
icc_draws   <- tau00_draws / (tau00_draws + pi^2 / 3)
icc_summary <- c(median = median(icc_draws), quantile(icc_draws, c(0.025, 0.975)))
cat("\nUnconditional ICC (latent-threshold method):\n")
print(icc_summary)

#### MODEL 1 — reward only (baseline for the comparison) ####
cat("\n#### Fitting Model 1 (maximal RE): stay ~ reward + (reward | participant) ####\n")
fit1_max <- brm(stay ~ reward + (reward | participant), data = df, family = bernoulli(),
                 chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
                 refresh = 0, backend = "rstan")
print(summary(fit1_max))
saveRDS(fit1_max, file.path(FITS_DIR, "fit1_maximal.rds"))

fit1_converged <- model_has_converged(fit1_max)
cat(sprintf("\nModel 1 maximal RE converged: %s\n", fit1_converged))

if (fit1_converged) {
  fit1 <- fit1_max
} else {
  cat("\n#### Falling back: Model 1 with (1 | participant) only ####\n")
  fit1_fallback <- brm(stay ~ reward + (1 | participant), data = df, family = bernoulli(),
                        chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
                        refresh = 0, backend = "rstan")
  print(summary(fit1_fallback))
  saveRDS(fit1_fallback, file.path(FITS_DIR, "fit1_fallback.rds"))
  fit1 <- fit1_fallback
}
saveRDS(fit1, file.path(FITS_DIR, "fit1_final.rds"))

#### MODEL 2 — reward x K_c cross-level interaction (key model) ####
cat("\n#### Fitting Model 2 (maximal RE): stay ~ reward * K_c + (reward | participant) ####\n")
fit2_max <- brm(stay ~ reward * K_c + (reward | participant), data = df, family = bernoulli(),
                 chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
                 refresh = 0, backend = "rstan")
print(summary(fit2_max))
saveRDS(fit2_max, file.path(FITS_DIR, "fit2_maximal.rds"))

fit2_converged <- model_has_converged(fit2_max)
cat(sprintf("\nModel 2 maximal RE converged: %s\n", fit2_converged))

if (fit2_converged) {
  fit2 <- fit2_max
} else {
  cat("\n#### Falling back: Model 2 with (1 | participant) only ####\n")
  fit2_fallback <- brm(stay ~ reward * K_c + (1 | participant), data = df, family = bernoulli(),
                        chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
                        refresh = 0, backend = "rstan")
  print(summary(fit2_fallback))
  saveRDS(fit2_fallback, file.path(FITS_DIR, "fit2_fallback.rds"))
  fit2 <- fit2_fallback
}
saveRDS(fit2, file.path(FITS_DIR, "fit2_final.rds"))

#### MODEL COMPARISON: Model 1 vs Model 2 ####
cat("\n#### Model comparison ####\n")
fit1 <- add_criterion(fit1, "loo")
fit2 <- add_criterion(fit2, "loo")
comparison <- loo_compare(fit1, fit2)
print(comparison)

r2_fit1_marginal    <- bayes_R2(fit1, re_formula = NA)
r2_fit1_conditional <- bayes_R2(fit1)
r2_fit2_marginal    <- bayes_R2(fit2, re_formula = NA)
r2_fit2_conditional <- bayes_R2(fit2)

cat("\nModel 1 bayes_R2 (marginal):\n");    print(r2_fit1_marginal)
cat("Model 1 bayes_R2 (conditional):\n");   print(r2_fit1_conditional)
cat("Model 2 bayes_R2 (marginal):\n");      print(r2_fit2_marginal)
cat("Model 2 bayes_R2 (conditional):\n");   print(r2_fit2_conditional)

#### INTERACTION EFFECT (primary result) ####
interaction_draws <- as_draws_df(fit2)$`b_reward:K_c`
interaction_summary <- c(median = median(interaction_draws), quantile(interaction_draws, c(0.025, 0.975)))
interaction_pd <- max(mean(interaction_draws > 0), mean(interaction_draws < 0))

cat("\nInteraction (reward x K_c) posterior summary:\n")
print(interaction_summary)
cat(sprintf("Probability of direction (pd): %.4f\n", interaction_pd))

saveRDS(fit0, file.path(FITS_DIR, "fit0.rds"))

results <- list(
  icc_draws            = icc_summary,
  fit1_converged       = fit1_converged,
  fit2_converged       = fit2_converged,
  loo_comparison       = comparison,
  r2_fit1_marginal     = r2_fit1_marginal,
  r2_fit1_conditional  = r2_fit1_conditional,
  r2_fit2_marginal     = r2_fit2_marginal,
  r2_fit2_conditional  = r2_fit2_conditional,
  interaction_summary  = interaction_summary,
  interaction_pd       = interaction_pd
)
saveRDS(results, file.path(FITS_DIR, "results_summary.rds"))

sink()
cat(sprintf("\nDone. Log written to %s, fitted models saved to %s/\n", LOG_FILE, FITS_DIR))
