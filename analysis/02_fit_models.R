rm(list = ls())

#### SETUP ####
library(dplyr)
library(lme4)
library(lmerTest)
library(performance)
library(parameters)
library(marginaleffects)

# r2_pseudo(): course helper for the proportional reduction in a random-effect
# variance component between two nested models (used for GLMM effect sizes,
# since performance::r2()/r2_nakagawa() are not used for GLMMs in this course).
source("https://github.com/mattansb/Hierarchical-Linear-Models-foR-Psychologists/raw/refs/heads/main/helpers.R")

DATA_CSV <- "data/trial_level.csv"
LOG_FILE <- "analysis/model_log.txt"
FITS_DIR <- "analysis/fits"
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)

# Keep the sink guarded with a finally block, so an error partway through
# still restores console output instead of leaving it silently redirected.
tryCatch({

  df <- read.csv(DATA_CSV) |>
    mutate(participant = factor(participant))

  cat(sprintf("Loaded %d trials from %d subjects.\n\n", nrow(df), n_distinct(df$participant)))

  #### CHECK HETEROGENEITY OF reward_oneback (is it worth splitting WP/BP?) ####
  # Same logic as the course's within-person-fluctuation check: treat the
  # Level 1 predictor as an outcome and fit an empty random-intercept model to
  # get its between-person ICC. A near-zero ICC means there is no meaningful
  # between-person variance to separate out, so reward_oneback can enter the
  # models uncentred/unsplit.
  cat("#### Checking reward_oneback's own ICC (justifies leaving it unsplit) ####\n")
  fit_reward_icc <- glmer(
    reward_oneback ~ 1 + (1 | participant),
    family = binomial(link = "logit"),
    data = df
  )
  print(icc(fit_reward_icc))

  #### MODEL 0 — empty (unconditional), for the ICC ####
  cat("#### Fitting Model 0: stay ~ 1 + (1 | participant) ####\n")
  fit0 <- glmer(
    stay ~ 1 + (1 | participant),
    family = binomial(link = "logit"),
    data = df
  )
  print(model_parameters(fit0))

  cat("\nUnconditional ICC (performance::icc(), Model 0):\n")
  print(icc(fit0))
  cat(sprintf("\nLevel-1 variance is fixed by the logit link at pi^2/3 = %.4f\n", pi^2 / 3))

  #### MODEL 1 — reward_oneback only (baseline for the comparison) ####
  # Maximal random-effects structure (Barr et al., 2013): random intercept and
  # slope for reward_oneback, allowed to correlate. If this is singular or fails
  # to converge, apply the remedies in the order taught (control tweaks, then
  # optimizer, then dropping only the random *correlation* via || - never
  # dropping a random slope/intercept outright).
  cat("\n#### Fitting Model 1 (maximal RE): stay ~ reward_oneback + (reward_oneback | participant) ####\n")
  fit1 <- glmer(
    stay ~ reward_oneback + (reward_oneback | participant),
    family = binomial(link = "logit"),
    data = df
  )
  cat(sprintf("Model 1 singular: %s | converged: %s\n", check_singularity(fit1), is.null(fit1@optinfo$conv$lme4$messages)))

  if (check_singularity(fit1) || !is.null(fit1@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 1 with glmerControl(calc.derivs = FALSE)...\n")
    fit1 <- glmer(
      stay ~ reward_oneback + (reward_oneback | participant),
      family = binomial(link = "logit"),
      data = df,
      control = glmerControl(calc.derivs = FALSE)
    )
  }
  if (check_singularity(fit1) || !is.null(fit1@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 1 with the bobyqa optimizer...\n")
    fit1 <- glmer(
      stay ~ reward_oneback + (reward_oneback | participant),
      family = binomial(link = "logit"),
      data = df,
      control = glmerControl("bobyqa")
    )
  }
  if (check_singularity(fit1) || !is.null(fit1@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 1 dropping the random intercept-slope correlation (||)...\n")
    fit1 <- glmer(
      stay ~ reward_oneback + (reward_oneback || participant),
      family = binomial(link = "logit"),
      data = df
    )
  }
  print(model_parameters(fit1, exponentiate = TRUE))

  #### MODEL 2 — reward_oneback x K_c cross-level interaction (key model) ####
  cat("\n#### Fitting Model 2 (maximal RE): stay ~ reward_oneback * K_c + (reward_oneback | participant) ####\n")
  fit2 <- glmer(
    stay ~ reward_oneback * K_c + (reward_oneback | participant),
    family = binomial(link = "logit"),
    data = df
  )
  cat(sprintf("Model 2 singular: %s | converged: %s\n", check_singularity(fit2), is.null(fit2@optinfo$conv$lme4$messages)))

  if (check_singularity(fit2) || !is.null(fit2@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 2 with glmerControl(calc.derivs = FALSE)...\n")
    fit2 <- glmer(
      stay ~ reward_oneback * K_c + (reward_oneback | participant),
      family = binomial(link = "logit"),
      data = df,
      control = glmerControl(calc.derivs = FALSE)
    )
  }
  if (check_singularity(fit2) || !is.null(fit2@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 2 with the bobyqa optimizer...\n")
    fit2 <- glmer(
      stay ~ reward_oneback * K_c + (reward_oneback | participant),
      family = binomial(link = "logit"),
      data = df,
      control = glmerControl("bobyqa")
    )
  }
  if (check_singularity(fit2) || !is.null(fit2@optinfo$conv$lme4$messages)) {
    cat("Retrying Model 2 dropping the random intercept-slope correlation (||)...\n")
    fit2 <- glmer(
      stay ~ reward_oneback * K_c + (reward_oneback || participant),
      family = binomial(link = "logit"),
      data = df
    )
  }
  print(model_parameters(fit2, exponentiate = TRUE))

  #### MODEL COMPARISON: Model 1 vs Model 2 ####
  cat("\n#### Model comparison (likelihood-ratio test; GLMMs are ML only, no refit needed) ####\n")
  comparison <- anova(fit1, fit2)
  print(comparison)

  cat("\nPseudo R2 (reduction in random-effect variance, Model 2 vs Model 1):\n")
  print(r2_pseudo(fit2, fit1))

  #### FIXED EFFECTS: odds ratios and average marginal effects ####
  cat("\n#### Model 2 fixed effects: odds ratios (model_parameters) ####\n")
  print(model_parameters(fit2, exponentiate = TRUE))

  cat("\n#### Model 2: average marginal effect of reward_oneback on P(stay) ####\n")
  print(avg_slopes(fit2, variables = "reward_oneback", re.form = NULL))

  cat("\n#### Model 2: average odds ratio for reward_oneback ####\n")
  print(avg_comparisons(fit2, variables = "reward_oneback", comparison = "lnor", transform = "exp"))

  #### ROBUSTNESS: block-within-participant grouping (Model 2, three levels) ####
  # Trials come in six blocks per person; add a grouping factor for block
  # within participant (multiple random grouping variables, week 7 syntax:
  # (1 | childid:schoolid)) to see whether the interaction survives once
  # within-session drift gets its own variance term.
  cat("\n#### Model 2 with an added (1 | participant:block_number) term ####\n")
  fit2_block <- glmer(
    stay ~ reward_oneback * K_c + (reward_oneback | participant) + (1 | participant:block_number),
    family = binomial(link = "logit"),
    data = df
  )
  print(model_parameters(fit2_block, exponentiate = TRUE))

  #### CONVERGENCE / SINGULARITY SUMMARY ####
  cat("\n#### Convergence diagnostics ####\n")
  cat(sprintf("Model 1 singular: %s\n", check_singularity(fit1)))
  cat(sprintf("Model 2 singular: %s\n", check_singularity(fit2)))
  print(check_convergence(fit1))
  print(check_convergence(fit2))

  results <- list(
    fit_reward_icc = fit_reward_icc,
    fit0 = fit0,
    fit1 = fit1,
    fit2 = fit2,
    fit2_block = fit2_block,
    icc0 = icc(fit0),
    comparison = comparison,
    r2_pseudo = r2_pseudo(fit2, fit1)
  )
  saveRDS(results, file.path(FITS_DIR, "results_summary.rds"))

  cat(sprintf("\nDone. Log written to %s\n", LOG_FILE))

}, finally = {
  if (sink.number() > 0) sink()
})
