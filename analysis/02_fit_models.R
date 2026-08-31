rm(list = ls())

#### SETUP ####
library(dplyr)
library(readr)
library(brms)
library(loo)
library(posterior)
library(bayesplot)
library(bayestestR)
library(marginaleffects)
library(performance)
library(parameters)

DATA_CSV <- "data/trial_level.csv"
FITS_DIR <- "analysis/fits"
LOG_FILE <- "analysis/model_log.txt"
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

CHAINS <- 4
ITER   <- 4000
WARMUP <- 1500
CORES  <- 4
SEED   <- 20260829 # fixed so the reported numbers are reproducible

log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)
on.exit({sink(); close(log_con)}, add = TRUE) # so an error can't leave the console muted

df <- read_csv(DATA_CSV, show_col_types = FALSE) |>
  mutate(participant = factor(participant),
         block_id    = factor(paste(participant, block_number, sep = "_")))

cat(sprintf("Loaded %d trials from %d subjects.\n", nrow(df), n_distinct(df$participant)))
cat(sprintf("SD of K between subjects: %.3f\n\n", sd(distinct(df, participant, K)$K)))

# Weakly-informative priors. brms' default for population-level slopes is an
# IMPROPER FLAT prior, which we do not want: it leaves the interaction
# unregularised and makes Bayes factors undefined. Normal(0, 1.5) on the logit
# scale is wide (it puts ~95% of the prior mass on odds ratios between about
# 1/19 and 19) but rules out the absurd.
priors <- c(
  prior(normal(0, 1.5),      class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(lkj(2),               class = "cor")
)

# Fits (or, on a rerun, loads) one model. The cache is keyed on the formula as
# well as the name, so editing a model can never silently load a stale fit.
fit_or_load <- function(formula, name, prior = priors) {
  path <- file.path(FITS_DIR, paste0(name, ".rds"))
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (identical(deparse(formula), deparse(cached$formula$formula))) {
      cat(sprintf("Loading cached fit: %s\n", path))
      return(cached)
    }
    cat(sprintf("Cached fit %s does not match the current formula - refitting.\n", path))
  }
  fit <- brm(formula, data = df, family = bernoulli(), prior = prior,
             chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
             seed = SEED, refresh = 0, backend = "rstan")
  saveRDS(fit, path)
  fit
}

# A model has converged if there are no post-warmup divergent transitions, every
# Rhat is below 1.01, and every effective sample size is at least 400 (Vehtari
# et al., 2021). Anything else falls back to the simpler random-effects
# structure, and both attempts are kept for reporting.
model_has_converged <- function(fit) {
  divergences <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  draws_smry  <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  cat(sprintf("  divergences = %d, max Rhat = %.4f, min ESS = %.0f\n",
              divergences, max(draws_smry$rhat, na.rm = TRUE),
              min(c(draws_smry$ess_bulk, draws_smry$ess_tail), na.rm = TRUE)))
  divergences == 0 &&
    max(draws_smry$rhat, na.rm = TRUE) < 1.01 &&
    min(c(draws_smry$ess_bulk, draws_smry$ess_tail), na.rm = TRUE) >= 400
}

post_summary <- function(x, label) {
  cat(sprintf("%-34s %6.3f [%6.3f, %6.3f]  pd = %.3f\n", label,
              median(x), quantile(x, 0.025), quantile(x, 0.975),
              max(mean(x > 0), mean(x < 0))))
}

#### MODEL 0 - unconditional (intercept-only), for the ICC ####

# Level 1:  logit(P(stay_ij)) = b_0j
# Level 2:              b_0j = gamma_00 + U_0j
# Composite: logit(P(stay_ij)) = gamma_00 + U_0j
#
# In a logistic HLM there is no estimated level 1 variance: the residual is
# fixed by the logit link at pi^2 / 3 = 3.29 on the latent scale, so the ICC is
# the latent-threshold ICC.

cat("\n#### Model 0: stay ~ 1 + (1 | participant) ####\n")
fit0 <- fit_or_load(stay ~ 1 + (1 | participant), "fit0",
                    prior = c(prior(student_t(3, 0, 2.5), class = "Intercept"),
                              prior(student_t(3, 0, 2.5), class = "sd")))
print(summary(fit0))
cat("Model 0 converged:", model_has_converged(fit0), "\n")

tau00_draws <- as_draws_df(fit0)$sd_participant__Intercept^2
icc_draws   <- tau00_draws / (tau00_draws + pi^2 / 3)
cat("\nUnconditional ICC (latent-threshold):\n")
post_summary(icc_draws, "ICC")

# Cross-check against performance::icc(), which uses the same decomposition:
print(performance::icc(fit0))

#### IS reward A "SMUSHED" LEVEL 1 PREDICTOR? ####

# A measured time-varying predictor usually carries between-person as well as
# within-person variance, and the two must be separated (person-mean centring)
# or their effects are conflated. We check how much of reward's variance is
# between-person by treating reward itself as the outcome of an empty model.
cat("\n#### Between-person variance in the level 1 predictor (reward) ####\n")
fit_rw <- fit_or_load(reward ~ 1 + (1 | participant), "fit_reward_icc",
                      prior = c(prior(student_t(3, 0, 2.5), class = "Intercept"),
                                prior(student_t(3, 0, 2.5), class = "sd")))
tau_rw    <- as_draws_df(fit_rw)$sd_participant__Intercept^2
post_summary(tau_rw / (tau_rw + pi^2 / 3), "ICC of reward")
# Reward was delivered by a fixed, choice-independent probability, so any
# between-person variation here is sampling noise. With a negligible ICC there
# is no between-person component to separate, and reward enters uncentred - its
# 0 point (an unrewarded previous trial) is already meaningful.

#### MODEL 1 - reward only ####

# Level 1:  logit(P(stay_ij)) = b_0j + b_1j * reward_ij
# Level 2:              b_0j = gamma_00 + U_0j
#                       b_1j = gamma_10 + U_1j
# Composite: logit(P(stay_ij)) = (gamma_00 + U_0j) + (gamma_10 + U_1j) * reward_ij
#
# U_0j and U_1j are allowed to covary (Hoffman's notation cannot show this, but
# the model estimates it). reward is a level 1 predictor, so per Barr et al.
# (2013) it gets a random slope: we want to generalise the reward effect beyond
# these 35 people.

cat("\n#### Model 1 (maximal RE): stay ~ reward + (reward | participant) ####\n")
fit1_max <- fit_or_load(stay ~ reward + (reward | participant), "fit1_maximal")
print(summary(fit1_max))

fit1_converged <- model_has_converged(fit1_max)
cat(sprintf("Model 1 maximal RE converged: %s\n", fit1_converged))

if (fit1_converged) {
  fit1 <- fit1_max
} else {
  cat("\n#### Falling back: Model 1 with (1 | participant) only ####\n")
  fit1 <- fit_or_load(stay ~ reward + (1 | participant), "fit1_fallback",
                      prior = priors[priors$class != "cor", ])
  print(summary(fit1))
}
saveRDS(fit1, file.path(FITS_DIR, "fit1_final.rds"))

#### MODEL 2 - reward x K_c cross-level interaction (key model) ####

# Level 1:  logit(P(stay_ij)) = b_0j + b_1j * reward_ij
# Level 2:              b_0j = gamma_00 + gamma_01 * K_c_j + U_0j
#                       b_1j = gamma_10 + gamma_11 * K_c_j + U_1j
# Composite: logit(P(stay_ij)) = (gamma_00 + gamma_01 * K_c_j + U_0j) +
#                                (gamma_10 + gamma_11 * K_c_j + U_1j) * reward_ij
#
# gamma_11 is the cross-level interaction: reward is measured at level 1, K at
# level 2, so K moderating the reward slope is a level 2 predictor of a level 1
# effect. K_c is grand-mean centred, so gamma_10 is the reward effect for a
# participant of average WM capacity.

cat("\n#### Model 2 (maximal RE): stay ~ reward * K_c + (reward | participant) ####\n")
fit2_max <- fit_or_load(stay ~ reward * K_c + (reward | participant), "fit2_maximal")
print(summary(fit2_max))

fit2_converged <- model_has_converged(fit2_max)
cat(sprintf("Model 2 maximal RE converged: %s\n", fit2_converged))

if (fit2_converged) {
  fit2 <- fit2_max
} else {
  cat("\n#### Falling back: Model 2 with (1 | participant) only ####\n")
  fit2 <- fit_or_load(stay ~ reward * K_c + (1 | participant), "fit2_fallback",
                      prior = priors[priors$class != "cor", ])
  print(summary(fit2))
}
saveRDS(fit2, file.path(FITS_DIR, "fit2_final.rds"))

#### IS THERE A THIRD LEVEL? (trials within blocks within people) ####

# Trials are collected in 6 blocks per participant, so the data could be treated
# as three-level. We test whether blocks carry any variance of their own.
cat("\n#### Three-level check: adding (1 | block_id) ####\n")
fit2_blk <- fit_or_load(stay ~ reward * K_c + (reward | participant) + (1 | block_id),
                        "fit2_block")
print(summary(fit2_blk))

#### MODEL COMPARISON: Model 1 vs Model 2 ####

cat("\n#### Model comparison ####\n")

## Trial-level LOO ---------------------------------------------------------
# NOTE: this leaves out one TRIAL at a time. A trial is cheap to predict from
# the other ~143 trials of the same person, so leave-one-trial-out CV is close
# to blind to a between-person moderator. We report it, but it is not the
# comparison the research question calls for.
fit1 <- add_criterion(fit1, "loo")
fit2 <- add_criterion(fit2, "loo")
print(loo_compare(fit1, fit2))
cat("\nPareto-k diagnostics:\n")
print(pareto_k_table(loo(fit1)))
print(pareto_k_table(loo(fit2)))

## Participant-level K-fold CV ---------------------------------------------
# Whole participants are held out, so the folds actually test whether knowing a
# new person's K helps predict their choices. This is the comparison that
# matches the hypothesis.
kf1 <- kfold(fit1, folds = "grouped", group = "participant", K = 10, seed = SEED)
kf2 <- kfold(fit2, folds = "grouped", group = "participant", K = 10, seed = SEED)
print(loo_compare(kf1, kf2))

## Bayes factor for the interaction ----------------------------------------
# Savage-Dickey density ratio at gamma_11 = 0. This is well defined only because
# we set a proper prior on the slopes.
cat("\nSavage-Dickey Bayes factor for gamma_11 = 0:\n")
print(bayesfactor_parameters(fit2, null = 0, parameters = "b_reward:K_c"))

## R2 -----------------------------------------------------------------------
cat("\nBayesian R2 (marginal = fixed effects only; conditional = + random effects):\n")
print(bayes_R2(fit1, re_formula = NA)); print(bayes_R2(fit1))
print(bayes_R2(fit2, re_formula = NA)); print(bayes_R2(fit2))

## Pseudo-R2 ----------------------------------------------------------------
# The effect size that fits a cross-level interaction: a level 2 moderator of a
# level 1 slope should explain variance in that RANDOM SLOPE. So we ask by how
# much adding K_c shrinks tau_11 (and tau_00, for the K_c main effect).
# Computed per draw, so it comes with a credible interval.
d1 <- as_draws_df(fit1)
d2 <- as_draws_df(fit2)
cat("\nPseudo-R2 (proportional reduction in each variance component, M1 -> M2):\n")
post_summary(1 - d2$sd_participant__reward^2    / d1$sd_participant__reward^2,    "  random slope (tau_11)")
post_summary(1 - d2$sd_participant__Intercept^2 / d1$sd_participant__Intercept^2, "  random intercept (tau_00)")

## Pseudo-standardised coefficients -----------------------------------------
cat("\nPseudo-standardised fixed effects (Model 2):\n")
print(standardize_parameters(fit2, method = "pseudo"))

#### FIXED AND RANDOM EFFECTS OF THE KEY MODEL ####

cat("\n#### Model 2 posterior summaries ####\n")
post_summary(d2$b_Intercept,     "gamma_00 (intercept)")
post_summary(d2$b_reward,        "gamma_10 (reward)")
post_summary(d2$b_K_c,           "gamma_01 (K_c)")
post_summary(d2$`b_reward:K_c`,  "gamma_11 (reward x K_c)")
post_summary(d2$sd_participant__Intercept, "sqrt(tau_00)")
post_summary(d2$sd_participant__reward,    "sqrt(tau_11)")
post_summary(d2$cor_participant__Intercept__reward, "cor(U_0, U_1)")

#### FOLLOW-UP ANALYSES ####

k_sd <- sd(distinct(df, participant, K)$K) # between-SUBJECT SD, not the SD over trials

## Simple slopes: the reward effect at low / average / high WM capacity ------
cat("\nSimple slopes of reward (log-odds) at -1 SD, mean, +1 SD of K:\n")
print(comparisons(fit2, variables = "reward",
                  newdata = datagrid(K_c = c(-k_sd, 0, k_sd)),
                  re_formula = NA, type = "link"))

cat("\nSame slopes as odds ratios:\n")
print(comparisons(fit2, variables = "reward",
                  newdata = datagrid(K_c = c(-k_sd, 0, k_sd)),
                  re_formula = NA, comparison = "lnor", transform = "exp"))

cat("\nDifference between the +1 SD and -1 SD reward slopes (= gamma_11 * 2 SD):\n")
post_summary(d2$`b_reward:K_c` * 2 * k_sd, "  slope difference (log-odds)")

## Average marginal effect --------------------------------------------------
# Predictions at U = 0 describe the MEDIAN participant, not the average one:
# because the logit link is non-linear, inv_logit(mean(eta)) != mean(inv_logit(eta)).
# Averaging over the estimated random effects gives the population-average
# effect of reward on the probability scale.
cat("\nAverage marginal effect of reward on P(stay), averaged over participants:\n")
print(avg_comparisons(fit2, variables = "reward", re_formula = NULL))

cat("\nAverage (marginal) odds ratio for reward:\n")
print(avg_comparisons(fit2, variables = "reward", re_formula = NULL,
                      comparison = "lnor", transform = "exp"))

#### POSTERIOR PREDICTIVE CHECKS ####

cat("\n#### Posterior predictive checks written to analysis/figures/ ####\n")
dir.create("analysis/figures", showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave("analysis/figures/ppc_bars.png",
                pp_check(fit2, type = "bars", ndraws = 500),
                width = 5, height = 3.5, dpi = 300, bg = "white")
ggplot2::ggsave("analysis/figures/ppc_by_participant.png",
                pp_check(fit2, type = "stat_grouped", stat = "mean",
                         group = "participant", ndraws = 500),
                width = 8, height = 6, dpi = 300, bg = "white")

cat("\n#### Session info ####\n")
print(sessionInfo())

cat(sprintf("\nDone. Log written to %s, fitted models saved to %s/\n", LOG_FILE, FITS_DIR))
