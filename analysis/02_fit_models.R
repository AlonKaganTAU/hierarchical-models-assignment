rm(list = ls())

#### SETUP ####
library(dplyr)
library(brms)
library(loo)
library(posterior)
library(bayesplot)

DATA_CSV <- "data/trial_level.csv"
FITS_DIR <- "analysis/fits"
LOG_FILE <- "analysis/model_log.txt"
dir.create(FITS_DIR, showWarnings = FALSE, recursive = TRUE)

CHAINS <- 4
ITER   <- 3500
WARMUP <- 1000
CORES  <- 4
SEED   <- 20260902

# Weakly informative priors on the logit scale. normal(0, 1.5) on the
# regression coefficients keeps odds ratios within a plausible range (roughly
# 1/20 to 20) without pulling the posterior towards any effect; brms' own
# default for class "b" is an improper flat prior, which is not what "weakly
# informative" means. Models without predictors or without a random slope take
# the matching subset.
PRIORS_UNCONDITIONAL <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd")
)

MODEL_PRIORS <- c(
  prior(normal(0, 1.5), class = "b"),
  PRIORS_UNCONDITIONAL,
  prior(lkj(2), class = "cor")
)

# Everything printed below also goes to LOG_FILE. Any sink left open by an
# earlier failed run is closed first, so the log always starts clean.
while (sink.number() > 0) sink()
log_con <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)

df <- read.csv(DATA_CSV) |>
  mutate(participant = factor(participant))

subjects <- distinct(df, participant, K)
k_sd     <- sd(subjects$K)  # between-subject SD, the scale for the +-1 SD contrasts

cat(sprintf("Loaded %d trials from %d subjects. Between-subject SD of K = %.3f\n\n",
            nrow(df), nrow(subjects), k_sd))

# Fits (or, on a rerun, loads) one model and caches it, so an interrupted run
# does not have to recompile and resample everything from scratch.
fit_or_load <- function(formula, name, priors = MODEL_PRIORS) {
  path <- file.path(FITS_DIR, paste0(name, ".rds"))
  if (file.exists(path)) {
    cat(sprintf("Loading cached fit: %s\n", path))
    return(readRDS(path))
  }
  fit <- brm(formula, data = df, family = bernoulli(), prior = priors,
             chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
             seed = SEED, refresh = 0, backend = "rstan")
  fit <- add_criterion(fit, "loo")
  saveRDS(fit, path)
  fit
}

# A model has converged if there were no post-warmup divergent transitions and
# every Rhat is below 1.01; anything else falls back to a simpler random
# effects structure (Barr et al., 2013), and both attempts are kept.
model_has_converged <- function(fit) {
  divergences <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  diagnostics <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  cat(sprintf("  divergent transitions: %d | max Rhat: %.4f | min bulk ESS: %.0f | min tail ESS: %.0f\n",
              divergences, max(diagnostics$rhat, na.rm = TRUE),
              min(diagnostics$ess_bulk, na.rm = TRUE), min(diagnostics$ess_tail, na.rm = TRUE)))
  divergences == 0 && max(diagnostics$rhat, na.rm = TRUE) < 1.01
}

# Posterior median and 95% credible interval for a vector of draws. The
# probability of direction is only meaningful for quantities that can take
# either sign, so it is optional.
describe_draws <- function(draws, pd = TRUE) {
  described <- c(median = median(draws),
                 lower  = unname(quantile(draws, 0.025)),
                 upper  = unname(quantile(draws, 0.975)))
  if (pd) described <- c(described, pd = max(mean(draws > 0), mean(draws < 0)))
  round(described, 4)
}

#### MODEL 0 - UNCONDITIONAL MODEL AND THE ICC ####
cat("\n#### Model 0: stay ~ 1 + (1 | participant) ####\n")
fit0 <- fit_or_load(stay ~ 1 + (1 | participant), "fit0", PRIORS_UNCONDITIONAL)
print(summary(fit0))
invisible(model_has_converged(fit0))

# In a logistic HLM the level-1 variance is not estimated but fixed at pi^2/3
# on the latent (logit) scale, which is what the ICC is computed against.
tau00_0   <- as_draws_df(fit0)$sd_participant__Intercept^2
icc_draws <- tau00_0 / (tau00_0 + pi^2 / 3)
cat("\nUnconditional ICC (latent-threshold):\n")
print(describe_draws(icc_draws, pd = FALSE))

# plogis(gamma_00) is a generalised mean, not the average probability of
# staying: the inverse logit of a mean is not the mean of inverse logits. The
# average of the model's per-trial predicted probabilities is the unbiased
# quantity to compare against the raw stay rate.
cat(sprintf("\nplogis(gamma_00) = %.3f | average predicted P(stay) = %.3f | observed stay rate = %.3f\n",
            plogis(median(as_draws_df(fit0)$b_Intercept)),
            mean(posterior_epred(fit0)), mean(df$stay)))

# The same logic applied to the time-varying predictor: how much of the
# variance in reward_oneback is between-person? If this is near zero the
# predictor carries essentially only within-person information, and its
# within- and between-person effects cannot be smushed.
cat("\n#### Model P: reward_oneback ~ 1 + (1 | participant) (predictor ICC) ####\n")
fit_pred <- fit_or_load(reward_oneback ~ 1 + (1 | participant), "fit_predictor_icc",
                        PRIORS_UNCONDITIONAL)
tau_pred      <- as_draws_df(fit_pred)$sd_participant__Intercept^2
icc_pred_draws <- tau_pred / (tau_pred + pi^2 / 3)
cat("\nICC of reward_oneback (latent-threshold):\n")
print(describe_draws(icc_pred_draws, pd = FALSE))

#### MODEL 1 - LEVEL-1 PREDICTOR (REWARD) ####
cat("\n#### Model 1: stay ~ reward_oneback + (reward_oneback | participant) ####\n")
fit1_max <- fit_or_load(stay ~ reward_oneback + (reward_oneback | participant), "fit1_maximal")
print(summary(fit1_max))

fit1_converged <- model_has_converged(fit1_max)
if (fit1_converged) {
  fit1 <- fit1_max
} else {
  cat("\n#### Falling back: Model 1 with (1 | participant) only ####\n")
  fit1 <- fit_or_load(stay ~ reward_oneback + (1 | participant), "fit1_random_intercept",
                      c(prior(normal(0, 1.5), class = "b"), PRIORS_UNCONDITIONAL))
  print(summary(fit1))
}

# Is the random reward slope itself warranted? Compare against the same fixed
# effects with random intercepts only.
cat("\n#### Model 1r (random intercepts only), for the random-slope comparison ####\n")
fit1_ri <- fit_or_load(stay ~ reward_oneback + (1 | participant), "fit1_random_intercept",
                       c(prior(normal(0, 1.5), class = "b"), PRIORS_UNCONDITIONAL))
cat("\nLOO comparison, random slope vs random intercept only:\n")
print(loo_compare(fit1, fit1_ri))

#### MODEL 2 - ADDING THE LEVEL-2 PREDICTOR (K) ####
cat("\n#### Model 2: stay ~ reward_oneback + K_c + (reward_oneback | participant) ####\n")
fit2 <- fit_or_load(stay ~ reward_oneback + K_c + (reward_oneback | participant), "fit2_main_effects")
print(summary(fit2))
invisible(model_has_converged(fit2))

#### MODEL 3 - CROSS-LEVEL INTERACTION (KEY MODEL) ####
cat("\n#### Model 3: stay ~ reward_oneback * K_c + (reward_oneback | participant) ####\n")
fit3_max <- fit_or_load(stay ~ reward_oneback * K_c + (reward_oneback | participant), "fit3_maximal")
print(summary(fit3_max))

fit3_converged <- model_has_converged(fit3_max)
if (fit3_converged) {
  fit3 <- fit3_max
} else {
  cat("\n#### Falling back: Model 3 with (1 | participant) only ####\n")
  fit3 <- fit_or_load(stay ~ reward_oneback * K_c + (1 | participant), "fit3_fallback",
                      c(prior(normal(0, 1.5), class = "b"), PRIORS_UNCONDITIONAL))
  print(summary(fit3))
}
saveRDS(fit3, file.path(FITS_DIR, "fit3_final.rds"))

#### MODEL COMPARISONS ####
cat("\n#### Model comparisons (approximate leave-one-out cross-validation) ####\n")
cat("\nModel 1 vs Model 2 (adding the level-2 main effect of K):\n")
print(loo_compare(fit1, fit2))
cat("\nModel 2 vs Model 3 (adding the cross-level interaction):\n")
print(loo_compare(fit2, fit3))

cat("\nBayes R2 (marginal = fixed effects only; conditional = fixed + random):\n")
r2_table <- rbind(
  bayes_R2(fit1, re_formula = NA), bayes_R2(fit1),
  bayes_R2(fit2, re_formula = NA), bayes_R2(fit2),
  bayes_R2(fit3, re_formula = NA), bayes_R2(fit3)
)
rownames(r2_table) <- c("M1 marginal", "M1 conditional", "M2 marginal",
                        "M2 conditional", "M3 marginal", "M3 conditional")
print(round(r2_table, 4))

#### PSEUDO-R2 FOR THE VARIANCE COMPONENT EACH PREDICTOR TARGETS ####
# A level-2 main effect should explain between-person variance in the
# intercept (tau00); a cross-level interaction should explain between-person
# variance in the level-1 slope it moderates (tau11). Draws from the two
# models are independent, so pairing them gives a Monte Carlo approximation of
# the distribution of the proportion reduction.
tau00_1 <- as_draws_df(fit1)$sd_participant__Intercept^2
tau00_2 <- as_draws_df(fit2)$sd_participant__Intercept^2
tau11_2 <- as_draws_df(fit2)$sd_participant__reward_oneback^2
tau11_3 <- as_draws_df(fit3)$sd_participant__reward_oneback^2

# The point estimate is the class formula applied to the posterior medians of
# the two variance components; the draw-by-draw interval is also printed, but
# because the two models are estimated independently it is very wide and is
# only a rough indication of how uncertain these numbers are.
cat(sprintf("\nPseudo-R2 for K_c on tau00 (Model 1 -> Model 2): %.3f\n",
            1 - median(tau00_2) / median(tau00_1)))
print(describe_draws((tau00_1 - tau00_2) / tau00_1, pd = FALSE))

cat(sprintf("\nPseudo-R2 for reward_oneback x K_c on tau11 (Model 2 -> Model 3): %.3f\n",
            1 - median(tau11_3) / median(tau11_2)))
print(describe_draws((tau11_2 - tau11_3) / tau11_2, pd = FALSE))

#### PSEUDO-STANDARDISED COEFFICIENTS ####
# Level-1 effects are scaled by the within-person SDs and level-2 effects by
# the between-person SDs. On the latent logit scale the within-person outcome
# SD is fixed at pi/sqrt(3); the between-person outcome SD is sqrt(tau00) from
# the unconditional model, whose draws are paired with Model 3's in the same
# Monte Carlo sense as the pseudo-R2 above.
draws3     <- as_draws_df(fit3)
sd_within  <- pi / sqrt(3)
sd_between <- sqrt(tau00_0)

cat("\nPseudo-beta, within-person effect of reward_oneback:\n")
print(describe_draws(draws3$b_reward_oneback * sd(df$reward_oneback_wp) / sd_within))
cat("\nPseudo-beta, between-person effect of K:\n")
print(describe_draws(draws3$b_K_c * k_sd / sd_between))

#### FIXED EFFECTS OF THE KEY MODEL, ON THE ODDS-RATIO SCALE TOO ####
cat("\nModel 3 fixed effects, log-odds and odds ratios:\n")
log_odds <- rbind(
  Intercept            = describe_draws(draws3$b_Intercept),
  reward_oneback       = describe_draws(draws3$b_reward_oneback),
  K_c                  = describe_draws(draws3$b_K_c),
  `reward_oneback:K_c` = describe_draws(draws3$`b_reward_oneback:K_c`)
)
odds_ratios <- rbind(
  Intercept            = describe_draws(exp(draws3$b_Intercept), pd = FALSE),
  reward_oneback       = describe_draws(exp(draws3$b_reward_oneback), pd = FALSE),
  K_c                  = describe_draws(exp(draws3$b_K_c), pd = FALSE),
  `reward_oneback:K_c` = describe_draws(exp(draws3$`b_reward_oneback:K_c`), pd = FALSE)
)
colnames(odds_ratios) <- paste0("OR_", colnames(odds_ratios))
print(cbind(log_odds, odds_ratios))

#### FOLLOW-UP: SIMPLE SLOPES OF REWARD AT +-1 SD OF K ####
# Population-level (fixed effects only) predictions, mirroring Figure 1.
for (k_level in c(-k_sd, k_sd)) {
  slope   <- draws3$b_reward_oneback + draws3$`b_reward_oneback:K_c` * k_level
  p_lose  <- plogis(draws3$b_Intercept + draws3$b_K_c * k_level)
  p_win   <- plogis(draws3$b_Intercept + draws3$b_K_c * k_level + slope)

  cat(sprintf("\nK_c = %+.2f (%s1 SD of K):\n", k_level, ifelse(k_level < 0, "-", "+")))
  cat("  simple slope (log-odds): "); print(describe_draws(slope))
  cat("  simple slope (OR):       "); print(describe_draws(exp(slope), pd = FALSE))
  cat("  P(stay | no reward):     "); print(describe_draws(p_lose, pd = FALSE))
  cat("  P(stay | reward):        "); print(describe_draws(p_win, pd = FALSE))
  cat("  difference in P(stay):   "); print(describe_draws(p_win - p_lose))
}

# The interaction expressed on the probability scale: how much smaller is the
# reward effect at +1 SD of K than at -1 SD?
slope_diff <- (plogis(draws3$b_Intercept + draws3$b_K_c * k_sd + draws3$b_reward_oneback + draws3$`b_reward_oneback:K_c` * k_sd) -
               plogis(draws3$b_Intercept + draws3$b_K_c * k_sd)) -
              (plogis(draws3$b_Intercept - draws3$b_K_c * k_sd + draws3$b_reward_oneback - draws3$`b_reward_oneback:K_c` * k_sd) -
               plogis(draws3$b_Intercept - draws3$b_K_c * k_sd))
cat("\nDifference between the two simple reward effects, probability scale:\n")
print(describe_draws(slope_diff))

#### ROBUSTNESS CHECKS ####
# (a) Trials are also nested in blocks; does a block-within-person random
#     intercept change the interaction?
cat("\n#### Robustness (a): adding a block-within-person random intercept ####\n")
fit3_block <- fit_or_load(
  stay ~ reward_oneback * K_c + (reward_oneback | participant) + (1 | participant:block_number),
  "fit3_block"
)
print(summary(fit3_block))
cat("\nInteraction with the block level added:\n")
print(describe_draws(as_draws_df(fit3_block)$`b_reward_oneback:K_c`))

# (b) Separating the within- and between-person parts of the time-varying
#     predictor, so that the level-1 effect is purely within-person.
cat("\n#### Robustness (b): within- and between-person reward separated ####\n")
fit3_split <- fit_or_load(
  stay ~ reward_oneback_wp * K_c + reward_oneback_pm_c + (reward_oneback_wp | participant),
  "fit3_within_between"
)
print(summary(fit3_split))
cat("\nInteraction with reward split into its within- and between-person parts:\n")
print(describe_draws(as_draws_df(fit3_split)$`b_reward_oneback_wp:K_c`))

#### MODEL CHECKING ####
# A logistic HLM has no level-1 residual variance to inspect, so the checks
# that matter are whether the model reproduces the data and whether the random
# effects look normal (assumptions covered in class).
cat("\n#### Model checking (Model 3) ####\n")

replicated <- posterior_predict(fit3, ndraws = 1000)

cat(sprintf("Overall stay rate: observed %.3f, predicted %.3f [%.3f, %.3f]\n",
            mean(df$stay), median(rowMeans(replicated)),
            quantile(rowMeans(replicated), 0.025), quantile(rowMeans(replicated), 0.975)))

observed_by_subject   <- tapply(df$stay, df$participant, mean)
replicated_by_subject <- apply(replicated, 1, \(y) tapply(y, df$participant, mean))
predictive_interval   <- apply(replicated_by_subject, 1, quantile, c(0.025, 0.975))
covered <- observed_by_subject >= predictive_interval[1, ] &
           observed_by_subject <= predictive_interval[2, ]
cat(sprintf("Per-subject stay rates inside their 95%% predictive interval: %d of %d\n",
            sum(covered), length(covered)))

random_effects <- ranef(fit3)$participant[, "Estimate", ]
cat("\nShapiro-Wilk on the posterior mean random effects (normality assumption):\n")
print(shapiro.test(random_effects[, "Intercept"]))
print(shapiro.test(random_effects[, "reward_oneback"]))

#### SAVE A COMPACT RESULTS OBJECT ####
results <- list(
  k_sd            = k_sd,
  icc             = describe_draws(icc_draws, pd = FALSE),
  icc_predictor   = describe_draws(icc_pred_draws, pd = FALSE),
  r2              = r2_table,
  pseudo_r2_tau00 = 1 - median(tau00_2) / median(tau00_1),
  pseudo_r2_tau11 = 1 - median(tau11_3) / median(tau11_2),
  interaction     = describe_draws(draws3$`b_reward_oneback:K_c`),
  loo_1_2         = loo_compare(fit1, fit2),
  loo_2_3         = loo_compare(fit2, fit3)
)
saveRDS(results, file.path(FITS_DIR, "results_summary.rds"))

sink()
cat(sprintf("\nDone. Log written to %s, fitted models saved to %s/\n", LOG_FILE, FITS_DIR))
