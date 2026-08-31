# Working Memory Capacity and Automatic Reward-Driven Choice Under Known Non-Contingency

Final assignment, Hierarchical Models (Tel Aviv University). Tests whether working
memory capacity (K) moderates the win-stay lose-shift (WSLS) effect in an explicitly
uncontrollable three-armed bandit task, via a Bayesian cross-level interaction model
(`stay ~ reward * K_c + (reward | participant)`, family = bernoulli, brms).

## Repository layout

- `data/raw/` — raw per-subject CSVs (task + working-memory change-detection), fourth
  pilot, N = 35. Source: `illusion-of-control` repo, `data/ioc-all-fixed-pilot/`.
- `data/trial_level.csv` — processed trial-level data (participant, block_number,
  trial_number, reward, stay, K, K_c), all 35 subjects, no participant exclusions.
- `analysis/01_build_trial_level_data.R` — builds `data/trial_level.csv` from the raw
  data: working memory capacity via Cowan's K, trial pairing via `dplyr::lag` (checked
  against the raw `*_oneback` columns), and grand-mean centring of K on the analysed
  sample.
- `analysis/02_fit_models.R` — fits Models 0–2 in brms (unconditional/ICC, reward-only,
  reward × K_c) with explicit weakly-informative priors and the
  maximal-random-effects-then-fallback convergence logic, plus the unconditional ICC,
  the ICC of the level 1 predictor, a three-level (blocks within people) check,
  trial-level `loo` and participant-level `kfold`, a Savage–Dickey Bayes factor,
  `bayes_R2`, pseudo-R², pseudo-standardised coefficients, simple slopes, average
  marginal effects, and posterior predictive checks. Writes `analysis/model_log.txt`
  and caches fitted models to `analysis/fits/` (not tracked in git — rerun to regenerate).
- `analysis/03_make_plot.R` — the two-panel figure used in the write-up.
- `writeup/` — the submitted Word document.

## Reproducing

Run the three scripts in order from the repository root:

```r
source("analysis/01_build_trial_level_data.R")
source("analysis/02_fit_models.R")
source("analysis/03_make_plot.R")
```

Requires R with `dplyr`, `readr`, `purrr`, `stringr`, `tidyr`, `ggplot2`, `patchwork`,
`brms`, `loo`, `posterior`, `bayesplot`, `bayestestR`, `marginaleffects`, `tidybayes`
(rstan backend). Sampling is seeded (`SEED <- 20260829`), so the reported numbers
reproduce exactly. Model fitting takes a few minutes on 4 cores; the participant-level
`kfold` refits each model 10 times and takes considerably longer. A rerun reuses any
cached fits in `analysis/fits/` whose formula still matches.
