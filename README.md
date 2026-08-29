# Working Memory Capacity Moderates Automatic Reward-Driven Choice Under Known Non-Contingency

Final assignment, Hierarchical Models (Tel Aviv University). Tests whether working
memory capacity (K) moderates the win-stay lose-shift (WSLS) effect in an explicitly
uncontrollable three-armed bandit task, via a single Bayesian cross-level interaction
model (`stay ~ reward * K_c + (reward | participant)`, family = bernoulli, brms).

## Repository layout

- `data/raw/` — raw per-subject CSVs (task + working-memory change-detection), fourth
  pilot, N = 35. Source: `illusion-of-control` repo, `data/ioc-all-fixed-pilot/`.
- `data/trial_level.csv` — processed trial-level data (participant, block_number,
  trial_number, reward, stay, K, K_c), all 35 subjects, no exclusions applied.
- `analysis/01_build_trial_level_data.R` — builds `data/trial_level.csv` from the raw
  data (WM capacity via Cowan's K, trial pairing/validity via `dplyr::lag`).
- `analysis/02_fit_models.R` — fits Models 0-2 in brms (unconditional/ICC, reward-only,
  reward x K_c), with the maximal-random-effects-then-fallback convergence logic,
  `loo_compare`, and `bayes_R2`. Writes `analysis/model_log.txt` and caches fitted
  models to `analysis/fits/` (not tracked in git — rerun the script to regenerate).
- `analysis/03_make_plot.R` — simple-slopes figure used in the write-up.
- `writeup/` — the submitted Word document.

## Reproducing

Run the three scripts in order from the repository root:

```r
source("analysis/01_build_trial_level_data.R")
source("analysis/02_fit_models.R")
source("analysis/03_make_plot.R")
```

Requires R with `dplyr`, `tidyr`, `readr`, `purrr`, `stringr`, `ggplot2`, `brms`,
`loo`, `posterior`, `bayesplot` (rstan backend). Model fitting takes a few minutes on
4 cores; a rerun will reuse any cached fits already present in `analysis/fits/`.
