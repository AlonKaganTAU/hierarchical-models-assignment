# Does working memory capacity constrain automatic reward-driven choice under known non-contingency?

Final assignment, Hierarchical Models (1071-8710), Tel Aviv University.

The question is whether working memory capacity (Cowan's K, a level-2 trait) moderates
the trial-level win-stay lose-shift effect (level 1) in an explicitly uncontrollable
three-armed bandit — a cross-level interaction, fit as a single Bayesian logistic HLM
(`stay ~ reward_oneback * K_c + (reward_oneback | participant)`, brms).

## Repository layout

- `data/raw/` — raw per-subject CSVs, task and working-memory change-detection, fourth
  pilot, N = 35. Source: `illusion-of-control` repo, `data/ioc-all-fixed-pilot/`.
- `data/trial_level.csv` — processed trial-level data, one row per usable trial.
- `analysis/01_build_trial_level_data.R` — builds `data/trial_level.csv`: Cowan's K per
  subject, trial pairing and validity, and centring of each variable at its own level.
- `analysis/02_fit_models.R` — fits Models 0–3 plus the two robustness models, and
  computes the ICC, model comparisons, Bayes R², pseudo-R², pseudo-standardised
  coefficients, simple slopes and model checks. Writes `analysis/model_log.txt`.
- `analysis/03_make_plot.R` — Figure 1 of the write-up.
- `analysis/model_log.txt` — full output of `02`; every number in the write-up comes
  from this file.
- `writeup/` — the submitted Word document.

## Models

| | formula |
|---|---|
| Model 0 | `stay ~ 1 + (1 \| participant)` — unconditional, for the ICC |
| Model 1 | `stay ~ reward_oneback + (reward_oneback \| participant)` |
| Model 2 | `stay ~ reward_oneback + K_c + (reward_oneback \| participant)` |
| Model 3 | `stay ~ reward_oneback * K_c + (reward_oneback \| participant)` — key model |

Two robustness models add a block-within-person random intercept, and split
`reward_oneback` into its within- and between-person parts.

## Reproducing

Run the three scripts in order from the repository root:

```r
source("analysis/01_build_trial_level_data.R")
source("analysis/02_fit_models.R")
source("analysis/03_make_plot.R")
```

Requires `dplyr`, `tidyr`, `readr`, `stringr`, `ggplot2`, `patchwork`, `brms`, `loo`,
`posterior` and `bayesplot`, with the rstan backend. Sampling takes roughly half an hour
on four cores; fitted models are cached in `analysis/fits/` (not tracked in git), so a
rerun reuses whatever is already there. `brm()` is called with a fixed `seed`, so a fresh
run reproduces the numbers in the write-up.
