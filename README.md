# Working Memory Capacity Moderates Automatic Reward-Driven Choice Under Known Non-Contingency

Final assignment, Hierarchical Models (Tel Aviv University). Tests whether working
memory capacity (K) moderates the win-stay lose-shift (WSLS) effect in an explicitly
uncontrollable three-armed bandit task, via a cross-level interaction model
(`stay ~ reward_oneback * K_c + (reward_oneback | participant)`, family =
binomial(link = "logit"), fitted with `glmer()`).

## Repository layout

- `data/raw/` — raw per-subject CSVs (task + working-memory change-detection), fourth
  pilot, N = 35. Source: `illusion-of-control` repo, `data/ioc-all-fixed-pilot/`.
- `data/trial_level.csv` — processed trial-level data (participant, block_number,
  trial_number, reward_oneback, stay, K, K_c), all 35 subjects, no exclusions applied.
- `analysis/01_build_trial_level_data.R` — builds `data/trial_level.csv` from the raw
  data (WM capacity via Cowan's K, trial pairing/validity via `dplyr::lag`, checked
  against the raw `*_oneback` columns). Centres K_c with `datawizard::center()` on
  the analysis sample (one row per participant, joined back onto the trial-level
  data), not on the raw WM-file sample or the trial-level frame.
- `analysis/02_fit_models.R` — fits Models 0-2 with `glmer()` (unconditional/ICC,
  reward_oneback-only, reward_oneback x K_c), maximal random-effects structure with
  the taught convergence remedies (control tweaks, optimizer, dropping only the
  random correlation via `||`), model comparison by likelihood-ratio test
  (`anova()`), and `r2_pseudo()` for the random-effect-variance effect size. Also
  reports reward_oneback's own (empty-model) ICC, justifying why it enters the
  models unsplit into within-/between-person components, and a block-within-
  participant robustness check (`(1 | participant:block_number)`, week 7 syntax).
  Writes `analysis/model_log.txt` and caches fitted models to `analysis/fits/` (not
  tracked in git — rerun the script to regenerate).
- `analysis/03_make_plot.R` — simple-slopes figure (`marginaleffects::plot_predictions()`,
  chaining ggplot layers directly onto its output) used in the write-up. K's SD for
  the ±1 SD probing values is computed over the 35 participants, not the 5,025 trials.
- `writeup/` — the submitted Word document.

## Reproducing

Run the three scripts in order from the repository root:

```r
source("analysis/01_build_trial_level_data.R")
source("analysis/02_fit_models.R")
source("analysis/03_make_plot.R")
```

Requires R with `dplyr`, `tidyr`, `readr`, `stringr`, `datawizard`, `ggplot2`, `lme4`,
`lmerTest`, `performance`, `parameters`, `marginaleffects`. `02_fit_models.R` also
sources the course's `r2_pseudo()` helper directly from GitHub (requires network
access) — see the `source()` call at the top of the script.
