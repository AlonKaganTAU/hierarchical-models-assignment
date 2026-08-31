rm(list = ls())

#### SETUP ####
library(dplyr)
library(readr)
library(purrr)
library(stringr)

DATA_TASK <- "data/raw/task"
DATA_WM   <- "data/raw/wm"
OUT_CSV   <- "data/trial_level.csv"

# Source: illusion-of-control repo, analysis/wm_wsls/wm_wsls.R (fourth pilot,
# fixed inter-trial interval, N = 35). This script reuses that pipeline's WM
# capacity computation and trial-pairing logic, but keeps all 35 subjects (no
# understood/felt exclusion) and returns trial-level data rather than
# per-subject WSLS betas: the assignment models the trial-level DV directly,
# with reward x K as a cross-level interaction.
#
# Output columns, by level of the hierarchy:
#   participant   - the random grouping variable
#   block_number  - 6 blocks per participant (level 1.5; see 02 for the
#                   three-level check that justifies ignoring it)
#   trial_number  - position within block
#   reward (L1)   - outcome of the PREVIOUS trial, 0/1 (= reward_oneback in the
#                   raw files); this is the level 1 predictor
#   stay   (L1)   - whether the same machine was chosen again, 0/1 (the DV)
#   K      (L2)   - working memory capacity, one value per participant
#   K_c    (L2)   - K grand-mean centred on the analysed sample

#### WM CAPACITY (K) PER SUBJECT ####
# Cowan's K = set size * (hit rate + correct rejection rate - 1), where the
# "different" trials give hits and the "same" trials give correct rejections.
# K is averaged over the two experimental set sizes (4 and 8).
wm_files <- list.files(DATA_WM, pattern = "^ioc-wm_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

wm <- wm_files |>
  map(\(f) {
    pid <- str_extract(f, "[a-f0-9]{24}")
    exp <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
      filter(trial_name == "test_squares", block_type == "exp") |>
      mutate(set_size = as.numeric(set_size),
             acc_bool = tolower(accuracy) == "true")

    k_by_set_size <- map_dbl(c(4, 8), \(ss) {
      trials    <- exp |> filter(set_size == ss)
      different <- trials |> filter(condition == "d")
      same      <- trials |> filter(condition == "s")
      if (nrow(different) == 0 || nrow(same) == 0) return(NA_real_)
      ss * (mean(different$acc_bool) + mean(same$acc_bool) - 1)
    })
    tibble(participant = pid, K = mean(k_by_set_size, na.rm = TRUE))
  }) |>
  list_rbind()

#### TRIAL-LEVEL STAY / REWARD PAIRS ####
task_files <- list.files(DATA_TASK, pattern = "^ioc-all_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

df_raw <- task_files |>
  map(\(f) {
    pid <- str_extract(f, "[a-f0-9]{24}")
    read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
      mutate(participant = pid)
  }) |>
  list_rbind()

# Trials are paired within block, so the first trial of every block drops out
# (it has no preceding trial). A pair is kept only when BOTH trials carried a
# valid, non-timed-out response - otherwise "stay" is undefined.
df <- df_raw |>
  filter(task == "gambling_choice", block_number != "training") |>
  mutate(trial_number    = as.numeric(trial_number),
         reward          = as.numeric(reward),
         reward_oneback  = as.numeric(reward_oneback),
         is_choice_valid = as.logical(is_choice_valid)) |>
  arrange(participant, block_number, trial_number) |>
  group_by(participant, block_number) |>
  mutate(choice_prev = lag(choice_key),
         valid_prev  = lag(is_choice_valid),
         reward_prev = lag(reward)) |>
  ungroup()

# The raw files already carry choice_key_oneback / reward_oneback. We derive the
# lags ourselves (so the pairing rule is explicit and auditable) and then check
# the two agree - if they ever disagree, the trial ordering assumption is wrong.
lag_check <- df |>
  filter(is_choice_valid, valid_prev) |>
  summarise(reward_mismatch = sum(reward_prev != reward_oneback, na.rm = TRUE),
            choice_mismatch = sum(choice_prev != choice_key_oneback, na.rm = TRUE))
cat(sprintf("Lag check vs. raw *_oneback columns: %d reward mismatches, %d choice mismatches\n",
            lag_check$reward_mismatch, lag_check$choice_mismatch))

# Staying is only meaningful if the previous machine is still on offer. In this
# pilot no machine was ever withheld, but we check rather than assume.
cat(sprintf("Trials where the previous choice was unavailable: %d\n",
            sum(df$unavailable_keys != "[]", na.rm = TRUE)))

df <- df |>
  filter(is_choice_valid, valid_prev, !is.na(reward_prev)) |>
  mutate(stay = as.integer(choice_key == choice_prev)) |>
  select(participant, block_number, trial_number, reward = reward_prev, stay)

#### MERGE WITH WM CAPACITY ####
# One WM file has no matching task file, so the inner join drops it. K must
# therefore be centred AFTER the join - centring on the raw WM table would
# centre on a person who is not in the analysis, and the fixed intercept would
# no longer sit at the analysed sample's mean K.
df <- df |>
  inner_join(wm, by = "participant") |>
  filter(!is.na(K))

K_by_subject <- df |> distinct(participant, K)
df <- df |> mutate(K_c = K - mean(K_by_subject$K))

#### SANITY CHECKS ####
n_subjects <- n_distinct(df$participant)
cat(sprintf("WM files: %d; task files: %d; subjects in analysis: %d\n",
            nrow(wm), length(task_files), n_subjects))

trials_per_subject <- df |> count(participant, name = "n_trials")
cat(sprintf("Trials per subject: min = %d, max = %d, mean = %.1f\n",
            min(trials_per_subject$n_trials), max(trials_per_subject$n_trials),
            mean(trials_per_subject$n_trials)))

cat(sprintf("K: mean = %.2f, SD = %.2f, range %.2f to %.2f (%d subjects below 0)\n",
            mean(K_by_subject$K), sd(K_by_subject$K),
            min(K_by_subject$K), max(K_by_subject$K), sum(K_by_subject$K < 0)))

cat(sprintf("Reward rate: %.3f overall; between-subject range %.3f to %.3f\n",
            mean(df$reward),
            min(tapply(df$reward, df$participant, mean)),
            max(tapply(df$reward, df$participant, mean))))

cat(sprintf("Missing reward/stay in final data: %d\n", sum(is.na(df$reward) | is.na(df$stay))))

#### WRITE OUTPUT ####
write_csv(df, OUT_CSV)
cat(sprintf("\nSaved -> %s (%d rows, %d subjects)\n", OUT_CSV, nrow(df), n_subjects))
