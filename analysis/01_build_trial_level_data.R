rm(list = ls())

#### SETUP ####
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(datawizard)

DATA_TASK <- "data/raw/task"
DATA_WM   <- "data/raw/wm"
OUT_CSV   <- "data/trial_level.csv"

# Source: illusion-of-control repo, analysis/wm_wsls/wm_wsls.R (fourth pilot,
# fixed inter-trial interval, N = 35). This script reuses that pipeline's WM
# K computation and trial-pairing logic, but keeps all 35 subjects (no
# understood/felt exclusion) and returns trial-level data rather than
# per-subject WSLS betas, since the assignment models the trial-level DV
# directly with reward_oneback x K as a cross-level interaction.

#### WM CAPACITY (K) PER SUBJECT ####
wm_files <- list.files(DATA_WM, pattern = "^ioc-wm_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

wm <- do.call(rbind, lapply(wm_files, function(f) {
  pid <- str_extract(f, "[a-f0-9]{24}")
  exp <- read.csv(f, colClasses = "character") |>
    filter(trial_name == "test_squares", block_type == "exp") |>
    mutate(set_size = as.numeric(set_size),
           acc_bool = tolower(accuracy) == "true")

  k_by_set_size <- sapply(c(4, 8), function(ss) {
    t <- exp |> filter(set_size == ss)
    d <- t |> filter(condition == "d")
    s <- t |> filter(condition == "s")
    if (nrow(d) == 0 || nrow(s) == 0) return(NA_real_)
    ss * (mean(d$acc_bool) + mean(s$acc_bool) - 1)
  })
  tibble(participant = pid, K = mean(k_by_set_size, na.rm = TRUE))
}))

#### TRIAL-LEVEL STAY / REWARD PAIRS ####
task_files <- list.files(DATA_TASK, pattern = "^ioc-all_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

df_raw <- do.call(rbind, lapply(task_files, function(f) {
  pid <- str_extract(f, "[a-f0-9]{24}")
  read.csv(f, colClasses = "character") |>
    mutate(participant = pid)
}))

df <- df_raw |>
  filter(task == "gambling_choice", block_number != "training") |>
  mutate(trial_number    = as.numeric(trial_number),
         reward          = as.numeric(reward),
         is_choice_valid = as.logical(is_choice_valid)) |>
  arrange(participant, block_number, trial_number) |>
  group_by(participant, block_number) |>
  mutate(choice_prev       = lag(choice_key),
         valid_prev        = lag(is_choice_valid),
         reward_prev       = lag(reward)) |>
  ungroup() |>
  filter(is_choice_valid, valid_prev, !is.na(reward_prev)) |>
  mutate(stay = as.integer(choice_key == choice_prev))

#### VALIDATE LAG-DERIVED VALUES AGAINST THE RAW *_oneback COLUMNS ####
# The raw files also ship choice_key_oneback/reward_oneback columns, but these
# run across block boundaries (e.g. block 1 trial 1's reward_oneback is
# actually carried over from the training block), so they are not what we
# want as the predictor -- our own group_by(block)-then-lag() correctly
# resets to NA at each block's first trial. We only use the raw columns here
# to sanity-check agreement on the trials where both are defined (i.e. not
# each block's first trial, which we already dropped above).
choice_mismatch <- sum(df$choice_prev != df$choice_key_oneback, na.rm = TRUE)
reward_mismatch <- sum(df$reward_prev != as.numeric(df$reward_oneback), na.rm = TRUE)
cat(sprintf("Lag-derived vs. raw *_oneback mismatches: choice = %d, reward = %d (of %d rows)\n",
            choice_mismatch, reward_mismatch, nrow(df)))

df <- df |>
  select(participant, block_number, trial_number, reward_oneback = reward_prev, stay)

#### MERGE WITH WM CAPACITY ####
# Centre K on the analysis sample (one row per participant), not on the
# trial-level frame -- K is repeated once per trial, so centring after the
# join to df would weight the mean by each participant's trial count rather
# than treating each participant equally.
subjects <- tibble(participant = unique(df$participant)) |>
  inner_join(wm, by = "participant") |>
  filter(!is.na(K)) |>
  mutate(K_c = center(K))

df <- df |>
  inner_join(subjects, by = "participant")

#### SANITY CHECKS ####
n_subjects <- n_distinct(df$participant)
cat(sprintf("Subjects: %d\n", n_subjects))

trials_per_subject <- df |> count(participant, name = "n_trials")
cat(sprintf("Trials per subject: min = %d, max = %d, mean = %.1f\n",
            min(trials_per_subject$n_trials), max(trials_per_subject$n_trials), mean(trials_per_subject$n_trials)))

cat(sprintf("K (analysis sample, one value per subject): range %.2f to %.2f, mean %.2f, SD %.2f\n",
            min(subjects$K), max(subjects$K), mean(subjects$K), sd(subjects$K)))
cat(sprintf("K_c mean over subjects (should be 0): %.12f\n", mean(subjects$K_c)))

cat(sprintf("Missing reward_oneback/stay in final data: %d\n", sum(is.na(df$reward_oneback) | is.na(df$stay))))

#### WRITE OUTPUT ####
write_csv(df, OUT_CSV)
cat(sprintf("\nSaved -> %s (%d rows, %d subjects)\n", OUT_CSV, nrow(df), n_subjects))
