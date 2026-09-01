rm(list = ls())

#### SETUP ####
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)

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

wm <- wm_files |>
  map(\(f) {
    pid <- str_extract(f, "[a-f0-9]{24}")
    exp <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
      filter(trial_name == "test_squares", block_type == "exp") |>
      mutate(set_size = as.numeric(set_size),
             acc_bool = tolower(accuracy) == "true")

    k_by_set_size <- map_dbl(c(4, 8), \(ss) {
      t <- exp |> filter(set_size == ss)
      d <- t |> filter(condition == "d")
      s <- t |> filter(condition == "s")
      if (nrow(d) == 0 || nrow(s) == 0) return(NA_real_)
      ss * (mean(d$acc_bool) + mean(s$acc_bool) - 1)
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

df <- df_raw |>
  filter(task == "gambling_choice", block_number != "training") |>
  mutate(trial_number   = as.numeric(trial_number),
         reward         = as.numeric(reward),
         is_choice_valid = as.logical(is_choice_valid)) |>
  arrange(participant, block_number, trial_number) |>
  group_by(participant, block_number) |>
  mutate(choice_prev       = lag(choice_key),
         valid_prev        = lag(is_choice_valid),
         reward_prev       = lag(reward)) |>
  ungroup() |>
  filter(is_choice_valid, valid_prev, !is.na(reward_prev)) |>
  mutate(stay = as.integer(choice_key == choice_prev)) |>
  select(participant, block_number, trial_number, reward_oneback = reward_prev, stay)

#### MERGE WITH WM CAPACITY ####
df <- df |>
  inner_join(wm, by = "participant") |>
  filter(!is.na(K)) |>
  mutate(K_c = K - mean(wm$K, na.rm = TRUE))

#### SANITY CHECKS ####
n_subjects <- n_distinct(df$participant)
cat(sprintf("Subjects: %d\n", n_subjects))

trials_per_subject <- df |> count(participant, name = "n_trials")
cat(sprintf("Trials per subject: min = %d, max = %d, mean = %.1f\n",
            min(trials_per_subject$n_trials), max(trials_per_subject$n_trials), mean(trials_per_subject$n_trials)))

cat(sprintf("K range: %.2f to %.2f, missing = %d\n",
            min(wm$K, na.rm = TRUE), max(wm$K, na.rm = TRUE), sum(is.na(wm$K))))

cat(sprintf("Missing reward_oneback/stay in final data: %d\n", sum(is.na(df$reward_oneback) | is.na(df$stay))))

#### WRITE OUTPUT ####
write_csv(df, OUT_CSV)
cat(sprintf("\nSaved -> %s (%d rows, %d subjects)\n", OUT_CSV, nrow(df), n_subjects))
