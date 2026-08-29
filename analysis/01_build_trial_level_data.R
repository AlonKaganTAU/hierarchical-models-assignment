rm(list = ls())

#### SETUP ####
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

DATA_TASK <- "data/raw/task"
DATA_WM   <- "data/raw/wm"
OUT_CSV   <- "data/trial_level.csv"

# Source: illusion-of-control repo, analysis/wm_wsls/wm_wsls.R (fourth pilot,
# fixed inter-trial interval, N = 35). This script reuses that pipeline's WM
# capacity computation and trial-pairing logic, but keeps all 35 subjects (no
# understood/felt exclusion) and returns trial-level data rather than
# per-subject WSLS betas, since the assignment models the trial-level outcome
# directly with a reward_oneback x K cross-level interaction.

#### WORKING MEMORY CAPACITY (COWAN'S K) ####
# One K per set size, K = set_size * (hit rate + correct rejection rate - 1),
# where condition "d" = change trials (hits) and "s" = no-change trials
# (correct rejections). A participant's K is the mean over the two set sizes.
wm_files <- list.files(DATA_WM, pattern = "^ioc-wm_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

wm <- read_csv(wm_files, id = "file", show_col_types = FALSE,
               col_types = cols(.default = col_character())) |>
  mutate(participant = str_extract(file, "[a-f0-9]{24}")) |>
  filter(trial_name == "test_squares", block_type == "exp") |>
  mutate(set_size = as.numeric(set_size),
         correct  = tolower(accuracy) == "true") |>
  summarise(accuracy = mean(correct), .by = c(participant, set_size, condition)) |>
  pivot_wider(names_from = condition, values_from = accuracy) |>
  mutate(k_set_size = set_size * (d + s - 1)) |>
  summarise(K = mean(k_set_size), .by = participant)

#### TRIAL-LEVEL STAY / REWARD PAIRS ####
task_files <- list.files(DATA_TASK, pattern = "^ioc-all_[a-f0-9]{24}_SESSION.*\\.csv$", full.names = TRUE)

df_raw <- read_csv(task_files, id = "file", show_col_types = FALSE,
                   col_types = cols(.default = col_character())) |>
  mutate(participant = str_extract(file, "[a-f0-9]{24}"))

# The lag is taken before any filtering, so choice_oneback and reward_oneback
# always refer to the trial that physically preceded the current one within the
# same block. Trials are kept only when both they and their predecessor carry a
# valid (non-timed-out) response.
df <- df_raw |>
  filter(task == "gambling_choice", block_number != "training") |>
  mutate(block_number    = as.numeric(block_number),
         trial_number    = as.numeric(trial_number),
         reward          = as.numeric(reward),
         is_choice_valid = as.logical(is_choice_valid)) |>
  arrange(participant, block_number, trial_number) |>
  mutate(choice_oneback = lag(choice_key),
         valid_oneback  = lag(is_choice_valid),
         reward_oneback = lag(reward),
         .by = c(participant, block_number)) |>
  filter(is_choice_valid, valid_oneback, !is.na(reward_oneback)) |>
  mutate(stay = as.integer(choice_key == choice_oneback)) |>
  select(participant, block_number, trial_number, reward_oneback, stay)

#### ADD LEVEL-2 PREDICTOR AND CENTER EACH VARIABLE AT ITS OWN LEVEL ####
# K is time-invariant (level 2) and is centered on the mean of the analysed
# sample, so that 0 = average capacity in this study. reward_oneback is a
# measured time-varying (level 1) predictor, so it is additionally split into a
# person mean (between-person, level 2) and a person-mean-centered deviation
# (within-person, level 1) to keep the two sources of variance from being
# smushed into a single coefficient.
subjects <- df |>
  distinct(participant) |>
  inner_join(wm, by = "participant") |>
  mutate(K_c = K - mean(K))

df <- df |>
  inner_join(subjects, by = "participant") |>
  mutate(reward_oneback_pm = mean(reward_oneback),
         reward_oneback_wp = reward_oneback - reward_oneback_pm,
         .by = participant)

reward_pm_grand <- df |>
  distinct(participant, reward_oneback_pm) |>
  pull(reward_oneback_pm) |>
  mean()

df <- df |>
  mutate(reward_oneback_pm_c = reward_oneback_pm - reward_pm_grand)

#### SANITY CHECKS ####
cat(sprintf("Subjects: %d, trials: %d\n", nrow(subjects), nrow(df)))

trials_per_subject <- count(df, participant, name = "n_trials")
cat(sprintf("Trials per subject: min = %d, max = %d, mean = %.1f\n",
            min(trials_per_subject$n_trials), max(trials_per_subject$n_trials),
            mean(trials_per_subject$n_trials)))

cat(sprintf("Blocks per subject: %s\n",
            paste(unique(summarise(df, b = n_distinct(block_number), .by = participant)$b), collapse = ", ")))

cat(sprintf("K: range %.2f to %.2f, mean %.2f, between-subject SD %.2f\n",
            min(subjects$K), max(subjects$K), mean(subjects$K), sd(subjects$K)))

cat(sprintf("K_c mean over subjects (should be 0): %.12f\n", mean(subjects$K_c)))

cat(sprintf("Person reward rate: range %.2f to %.2f, between-subject SD %.3f\n",
            min(df$reward_oneback_pm), max(df$reward_oneback_pm),
            sd(distinct(df, participant, reward_oneback_pm)$reward_oneback_pm)))

cat(sprintf("Missing values in modelled columns: %d\n",
            sum(is.na(df$reward_oneback) | is.na(df$stay) | is.na(df$K_c))))

#### WRITE OUTPUT ####
write_csv(df, OUT_CSV)
cat(sprintf("\nSaved -> %s (%d rows)\n", OUT_CSV, nrow(df)))
