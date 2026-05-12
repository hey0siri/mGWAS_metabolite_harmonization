library(dplyr)
library(tidyverse)
library(purrr)

base_dir <- getwd()
utils_path <- file.path(base_dir, "scripts/helpers_metabolite_utils.R") 
source(utils_path)

#TODO: figure out the intermediate files/the biochemical names
# =========================
# Load data
# =========================
data_path <- file.path(base_dir, "original_data/compileduids.RData")
load(data_path)

comets.clean <- mastermetid %>%
  select(-uidsource, -main_class, -chemical_id, -comp_id, -metid) %>%
  mutate(
    original_index = row_number(),
    original = TRUE,
    contains_original = TRUE,
    contains_derived = FALSE,
    n_source_rows = 1,
    multiple_identifiers = FALSE,
    multiple_isomers = FALSE
  )

# # =========================
# # Identify concatenated identifiers
# # =========================
# columns_to_split <- c(
#   "uid_01", "hmdb_id", "biochemical"
# )
# 
# for (col in columns_to_split){
#   investigate_concatenated_rows(comets.clean, col, punctuation = c("#"))
# }
# 
# investigate_concatenated_rows(comets.clean, "biochemical", punctuation = c(" or "))
# investigate_concatenated_rows(comets.clean, "uid_01", punctuation = c("_"))

# =========================
# Parallel (i.e. across all 3 cols) splitting for "#" concatenation
# =========================

columns_to_split <- c("uid_01", "hmdb_id", "biochemical")

# Rows that need splitting
rows_with_hash <- comets.clean %>%
  filter(stringr::str_detect(uid_01, fixed("#")))

# Rows that do not
rows_without_hash <- comets.clean %>%
  filter(!stringr::str_detect(uid_01, fixed("#")))

# Perform aligned split
rows_split <- expand_mixed_delimiters(
  rows_with_hash,
  cols = columns_to_split,
  delim = "#"
)
# Check for # delim in remaining rows in hmdb IDs
# No concatenated hmdbs with #, so skipping split
rows_without_hash.1 <- rows_without_hash %>%
  filter(!stringr::str_detect(hmdb_id, fixed("#")))

# Check for _ delim in remaining rows in hmdb IDs
# No concatenated hmdbs with _, so skipping split
rows_without_hash.2 <- rows_without_hash %>%
  filter(!stringr::str_detect(hmdb_id, fixed("_")))

# Recombine
comets.expanded_hash <- bind_rows(rows_without_hash, rows_split)

# =========================
# Handle VAR and B suffixes
# =========================
# Skipping over, no var variants in uids for comets
var_report <- investigate_suffix_variants(
  comets.expanded_hash,
  id_column = "uid_01",
  suffix_regex = "_VAR\\d+$",
  metabolite_col_name = "biochemical",
  compare_columns = c("hmdb_id", "biochemical")
)

# Variants with _BX after uid
b_report <- investigate_suffix_variants(
  comets.expanded_hash,
  id_column = "uid_01",
  suffix_regex = "_B\\d+$",
  metabolite_col_name = "biochemical",
  compare_columns = c("hmdb_id", "biochemical")
)

b_safe <- b_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

comets.b_processed <- comets.expanded_hash %>%
  mutate(uid_01 = ifelse(
    stringr::str_replace(uid_01, "_[A-Z]+\\d+$", "") %in% b_safe,
    stringr::str_replace(uid_01, "_[A-Z]+\\d+$", ""),
    uid_01
  )) %>% unique()

# =========================
# Parallel (i.e. across all 3 cols) splitting for "_" concatenation
# =========================
#TODO: ask abotu what to do for certain cases
# Check for _ delim in remaining rows in uids
rows_with_underscore <- comets.b_processed %>%
  filter(stringr::str_detect(uid_01, fixed("_")))

#once the logic for splitting is solidified above, will need to split here and then do or concat splitting

# =========================
# Splitting for " or " concatenation
# =========================
hash_original_indices <- rows_with_hash$original_index

rows_with_or <- comets.b_processed %>%
  filter(stringr::str_detect(biochemical, "\\s+or\\s+(?![^()]*\\))"),
         (!(original_index %in% hash_original_indices & contains_original)))

expand_or_outside_parentheses <- function(df, column) {
  
  or_outside_parens <- "\\s+or\\s+(?![^()]*\\))"
  
  df %>%
    mutate({{ column }} := as.character({{ column }})) %>%
    tidyr::separate_rows(
      {{ column }},
      sep = or_outside_parens
    ) %>%
    mutate({{ column }} := stringr::str_trim({{ column }}))
}

rows_with_or <- comets.b_processed %>%
  filter(stringr::str_detect(biochemical, "\\s+or\\s+(?![^()]*\\))"),
         (!(original_index %in% hash_original_indices & contains_original)))

rows_without_or <- comets.b_processed %>%
  filter(!stringr::str_detect(biochemical, "\\s+or\\s+(?![^()]*\\))") | 
           (original_index %in% hash_original_indices & contains_original))

expanded_or_rows <- expand_or_outside_parentheses(
  rows_with_or,
  biochemical
) 

rows_with_or <- rows_with_or %>%
  mutate(
    contains_derived = TRUE
  )

comets.expanded_or <- bind_rows(
  rows_without_or,
  rows_with_or,
  expanded_or_rows
)

# comets.expanded_or <- expand_or_outside_parentheses(
#   comets.b_processed,
#   biochemical
# )


# =========================
# Biochemical name normalization
# =========================
comets.expanded.1 <- biochemical_name_cleaning(comets.expanded_or, "biochemical")


# =========================
# HMDB standardization
# =========================
comets.hmdb_cleaned.1 <- standardize_hmdb(comets.expanded.1, "uid_01")
comets.hmdb_cleaned.2 <- standardize_hmdb(comets.hmdb_cleaned.1, "hmdb_id")

# =========================
# Final output
# =========================
comets.final <- comets.hmdb_cleaned.2 %>% distinct()
output_path <- file.path(base_dir, "output/comets_preprocessed_revised.csv")
write.csv(comets.final, output_path, row.names = FALSE)

