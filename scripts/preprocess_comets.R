library(dplyr)
library(tidyverse)
library(purrr)

source("helpers_metabolite_utils.R")

# =========================
# Load data
# =========================
load(file.path("compileduids.RData"))

comets.clean <- mastermetid %>%
  select(-uidsource, -main_class, -chemical_id, -comp_id, -metid)

# =========================
# Identify concatenated identifiers
# =========================
columns_to_split <- c(
  "uid_01", "hmdb_id", "biochemical"
)

for (col in columns_to_split){
  investigate_concatenated_rows(comets.clean, col, punctuation = c("#"))
}

investigate_concatenated_rows(comets.clean, "biochemical", punctuation = c(" or "))

# =========================
# Sequential splitting for "#" concatenation
# =========================
comets.split_dfs <- vector("list", length(columns_to_split) + 1)
comets.split_dfs[[1]] <- comets.clean

for (i in seq_along(columns_to_split)) {
  col <- columns_to_split[[i]]
  comets.split_dfs[[i + 1]] <- expand_dataframes(comets.split_dfs[[i]], col, punctuation = c("#"))
}

# =========================
# Splitting for " or " concatenation
# =========================
expand_or_outside_parentheses <- function(
    input_df,
    input_column,
    final_suffix = "_final"
) {
  
  final_column <- paste0(input_column, final_suffix)
  
  # "or" only when NOT inside parentheses
  or_outside_parens <- "\\s+or\\s+(?![^()]*\\))"
  
  input_df %>%
    mutate(
      !!input_column := as.character(.data[[input_column]]),
      !!final_column := .data[[input_column]]
    ) %>%
    tidyr::separate_rows(
      all_of(final_column),
      sep = or_outside_parens
    ) %>%
    mutate(
      !!final_column := stringr::str_trim(.data[[final_column]])
    )
}

comets.expanded.or_split <- expand_or_outside_parentheses(comets.split_dfs[[4]], "biochemical_final")

# =========================
# Biochemical name normalization
# =========================
comets.expanded <- biochemical_name_cleaning(comets.expanded.or_split, "biochemical_final_final")


# TODO: investigate cases where biochemical name root is the same, but differs with (1) or (2)


# =========================
# Handle VAR and B suffixes
# =========================
var_report <- investigate_suffix_variants(comets.expanded, "uid_01_final", "_VAR\\d+$", "biochemical_final_final_final", c("hmdb_id_final","biochemical_final_final_final"))
b_safe <- var_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

comets.expanded.var_variants <- comets.expanded %>%
  mutate(uid_01_final_1 = ifelse(
    stringr::str_replace(uid_01_final, "_[A-Z]+\\d+$", "") %in% b_safe,
    stringr::str_replace(uid_01_final, "_[A-Z]+\\d+$", ""),
    uid_01_final
  ))

b_report <- investigate_suffix_variants(comets.expanded.var_variants, "uid_01_final", "_B\\d+$", "biochemical_final_final_final", c("hmdb_id_final","biochemical_final_final_final"))
b2_safe <- b_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

comets.expanded.b_variants <- comets.expanded.var_variants %>%
  mutate(uid_01_final_2 = ifelse(
    stringr::str_replace(uid_01_final_1, "_[A-Z]+\\d+$", "") %in% b2_safe,
    stringr::str_replace(uid_01_final_1, "_[A-Z]+\\d+$", ""),
    uid_01_final_1
  ))

# filter rows that have multiple hmdb ids concatenated
# TODO: save and investigate
comets.mult_rows_concat <- comets.expanded.b_variants %>% filter(
  grepl("_", uid_01_final_2)
)

# =========================
# HMDB standardization
# =========================
comets.hmdb_cleaned.1 <- standardize_hmdb(comets.expanded.b_variants, "uid_01_final_2")
comets.hmdb_cleaned.2 <- standardize_hmdb(comets.hmdb_cleaned.1, "hmdb_id_final")

# =========================
# Final output
# =========================
comets.final <- comets.hmdb_cleaned.2 %>% distinct()

write.csv(comets.final, "comets_preprocessed.csv", row.names = FALSE)

