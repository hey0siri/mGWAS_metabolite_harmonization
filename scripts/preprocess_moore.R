library(dplyr)
library(tidyverse)
library(purrr)

source("helpers_metabolite_utils.R")

# =========================
# Load data
# =========================
moore_list <- read.csv("fromSteveMoore_metabolite_list.csv")
load(file.path("compileduids.RData"))

moore_list.clean <- moore_list %>%
  select(-dataset_id, -study_name, -platform, -merge_variable,
         -input_chemical_id, -input_comp_id, -input_cas, -input_metabid)

# =========================
# Identify concatenated identifiers
# =========================
columns_to_split <- c(
  "metabolite_id",
  "input_hmdb_id",
  "input_kegg",
  "input_chebi",
  "input_chemspider",
  "input_inchikey"
)

for (col in columns_to_split) {
  if (col == "input_hmdb_id") {
    investigate_concatenated_rows(moore_list.clean, col, punctuation = c("|", ";", ",", "#", "_"))
  } else {
    investigate_concatenated_rows(moore_list.clean, col)
  }
}

investigate_concatenated_rows(moore_list.clean, "input_pubchem")

# =========================
# PubChem splitting logic
# =========================
rows_with_pubchem <- readr::read_csv("rows_with_delimiters_input_pubchem.csv")
rows_remaining <- readr::read_csv("rows_with_delimiters_input_pubchem_remaining.csv")

rows_with_pubchem <- rows_with_pubchem %>%
  mutate(decision = map(input_pubchem, decide_pubchem_split)) %>%
  unnest(decision)

rows_to_split <- rows_with_pubchem %>% filter(split_input_pubchem)
rows_to_keep <- rows_with_pubchem %>% filter(!split_input_pubchem) %>%
  mutate(input_pubchem_final = input_pubchem)

rows_to_split <- expand_dataframes(rows_to_split, "input_pubchem")

rows_to_split <- rows_to_split %>%
  mutate(input_pubchem = as.character(input_pubchem))

rows_to_keep <- rows_to_keep %>%
  mutate(input_pubchem = as.character(input_pubchem))

rows_remaining <- rows_remaining %>%
  mutate(input_pubchem = as.character(input_pubchem))

moore_pubchem_expanded <- bind_rows(rows_to_split, rows_to_keep, rows_remaining)

# =========================
# Sequential splitting of other identifiers
# =========================
moore_list.split_dfs <- vector("list", length(columns_to_split) + 1)
moore_list.split_dfs[[1]] <- moore_pubchem_expanded

for (i in seq_along(columns_to_split)) {
  col <- columns_to_split[[i]]
  if (col == "input_hmdb_id") {
    moore_list.split_dfs[[i + 1]] <- expand_dataframes(
      moore_list.split_dfs[[i]],
      col,
      punctuation = c("|", ";", ",", "#", "_")
    )
  } else {
    moore_list.split_dfs[[i + 1]] <- expand_dataframes(moore_list.split_dfs[[i]], col)
  }
}

moore_expanded <- moore_list.split_dfs[[length(moore_list.split_dfs)]]

# =========================
# Biochemical name normalization
# =========================
moore_expanded <- biochemical_name_cleaning(moore_expanded, "input_metabolite_name")

# =========================
# Handle VAR and B suffixes
# =========================
var_report <- investigate_suffix_variants(moore_expanded, "metabolite_id_final", "_VAR\\d+$")
b_safe <- var_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

moore_expanded <- moore_expanded %>%
  mutate(metabolite_id_final_1 = ifelse(
    stringr::str_replace(metabolite_id_final, "_[A-Z]+\\d+$", "") %in% b_safe,
    stringr::str_replace(metabolite_id_final, "_[A-Z]+\\d+$", ""),
    metabolite_id_final
  ))

b_report <- investigate_suffix_variants(moore_expanded, "metabolite_id_final_1", "_B\\d+$")
b2_safe <- b_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

moore_expanded <- moore_expanded %>%
  mutate(metabolite_id_final_2 = ifelse(
    stringr::str_replace(metabolite_id_final_1, "_[A-Z]+\\d+$", "") %in% b2_safe,
    stringr::str_replace(metabolite_id_final_1, "_[A-Z]+\\d+$", ""),
    metabolite_id_final_1
  ))

var_unsafe <- var_report %>% filter(!all_names_identical, !all_identifiers_identical) %>% pull(base_id)
b_unsafe <- b_report %>% filter(!all_names_identical, !all_identifiers_identical) %>% pull(base_id)

write.csv(var_unsafe, "moore_metab_id_var_variants.csv")
write.csv(b_unsafe, "moore_metab_id_b_variants.csv")

# =========================
# Descriptive term and multi-ID flags
# =========================
moore_expanded <- moore_expanded %>%
  mutate(
    has_descriptive_terms = detect_descriptive_terms(input_metabolite_name) |
      detect_descriptive_terms(metabolite_id_final_2),
    has_multiple_ids = detect_multiple_ids(metabolite_id_final_2)
  )

moore_list.descriptive_terms <- moore_expanded %>% filter(has_descriptive_terms)
moore_list.multiple_ids <- moore_expanded %>% filter(has_multiple_ids)

write.csv(moore_list.descriptive_terms, "moore_list_descriptive_terms.csv")
write.csv(moore_list.multiple_ids, "moore_list_multiple_ids.csv")

# =========================
# Manual metabolite ID cleaning based on inspection
# =========================
# manual clean: PUBCHEM52922056

# =========================
# HMDB standardization
# =========================
moore_expanded <- standardize_hmdb(moore_expanded, "input_hmdb_id_final")
moore_expanded <- standardize_hmdb(moore_expanded, "metabolite_id_final_2")

# =========================
# Final output
# =========================
moore_list.final <- moore_expanded %>% distinct()

write.csv(moore_list.final, "moore_list_preprocessed.csv", row.names = FALSE)
