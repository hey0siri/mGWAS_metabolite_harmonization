library(dplyr)
library(tidyverse)
library(httr)
library(jsonlite)
library(purrr)

# =========================
# Regex helpers
# =========================
escape_regex <- function(x) {
  if (x %in% c(".", "|", "(", ")", "[", "]", "{", "}", "^", "$", "*", "+", "#", "?")) {
    paste0("\\", x)
  } else {
    x
  }
}

# =========================
# HMDB standardization
# =========================
standardize_hmdb <- function(input_df, input_hmdb_col) {
  input_df %>%
    mutate(
      !!sym(input_hmdb_col) := ifelse(
        !is.na(.data[[input_hmdb_col]]) &
          .data[[input_hmdb_col]] != "" &
          nchar(.data[[input_hmdb_col]]) != 11,
        stringr::str_replace(.data[[input_hmdb_col]], "HMDB", "HMDB00"),
        .data[[input_hmdb_col]]
      )
    )
}

# =========================
# Concatenated identifier inspection
# =========================
investigate_concatenated_rows <- function(
    input_df,
    input_column,
    punctuation = c("|", ";", ",", "#")
) {
  regex_sep <- stringr::str_c(sapply(punctuation, escape_regex), collapse = "|")
  
  rows_with_delimiters <- input_df %>%
    mutate(!!input_column := as.character(.data[[input_column]])) %>%
    filter(grepl(regex_sep, .data[[input_column]]))
  
  rows_remaining <- input_df %>%
    mutate(!!input_column := as.character(.data[[input_column]])) %>%
    filter(!grepl(regex_sep, .data[[input_column]]))
  
  write.csv(rows_with_delimiters,
            paste0("rows_with_delimiters_", input_column, "_", regex_sep, ".csv"),
            row.names = FALSE)
  
  write.csv(rows_remaining,
            paste0("rows_with_delimiters_", input_column, "_", regex_sep, "_remaining.csv"),
            row.names = FALSE)
}

expand_dataframes <- function(
    input_df,
    input_column,
    punctuation = c("|", ";", ",", "#"),
    final_suffix = "_final"
) {
  regex_sep <- stringr::str_c(sapply(punctuation, escape_regex), collapse = "|")
  final_column <- paste0(input_column, final_suffix)
  
  input_df %>%
    mutate(
      !!input_column := as.character(.data[[input_column]]),
      !!final_column := .data[[input_column]]
    ) %>%
    tidyr::separate_rows(all_of(final_column), sep = regex_sep) %>%
    mutate(!!final_column := stringr::str_trim(.data[[final_column]]))
}

split_check.hmdb <- function(input_df, input_column) {
  hmdb_vals <- input_df %>%
    filter(grepl("HMDB", .data[[input_column]])) %>%
    pull(input_column)
  
  message(
    "Any rows with multiple HMDB IDs: ",
    any(stringr::str_count(hmdb_vals, "HMDB") > 1)
  )
}

# =========================
# PubChem helpers
# =========================
query_pubchem_cid <- function(cid) {
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    cid,
    "/property/InChIKey,IUPACName,CanonicalSMILES/JSON"
  )
  
  res <- tryCatch(GET(url), error = function(e) NULL)
  
  if (is.null(res) || httr::status_code(res) != 200) {
    return(tibble(
      pubchem_id = as.character(cid),
      inchikey = NA_character_,
      iupac_name = NA_character_,
      smiles = NA_character_
    ))
  }
  
  parsed <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
  props <- parsed$PropertyTable$Properties
  
  tibble(
    pubchem_id = as.character(cid),
    inchikey = props[["InChIKey"]],
    iupac_name = props[["IUPACName"]],
    smiles = props[["CanonicalSMILES"]]
  )
}

decide_pubchem_split <- function(input_pubchem) {
  pubchem_entries <- stringr::str_split(input_pubchem, ";")[[1]]
  
  inchikeys <- purrr::map_chr(pubchem_entries, ~ query_pubchem_cid(.x)$inchikey)
  
  if (any(is.na(inchikeys))) {
    return(tibble(split_input_pubchem = FALSE,
                  split_reason = "PUBCHEM_LOOKUP_FAILED"))
  }
  
  if (dplyr::n_distinct(inchikeys) == 1) {
    return(tibble(split_input_pubchem = TRUE,
                  split_reason = "IDENTICAL_INCHIKEY"))
  }
  
  cores <- purrr::map_chr(inchikeys, ~ stringr::str_split(.x, "-")[[1]][1])
  
  if (dplyr::n_distinct(cores) == 1) {
    return(tibble(split_input_pubchem = TRUE,
                  split_reason = "SAME_CORE_DIFFERENT_STEREO_OR_SALT"))
  }
  
  tibble(split_input_pubchem = FALSE,
         split_reason = "DIFFERENT_CONNECTIVITY")
}

# =========================
# Biochemical name cleaning
# =========================
biochemical_name_cleaning <- function(
    input_df,
    input_column,
    punctuation = c("-", "/", "¥", ".", ",", ":", ";", "(", ")", "|", "[", "]", "#"),
    final_suffix = "_final"
) {
  regex_sep <- stringr::str_c(sapply(punctuation, escape_regex), collapse = "|")
  final_column <- paste0(input_column, final_suffix)
  
  input_df %>%
    mutate(
      !!final_column := .data[[input_column]] %>%
        as.character() %>%
        stringr::str_trim() %>%
        stringr::str_replace_all(regex_sep, "_") %>%
        stringr::str_replace_all("\\*", "") %>%
        stringr::str_replace_all("_+", "_") %>%
        stringr::str_replace_all("^_|_$", "") %>%
        stringr::str_replace_all("\\s*_\\s*", "_") %>%
        stringr::str_to_lower()
    )
}

# =========================
# Suffix inspection
# =========================
strip_suffix <- function(x, suffix_regex) {
  stringr::str_replace(x, suffix_regex, "")
}

investigate_suffix_variants <- function(
    input_df,
    id_column = "metabolite_id_final",
    suffix_regex = "_[A-Z]+\\d+$",
    metabolite_col_name = "input_metabolite_name_final",
    compare_columns = c(
      "input_metabolite_name_final",
      "input_hmdb_id",
      "input_pubchem",
      "input_kegg",
      "input_chebi",
      "input_chemspider"
    )
) {
  input_df %>%
    mutate(
      base_id = stringr::str_replace(.data[[id_column]], suffix_regex, ""),
      has_suffix = .data[[id_column]] != base_id
    ) %>%
    filter(has_suffix) %>%
    group_by(base_id) %>%
    summarise(
      n_variants = n(),
      variant_ids = paste(unique(.data[[id_column]]), collapse = "; "),
      all_names_identical = dplyr::n_distinct(metabolite_col_name, na.rm = TRUE) <= 1,
      all_identifiers_identical = all(
        dplyr::across(all_of(compare_columns), ~ dplyr::n_distinct(.x, na.rm = TRUE) <= 1)
      ),
      .groups = "drop"
    )
}

# =========================
# Descriptor / ID detection
# =========================
detect_descriptive_terms <- function(x) {
  terms <- c("\\bplasma\\b", "\\bserum\\b", "\\burine\\b", "\\bblood\\b", "\\b[nmu]?mol\\b")
  pattern <- stringr::str_c(terms, collapse = "|")
  stringr::str_detect(tolower(x), pattern)
}

detect_multiple_ids <- function(x) {
  stringr::str_count(toupper(x), "(HMDB|PUBCHEM|CHEM|CSID|ME)[0-9]+") >= 2
}
