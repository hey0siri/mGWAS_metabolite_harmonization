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
  regex_sep <- "\\||;|,|\\#"
  
  input_df %>%
    mutate(
      !!sym(input_hmdb_col) := ifelse(
        !is.na(.data[[input_hmdb_col]]) &
          .data[[input_hmdb_col]] != "" &
          nchar(.data[[input_hmdb_col]]) != 11 &
          !stringr::str_detect(.data[[input_hmdb_col]], regex_sep),
        .data[[input_hmdb_col]] %>% stringr::str_replace("^HDMB", "HMDB") %>% stringr::str_replace("HMDB", "HMDB00"),
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

# expand_dataframes <- function(
#     input_df,
#     input_column,
#     punctuation = c("|", ";", ",", "#"),
#     final_suffix = "_final"
# ) {
#   browser()
#   regex_sep <- stringr::str_c(sapply(punctuation, escape_regex), collapse = "|")
#   final_column <- paste0(input_column, final_suffix)
#   
#   input_df %>%
#     mutate(
#       !!input_column := as.character(.data[[input_column]]),
#       !!final_column := .data[[input_column]]
#     ) %>%
#     tidyr::separate_rows(all_of(final_column), sep = regex_sep) %>%
#     mutate(!!final_column := stringr::str_trim(.data[[final_column]]))
# }

expand_dataframes <- function(
    input_df,
    input_column,
    punctuation = c("|", ";", ",", "#")
) {
  browser()
  regex_sep <- stringr::str_c(sapply(punctuation, escape_regex), collapse = "|")
  
  input_df <- input_df %>%
    dplyr::mutate(
      original_index = dplyr::row_number(),
      !!input_column := as.character(.data[[input_column]])
    )
  
  expanded <- input_df %>%
    tidyr::separate_rows(all_of(input_column), sep = regex_sep) %>%
    dplyr::mutate(!!input_column := stringr::str_trim(.data[[input_column]]))
  
  # Combine original + expanded
  dplyr::bind_rows(input_df, expanded)
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

expand_parallel_columns <- function(df, cols, delim = "#") {
  
  # Convert all character
  df <- df %>%
    mutate(across(all_of(cols), as.character))
  
  # Verify equal entries across all columns
  entry_len_check <- df %>%
    filter(stringr::str_detect(.data[[cols[1]]], delim)) %>%
    mutate(
      across(
        all_of(cols),
        ~ stringr::str_count(.x, fixed(delim)) + 1,
        .names = "n_{.col}"
      )
    )
  
  unequal_rows <- entry_len_check %>%
    filter(
      dplyr::n_distinct(c_across(starts_with("n_"))) > 1
    )
  
  if (nrow(unequal_rows) > 0) {
    stop("ERROR: Misaligned delimiter counts detected across columns. Inspect rows before splitting.")
  }
  
  # Parallel split
  df %>%
    tidyr::separate_longer_delim(
      cols = all_of(cols),
      delim = delim
    ) %>%
    mutate(across(all_of(cols), stringr::str_trim))
}

expand_mixed_delimiters <- function(df, cols, delim = "#") {
  
  expanded_df <- df %>%
    rowwise() %>%
    mutate(
      split_list = list({
  
        # extract current row
        row_vals <- cur_data()
        
        # split each column by delimiter
        splits <- map(cols, ~ {
          val <- row_vals[[.x]]
          str_split(val, fixed(delim))[[1]] %>% str_trim()
        })
        names(splits) <- cols
        
        lengths_vec <- map_int(splits, length)
        
        # Case 1: All equal length → parallel split
        if(length(unique(lengths_vec)) == 1){
          tibble::as_tibble(splits)
          
          # Case 2: Only one column has multiple values → cross with singletons
        } else if(sum(lengths_vec > 1) == 1){
          expand.grid(splits, stringsAsFactors = FALSE) %>% tibble::as_tibble()
          
          # Case 3: lengths mismatch → flag
        } else {
          stop(paste("Mismatched lengths for row:", cur_row()))
        }
        
      })
    ) %>%
    # Remove original columns BEFORE unnest
    select(-all_of(cols)) %>%
    unnest(split_list) %>%
    ungroup()
    
    expanded_df <- expanded_df %>%
      mutate(
        original = FALSE,
        contains_original = FALSE,
        contains_derived = TRUE
      )
    original_df <- df %>%
      mutate(
        contains_derived = TRUE
      )
    
    bind_rows(original_df, expanded_df)
    
}


# =========================
# Collapse rows with identical information
# =========================

has_delimiter <- function(x) {
  grepl("\\||;|,|\\#", x)
}

collapse_without_splitting <- function(df, grouping_cols, hmdb_col) {
  browser()
  
  regex_sep <- "\\||;|,|\\#"
  
  # flag rows with multiple IDs
  df <- df %>%
    dplyr::mutate(
      has_multiple_ids = stringr::str_detect(.data[[hmdb_col]], regex_sep)
    )
  
  multi_rows <- df %>%
    dplyr::filter(has_multiple_ids)
  
  single_rows <- df %>%
    dplyr::filter(!has_multiple_ids)
  
  # Collapse by take first non-NA value per column
  collapsed <- single_rows %>%
    dplyr::group_by(dplyr::across(all_of(grouping_cols))) %>%
    dplyr::summarise(
      dplyr::across(
        everything(),
        ~ {
          vals <- .[!is.na(.) & . != ""]
          if (length(vals) == 0) NA else vals[1]
        }
      ),
      .groups = "drop"
    )
  
  dplyr::bind_rows(collapsed, multi_rows) %>%
    dplyr::select(-has_multiple_ids)
}

is_consistent <- function(x) {
  vals <- unique(x[!is.na(x) & x != ""])
  length(vals) <= 1
}

#IGNORE
collapse_if_consistent <- function(df, group_cols, id_cols) {
  collapse_fn <- function(x) {
    vals <- x[!is.na(x) & x != ""]
    if (length(vals) == 0) return(NA)
    vals[1]
  }
  
  df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    
    dplyr::mutate(
      .consistent = purrr::pmap_lgl(
        dplyr::across(all_of(id_cols)),
        ~ TRUE  # placeholder, we compute below
      )
    ) %>%
    
    dplyr::group_modify(~ {
      data <- .x
      
      consistent_all <- all(sapply(data[id_cols], is_consistent))
      
      if (consistent_all) {
        # collapse
        data %>%
          dplyr::summarise(
            dplyr::across(everything(), collapse_fn),
            .groups = "drop"
          )
      } else {
        data
      }
    }) %>%
    dplyr::ungroup()
}

#IGNORE
collapse_if_consistent_orig_idx <- function(df, group_cols, id_cols) {
  df <- df %>%
    mutate(original_index = as.character(original_index))
  
  collapse_fn <- function(x) {
    vals <- x[!is.na(x) & x != ""]
    if (length(vals) == 0) return(NA)
    vals[1]
  }
  
  collapse_index <- function(x) {
    x %>%
      as.character() %>%
      stringr::str_split("\\|") %>%  # handle already concatenated
      unlist() %>%
      unique() %>%
      paste(collapse = "|")
  }
  
  df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    
    dplyr::group_modify(~ {
      data <- .x
      
      consistent_all <- all(sapply(data[id_cols], is_consistent))
      
      if (isTRUE(consistent_all)) {
        data %>%
          dplyr::summarise(
            
            # special column
            original_index = collapse_index(original_index),
            
            # everything else
            dplyr::across(
              -original_index,
              collapse_fn
            ),
            
            .groups = "drop"
          )
        
      } else {
        # keep rows as-is
        data
      }
    }) %>%
    dplyr::ungroup()
}

collapse_if_consistent_with_provenance <- function(df, group_cols, id_cols) {
  
  collapse_fn <- function(x) {
    vals <- x[!is.na(x) & x != ""]
    if (length(vals) == 0) return(NA)
    vals[1]
  }
  
  collapse_index <- function(x) {
    x %>%
      as.character() %>%
      stringr::str_split("\\|") %>%
      unlist() %>%
      unique()
  }
  
  df %>%
    mutate(original_index = as.character(original_index)) %>%
    group_by(across(all_of(group_cols))) %>%
    
    group_modify(~ {
      data <- .x
      
      orig_vals <- if ("contains_original" %in% names(data)) {
        data$contains_original
      } else {
        data$original == TRUE
      }
      
      deriv_vals <- if ("contains_derived" %in% names(data)) {
        data$contains_derived
      } else {
        data$original == FALSE
      }
      
      consistent_all <- all(sapply(data[id_cols], is_consistent))
      
      if (consistent_all) {
        idx_vec <- collapse_index(data$original_index)
        
        data %>%
          summarise(
            # provenance
            original_index = paste(idx_vec, collapse = "|"),
            contains_original = any(orig_vals, na.rm = TRUE),
            contains_derived  = any(deriv_vals, na.rm = TRUE),
            n_source_rows = length(idx_vec),
            
            # rest
            across(
              -c(original_index, original, contains_original, contains_derived),
              collapse_fn
            ),
            
            .groups = "drop"
          )
        
      } else {
        data
      }
    }) %>%
    ungroup()
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
  
  return (input_df %>%
    mutate(
      !!final_column := .data[[input_column]] %>%
        as.character() %>%
        stringr::str_trim() %>%
        stringr::str_replace_all("\\*", "") %>%
        stringr::str_replace_all("\\s*\\(\\d+\\)$", "") %>%
        stringr::str_replace_all(regex_sep, "_") %>%
        stringr::str_replace_all("_+", "_") %>%
        stringr::str_replace_all("^_|_$", "") %>%
        stringr::str_replace_all("\\s*_\\s*", "_") %>%
        stringr::str_to_lower()
    ))
}

clean_metabolite_names <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace("^_|_$", "") %>%
    
    stringr::str_remove_all(
      stringr::str_c(
        "(",
        paste(c(
          "plasma","serum","urine","blood",
          "mmol_l","umol_l","nmol_l","mol_l","mmol","umol","nmol",
          "signal_area","signal","area",
          "g_l","mg_dl", "_in_*_*dl",
          "s\\d+"
        ), collapse="|"),
        ")"
      )
    ) %>%
    
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace("^_|_$", "")
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
  terms <- c("\\bplasma\\b", "\\bserum\\b", "\\burine\\b", "\\bblood\\b", "\\b[nmu]?mol\\b", ".s5")
  pattern <- stringr::str_c(terms, collapse = "|")
  stringr::str_detect(tolower(x), pattern)
}

detect_multiple_ids <- function(x) {
  stringr::str_count(toupper(x), "(HMDB|PUBCHEM|CHEM|CSID|ME)[0-9]+") >= 2
}
