library(dplyr)
library(tidyverse)
library(purrr)
library(readr)

base_dir <- getwd()
utils_path <- file.path(base_dir, "scripts/helpers_metabolite_utils.R") 
source(utils_path)

# =========================
# Load data
# =========================
data_path <- file.path(base_dir, "original_data/fromSteveMoore_metabolite_list.csv")
moore_list <- read.csv(data_path, header = TRUE)

moore_list.clean <- moore_list %>%
  select(-dataset_id, -study_name, -platform, -merge_variable,
         -input_chemical_id, -input_comp_id, -input_cas, -input_metabid)
moore_list.processing <- moore_list.clean 


# =========================
# Clean names + HMDBS
# =========================
moore_list.processing.clean.1 <- standardize_hmdb(moore_list.processing, "input_hmdb_id")
rows_w_or_or_slash_name.orig <- moore_list.processing.clean.1 %>%
  filter(
    stringr::str_detect(input_metabolite_name, "\\s+or\\s+(?![^()]*\\))") |
      (
        stringr::str_detect(input_metabolite_name, "/") &
          !stringr::str_detect(input_metabolite_name, "\\([^)]*/[^)]*\\)")
      )
  )

moore_list.processing.clean.2 <- biochemical_name_cleaning(moore_list.processing.clean.1, "input_metabolite_name")
moore_list.processing.clean.3 <- moore_list.processing.clean.2 %>% distinct()

#TODO: need to append this at the end
id_cols.temp <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name")
moore_list.processing.clean <- moore_list.processing.clean.3 %>%
  dplyr::filter(
    !if_all(all_of(id_cols.temp), ~ is.na(.) | . == "")
  )

no_info_rows <-  moore_list.processing.clean.3 %>%
  dplyr::filter(
    if_all(all_of(id_cols.temp), ~ is.na(.) | . == "")
  )

group_cols <- c("input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")

mult_rows <- moore_list.processing.clean %>%
  dplyr::filter(
    Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
df_clean <- moore_list.processing.clean %>%
  dplyr::filter(
    !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
test_df <- collapse_if_consistent(df_clean, group_cols, id_cols)

group_cols <- c("metabolite_id")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")
test_df.2 <- collapse_if_consistent(test_df, group_cols, id_cols)


group_cols <- c("metabolite_id", "input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")
test_df.3 <- collapse_if_consistent(test_df.2, group_cols, id_cols)

test_df.3 <- test_df.3 %>% select(-.consistent)
moore_list.processing.clean <- dplyr::bind_rows(test_df.3, mult_rows)
moore_list.processing.clean <- moore_list.processing.clean %>% mutate("multiple_identifiers" = FALSE, "multiple_isomers" = FALSE, "original_index" = 1:nrow(moore_list.processing.clean), "original"=TRUE)

# =========================
# Parallel (i.e. across all 3 cols) splitting for "#" concatenation
# =========================

columns_to_split <- c(
  "metabolite_id",
  "input_hmdb_id",
  "input_kegg",
  "input_chebi",
  "input_chemspider",
  "input_inchikey"
)

# Rows that need splitting hmdb identifiers
rows_with_punct.hmdb <- moore_list.processing.clean %>%
  filter(stringr::str_detect(input_hmdb_id, "[;_#,|]"))


# Rows that need splitting pubchem identifiers
rows_with_punct.pubchem <- moore_list.processing.clean %>%
  filter(stringr::str_detect(input_pubchem, "[;_#,|]"))

# Rows that need splitting kegg identifiers
# TODO: some kegg identifiers do not even have kegg identifiers (i.e. adenine containing, serine and threonine metabolism)
rows_with_punct.kegg <- moore_list.processing.clean %>%
  filter(stringr::str_detect(input_kegg, "[;_#,|]"))

# Rows that need splitting inchikey identifiers -> NONE
rows_with_punct.inchikey <- moore_list.processing.clean %>%
  filter(stringr::str_detect(input_inchikey, "[;_#,|]"))

#TODO: formatting is a bit off, where only some reported in scientific notation
# want to check if cohorts do similar (MetLinkR doesn't have chemspider)
# need to standardize everything
# Rows that need splitting chemspider identifiers -> NONE
rows_with_punct.chemspider <- moore_list.processing.clean %>%
  filter(stringr::str_detect(input_chemspider, "[;_#,|]"))

# Looking at the rows that have multiple identifiers and performing manual splitting
# These metabolite ids were done via manual expection
# HMDB00072_HMDB00958, hmdb + kegg
# 
metabolite_id_mult_ids.parallel_splitting <- c("HMDB00072_HMDB00958")

#Adding rows
aconitate_mult_ids <- rows_with_punct.kegg %>% filter(metabolite_id == "HMDB00072_HMDB00958")
aconitate_trans <- aconitate_mult_ids %>% mutate(input_hmdb_id = "HMDB00958", input_kegg="C02341", "multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)
aconitate_cis <- aconitate_mult_ids %>% mutate(input_hmdb_id="HMDB00072", input_kegg="C00417", "multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)
aconitate_trans_indv <- aconitate_mult_ids %>% mutate(input_hmdb_id = "HMDB00958", input_kegg="C02341", input_metabolite_name="aconitate_trans", input_metabolite_name_final="aconitate_trans", "multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)
aconitate_cis_indv <- aconitate_mult_ids %>% mutate(input_hmdb_id="HMDB00072", input_kegg="C00417", input_metabolite_name="aconitate_cis", input_metabolite_name_final="aconitate_cis","multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)
aconitate_trans_base <- aconitate_mult_ids %>% mutate(input_hmdb_id = "HMDB00958", input_kegg="C02341", input_metabolite_name="aconitate", input_metabolite_name_final="aconitate","multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)
aconitate_cis_base <- aconitate_mult_ids %>% mutate(input_hmdb_id="HMDB00072", input_kegg="C00417", input_metabolite_name="aconitate", input_metabolite_name_final="aconitate","multiple_identifiers" = TRUE, "multiple_isomers" = TRUE, "original"=FALSE)

moore_list.split.1 <- moore_list.processing.clean %>% bind_rows(aconitate_trans, aconitate_cis, aconitate_trans_indv, aconitate_cis_indv, aconitate_trans_base, aconitate_cis_base)
aconitate_mult_ids.orig_idx <- aconitate_mult_ids$original_index
moore_list.split.1 <- moore_list.split.1 %>%
  mutate(
    multiple_identifiers = if_else(
      metabolite_id == "HMDB00072_HMDB00958" & 
        original_index %in% aconitate_mult_ids.orig_idx,
      TRUE,
      multiple_identifiers
    ),
    multiple_isomers = if_else(
      metabolite_id == "HMDB00072_HMDB00958" & 
        original_index %in% aconitate_mult_ids.orig_idx,
      TRUE,
      multiple_isomers
    ),
    original = if_else(
      metabolite_id == "HMDB00072_HMDB00958" & 
        original_index %in% aconitate_mult_ids.orig_idx,
      FALSE,
      original
    )
  )

# Isolating rows that still need splitting
rows_with_punct.hmdb.1 <- rows_with_punct.hmdb %>% filter(!(original_index==aconitate_mult_ids.orig_idx))
rows_with_punct.pubchem.1 <- rows_with_punct.pubchem %>% filter(!(original_index==aconitate_mult_ids.orig_idx))
rows_with_punct.kegg.1 <- rows_with_punct.kegg %>% filter(!(original_index==aconitate_mult_ids.orig_idx))

rows_wo_punct <- moore_list.processing.clean %>% filter(!(original_index %in% c(rows_with_punct.hmdb$original_index, rows_with_punct.kegg$original_index, rows_with_punct.pubchem$original_index)))

moore_list.processing.clean.1 <- moore_list.split.1 %>%
  dplyr::mutate(
    multiple_identifiers = original_index %in% c(
      rows_with_punct.hmdb$original_index,
      rows_with_punct.kegg$original_index,
      rows_with_punct.pubchem$original_index
    )
  )

#TODO: do we want to keep an orig copy of the rows w mult identifiers as well?
regex_sep <- "\\||;|,|\\#|_"
expanded.rows_with_punct.hmdb <- rows_with_punct.hmdb.1 %>%
  tidyr::separate_rows(all_of("input_hmdb_id"), sep = regex_sep) %>%
  dplyr::mutate(!!"input_hmdb_id" := stringr::str_trim(.data[["input_hmdb_id"]]), !!"multiple_identifiers" := TRUE, !!"original":=FALSE)

expanded.rows_with_punct.pubchem <- rows_with_punct.pubchem.1 %>%
  tidyr::separate_rows(all_of("input_pubchem"), sep = regex_sep) %>%
  dplyr::mutate(!!"input_pubchem" := stringr::str_trim(.data[["input_pubchem"]]), !!"multiple_identifiers" := TRUE, !!"original":=FALSE)

expanded.rows_with_punct.kegg <- rows_with_punct.kegg.1 %>%
  tidyr::separate_rows(all_of("input_kegg"), sep = regex_sep) %>%
  dplyr::mutate(!!"input_kegg" := stringr::str_trim(.data[["input_kegg"]]), !!"multiple_identifiers" := TRUE, !!"original":=FALSE)

expanded.moore_list <- dplyr::bind_rows(moore_list.processing.clean.1, expanded.rows_with_punct.hmdb, expanded.rows_with_punct.pubchem, expanded.rows_with_punct.kegg)

expanded.moore_list <- standardize_hmdb(expanded.moore_list, "input_hmdb_id")
expanded.moore_list <- expanded.moore_list %>% distinct()


expanded.moore_list <- expanded.moore_list %>%
  mutate(
    original_index = as.character(original_index),
    contains_original = original == TRUE,
    contains_derived  = original == FALSE
  )

#collapse
group_cols <- c("input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")

mult_rows <- expanded.moore_list %>%
  dplyr::filter(
    Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
df_clean <- expanded.moore_list %>%
  dplyr::filter(
    !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
expanded.test_df <- collapse_if_consistent_with_provenance(df_clean, group_cols, id_cols)

group_cols <- c("metabolite_id")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")
expanded.test_df.2 <- collapse_if_consistent_with_provenance(expanded.test_df, group_cols, id_cols)


group_cols <- c("metabolite_id", "input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")
expanded.test_df.3 <- collapse_if_consistent_with_provenance(expanded.test_df.2, group_cols, id_cols)

expanded.test_df.3 <- expanded.test_df.3 %>% mutate(input_hmdb_id = as.character(input_hmdb_id))
expanded.test_df.3 <- expanded.test_df.3 %>% mutate(input_pubchem = as.character(input_pubchem))
group_cols <- c("input_hmdb_id", "input_pubchem")
id_cols <- c("input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")
expanded.test_df.4 <- collapse_if_consistent_with_provenance(expanded.test_df.3, group_cols, id_cols)

group_cols <- c("input_pubchem")
id_cols <- c("input_hmdb_id", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")
expanded.test_df.5 <- collapse_if_consistent_with_provenance(expanded.test_df.4, group_cols, id_cols)

group_cols <- c("input_pubchem","input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")
expanded.test_df.6 <- collapse_if_consistent_with_provenance(expanded.test_df.5, group_cols, id_cols)

mult_rows <- mult_rows %>%
  mutate(original_index = as.character(original_index))
expanded.moore_list.1 <- dplyr::bind_rows(expanded.test_df.4, mult_rows)

# =========================
# Handle VAR and B suffixes
# =========================
var_report <- investigate_suffix_variants(expanded.moore_list.1, "metabolite_id", "_VAR\\d+$")
b_safe <- var_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

expanded.moore_list.2 <- expanded.moore_list.1 %>%
  mutate(metabolite_id = ifelse(
    stringr::str_replace(metabolite_id, "_[A-Z]+\\d+$", "") %in% b_safe,
    stringr::str_replace(metabolite_id, "_[A-Z]+\\d+$", ""),
    metabolite_id
  ))

b_report <- investigate_suffix_variants(expanded.moore_list.2, "metabolite_id", "_B\\d+$")
b2_safe <- b_report %>% filter(all_names_identical, all_identifiers_identical) %>% pull(base_id)

expanded.moore_list.3 <- expanded.moore_list.2 %>%
  mutate(metabolite_id = ifelse(
    stringr::str_replace(metabolite_id, "_[A-Z]+\\d+$", "") %in% b2_safe,
    stringr::str_replace(metabolite_id, "_[A-Z]+\\d+$", ""),
    metabolite_id
  ))

var_unsafe <- var_report %>% filter(!all_names_identical, !all_identifiers_identical) %>% pull(base_id)
b_unsafe <- b_report %>% filter(!all_names_identical, !all_identifiers_identical) %>% pull(base_id)

#Collapse once again
expanded.moore_list.3 <- expanded.moore_list.3 %>% distinct()

group_cols <- c("metabolite_id")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")

mult_rows.1 <- expanded.moore_list.3 %>%
  dplyr::filter(
    Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
df_clean.1 <- expanded.moore_list.3 %>%
  dplyr::filter(
    !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
expanded.bvar.test_df <- collapse_if_consistent_with_provenance(df_clean.1, group_cols, id_cols)

group_cols <- c("metabolite_id", "input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi")
expanded.bvar.test_df.2 <- collapse_if_consistent_with_provenance(expanded.bvar.test_df, group_cols, id_cols)

expanded.bvar.test_df.2 <- expanded.bvar.test_df.2 %>% mutate(input_hmdb_id = as.character(input_hmdb_id))
expanded.bvar.test_df.2 <- expanded.bvar.test_df.2 %>% mutate(input_pubchem = as.character(input_pubchem))

group_cols <- c("input_hmdb_id", "input_pubchem")
id_cols <- c("input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")
expanded.bvar.test_df.3 <- collapse_if_consistent_with_provenance(expanded.bvar.test_df.2, group_cols, id_cols)

mult_rows.1 <- mult_rows.1 %>%
  mutate(original_index = as.character(original_index))
expanded.moore_list.4 <- dplyr::bind_rows(expanded.bvar.test_df.3, mult_rows.1)

# =========================================================
# Investigating "or", "/" concatenation in biochemical names
# =========================================================
rows_w_or_name <- expanded.moore_list.4 %>% filter(stringr::str_detect(input_metabolite_name_final, "\\s+or\\s+(?![^()]*\\))"))
rows_w_slash_name <- expanded.moore_list.4 %>%
  filter(stringr::str_detect(input_metabolite_name, "/") &
           !stringr::str_detect(input_metabolite_name, "\\([^)]*/[^)]*\\)"))
rows_w_underscore_name <-  expanded.moore_list.4 %>% filter(stringr::str_detect(input_metabolite_name, "_"))


rows_w_or_or_slash_name <- expanded.moore_list.4 %>%
  filter(
    stringr::str_detect(input_metabolite_name_final, "\\s+or\\s+(?![^()]*\\))") |
      (
        stringr::str_detect(input_metabolite_name, "/") &
          !stringr::str_detect(input_metabolite_name, "\\([^)]*/[^)]*\\)")
      )
  )

rows_no_or_or_slash_name <- expanded.moore_list.4 %>%
  filter(
    !(
      stringr::str_detect(input_metabolite_name_final, "\\s+or\\s+(?![^()]*\\))") |
        (
          stringr::str_detect(input_metabolite_name, "/") &
            !stringr::str_detect(input_metabolite_name, "\\([^)]*/[^)]*\\)")
        )
    )
  )

# rows_w_or_or_slash_name was given input to Claude for manual splitting, output is being read below
rows_split_path <- file.path(base_dir, "original_data/rows_split.csv")
splitted_rows_w_or_or_slash_name.draft <- read.csv(rows_split_path)
splitted_rows_w_or_or_slash_name <- splitted_rows_w_or_or_slash_name.draft %>% group_by(Unnamed..0) %>%
  mutate(original_index = max(original_index), original=FALSE, contains_original=FALSE, contains_derived=TRUE, n_source_rows=1) %>%
  select(-Unnamed..0) %>%
  ungroup()

# rows_w_or_or_slash_name was given input to Claude for manual splitting, output is being read below
splitted_rows_w_or_or_slash_name.draft.1 <- read.csv("biochemical_names_split.csv")
splitted_rows_w_or_or_slash_name <- splitted_rows_w_or_or_slash_name.draft.1 %>%
  mutate(
    # fix boolean columns first
    across(c(contains_original, contains_derived, original, multiple_identifiers, multiple_isomers), ~ case_when(
      . %in% c("True", "TRUE", "true") ~ TRUE,
      . %in% c("False", "FALSE", "false") ~ FALSE,
      is.na(.) ~ NA,
      TRUE ~ as.logical(.)
    )),
    
    original_index_num = as.numeric(original_index)
  ) %>%
  group_by(Unnamed..0) %>%
  mutate(
    is_main = original_index_num == max(original_index_num, na.rm = TRUE),
    
    # preserve main row, fix split rows
    original = if_else(is_main, original, FALSE),
    contains_original = if_else(is_main, contains_original, FALSE),
    contains_derived  = if_else(is_main, contains_derived, TRUE),
    
    # reset counts only for split rows
    n_source_rows = if_else(is_main, n_source_rows, 1)
  ) %>%
  ungroup() %>%
  select(-original_index_num, -is_main)

splitted_rows_w_or_or_slash_name.1 <- splitted_rows_w_or_or_slash_name %>% group_by(Unnamed..0) %>%
  mutate(original_index = max(original_index)) %>% ungroup()

splitted_rows_w_or_or_slash_name.1 <- splitted_rows_w_or_or_slash_name.1 %>% select(-Unnamed..0)

expanded.moore_list.5 <- rbind(rows_no_or_or_slash_name, splitted_rows_w_or_or_slash_name.1)

# Collapse
expanded.moore_list.5 <- expanded.moore_list.5 %>% distinct()

group_cols <- c("metabolite_id")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "input_metabolite_name_final")

mult_rows.2 <- expanded.moore_list.5 %>%
  dplyr::filter(
    Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
df_clean.2 <- expanded.moore_list.5 %>%
  dplyr::filter(
    !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
expanded.biochem.test_df <- collapse_if_consistent_with_provenance(df_clean.2, group_cols, id_cols)

group_cols <- c("input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "metabolite_id")
expanded.biochem.test_df.2 <- collapse_if_consistent_with_provenance(expanded.biochem.test_df, group_cols, id_cols)

group_cols <- c("input_hmdb_id")
id_cols <- c("input_metabolite_name_final", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "metabolite_id")
expanded.biochem.test_df.3 <- collapse_if_consistent_with_provenance(expanded.biochem.test_df.2, group_cols, id_cols)

expanded.moore_list.6 <- bind_rows(expanded.biochem.test_df.3, mult_rows.2)


rows_with_descriptive_terms <-  expanded.moore_list.6 %>% filter(detect_descriptive_terms(input_metabolite_name_final))

expanded.moore_list.7 <- expanded.moore_list.6 %>%
  mutate(
    input_metabolite_name_final = clean_metabolite_names(input_metabolite_name_final)
  )

group_cols <- c("input_metabolite_name_final")
id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "metabolite_id")

mult_rows.3 <- expanded.moore_list.7 %>%
  dplyr::filter(
    Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
df_clean.3 <- expanded.moore_list.7 %>%
  dplyr::filter(
    !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
  )
expanded.biochemclean.test_df <- collapse_if_consistent_with_provenance(df_clean.3, group_cols, id_cols)



# =========================
# Final output
# =========================
moore_list.final <- expanded.biochemclean.test_df %>% distinct()
output_path <- file.path(base_dir, "output/moore_list_preprocessed_revised.csv")
write.csv(moore_list.final, output_path, row.names = FALSE)
