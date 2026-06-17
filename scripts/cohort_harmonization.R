library(dplyr)
library(tidyverse)
library(readxl)

base_dir <- getwd()
utils_path <- file.path(base_dir, "scripts/helpers_metabolite_utils.R") 
source(utils_path)

# =========================
# Load metabolite databases
# =========================
metabolite_database <- read.csv(file.path(base_dir, "output/draft_comets_moore_merged.csv"))

# =================================================
# Load cohort metabolites and column names
# =================================================
# Change this as needed
sample_cohort_file_name <- file.path(base_dir, "sample_cohorts/fromEGG_alspac_metabolite_list.xlsx")
sample_cohort_name <- "ALSPAC"
sample_cohort <- readxl::read_xlsx(file.path(base_dir, "sample_cohorts/fromEGG_alspac_metabolite_list.xlsx"))

#TODO: do not hard code this
sample_cohort_biochem_name <- "CHEMICAL_NAME"
sample_cohort_hmdb_name <- "HMDB"
sample_cohort_kegg_name <- "KEGG"
sample_cohort_pubchem_name <- "PUBCHEM"
sample_cohort_inchikey_name <- "INCHIKEY"
sample_cohort_chemspider_name <- "CHEMSPIDER"

# ========================================================
# Expand rows and standardize HMDB + biochemical names
# ========================================================
metabolite_database.1 <- biochemical_name_cleaning(metabolite_database, "MLR_Harmonized_Name_moore")
metabolite_database.1 <- biochemical_name_cleaning(metabolite_database.1, "MLR_Harmonized_Name_comets")
sample_cohort.clean.1 <- sample_cohort %>% filter(TYPE == "NAMED")
sample_cohort.clean.1 <- standardize_hmdb(sample_cohort.clean.1, sample_cohort_hmdb_name)


sample_cohort.clean.2 <- biochemical_name_cleaning(sample_cohort.clean.1, sample_cohort_biochem_name)
sample_cohort.clean.3 <- sample_cohort.clean.2 %>% distinct()
sample_cohort_biochem_name.orig <- sample_cohort_biochem_name
sample_cohort_biochem_name <- "CHEMICAL_NAME_final"

# Assigning original indices to keep track of splits
sample_cohort.clean.3 <- sample_cohort.clean.3 %>% mutate(original_index = 1:nrow(sample_cohort.clean.3))

regex_sep <- "\\||;|,|\\#|_"
# Rows that need splitting hmdb identifiers
rows_with_punct.hmdb <- sample_cohort.clean.3 %>%
  filter(stringr::str_detect(.data[[sample_cohort_hmdb_name]], regex_sep))


# Rows that need splitting pubchem identifiers
rows_with_punct.pubchem <- sample_cohort.clean.3 %>%
  filter(stringr::str_detect(.data[[sample_cohort_pubchem_name]], "[;_#,|]"))

# Rows that need splitting kegg identifiers
# TODO: some kegg identifiers do not even have kegg identifiers (i.e. adenine containing, serine and threonine metabolism)
rows_with_punct.kegg <- sample_cohort.clean.3 %>%
  filter(stringr::str_detect(.data[[sample_cohort_kegg_name]], "[;_#,|]"))

# Rows that need splitting inchikey identifiers -> NONE
rows_with_punct.inchikey <- sample_cohort.clean.3 %>%
  filter(stringr::str_detect(.data[[sample_cohort_inchikey_name]], "[;_#,|]"))


rows_with_all_identifiers <- sample_cohort.clean.3 %>% filter(original_index %in% rows_with_punct.hmdb$original_index 
                                                              & original_index %in% rows_with_punct.pubchem$original_index
                                                              & original_index %in% rows_with_punct.kegg$original_index)

rows_with_hmdb_kegg <- sample_cohort.clean.3 %>% filter(original_index %in% rows_with_punct.hmdb$original_index 
                                                              & original_index %in% rows_with_punct.kegg$original_index
                                                              & !(original_index %in% rows_with_all_identifiers$original_index))
rows_with_hmdb_pubchem <- sample_cohort.clean.3 %>% filter(original_index %in% rows_with_punct.hmdb$original_index 
                                                        & original_index %in% rows_with_punct.pubchem$original_index
                                                        & !(original_index %in% rows_with_all_identifiers$original_index))
rows_with_kegg_pubchem <- sample_cohort.clean.3 %>% filter(original_index %in% rows_with_punct.kegg$original_index 
                                                           & original_index %in% rows_with_punct.pubchem$original_index
                                                           & !(original_index %in% rows_with_all_identifiers$original_index))

rows_with_punct.hmdb.1 <- rows_with_punct.hmdb %>% filter(!(original_index %in% rows_with_all_identifiers$original_index) &
                                                          !(original_index %in% rows_with_hmdb_kegg$original_index) &
                                                          !(original_index %in% rows_with_hmdb_pubchem$original_index))
rows_with_punct.kegg.1 <- rows_with_punct.kegg %>% filter(!(original_index %in% rows_with_all_identifiers$original_index) &
                                                            !(original_index %in% rows_with_hmdb_kegg$original_index) &
                                                            !(original_index %in% rows_with_kegg_pubchem$original_index))
rows_with_punct.pubchem.1 <- rows_with_punct.pubchem %>% filter(!(original_index %in% rows_with_all_identifiers$original_index) &
                                                            !(original_index %in% rows_with_hmdb_pubchem$original_index) &
                                                            !(original_index %in% rows_with_kegg_pubchem$original_index))
rows_without_concatenation <- sample_cohort.clean.3 %>% filter(!(original_index %in%rows_with_punct.kegg$original_index) &
                                                                 !(original_index %in%rows_with_punct.hmdb$original_index) &
                                                                 !(original_index %in%rows_with_punct.pubchem$original_index))

# Creating separate long-formatted tables instead of parallel splitting
hmdb_lookup <- sample_cohort.clean.3 %>% tidyr::separate_rows(all_of(sample_cohort_hmdb_name), sep = regex_sep)
pubchem_lookup <- sample_cohort.clean.3 %>% tidyr::separate_rows(all_of(sample_cohort_pubchem_name), sep = regex_sep)
kegg_lookup <- sample_cohort.clean.3 %>% tidyr::separate_rows(all_of(sample_cohort_kegg_name), sep = regex_sep)


# ==============================================================
# Matchup with HMDB (btw metabolite database and sample cohort)
# ==============================================================
matching_hmdbs.metab_db <- metabolite_database %>% filter(input_hmdb_id %in% hmdb_lookup[[sample_cohort_hmdb_name]] | hmdb_id %in% hmdb_lookup[[sample_cohort_hmdb_name]])
matching_hmdbs.sample_cohort <- hmdb_lookup %>% filter(!is.na(.data[[sample_cohort_hmdb_name]]) & .data[[sample_cohort_hmdb_name]] != "NA" & 
                                                                   (.data[[sample_cohort_hmdb_name]] %in% metabolite_database$input_hmdb_id | 
                                                                      .data[[sample_cohort_hmdb_name]] %in% metabolite_database$hmdb_id))
hmdb_matches <- metabolite_database %>%
  mutate(match_hmdb = input_hmdb_id) %>%
  bind_rows(
    metabolite_database %>%
      mutate(match_hmdb = hmdb_id)
  ) %>%
  filter(!is.na(match_hmdb))

hmdb_matches <- hmdb_matches %>% select(
  match_hmdb,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_hmdbs.sample_cohort.1 <- matching_hmdbs.sample_cohort %>%
  left_join(
    hmdb_matches,
    by = setNames("match_hmdb", sample_cohort_hmdb_name),
    relationship = "many-to-many"
  )

matching_hmdbs.sample_cohort.1 <- matching_hmdbs.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_hmdbs.sample_cohort.1 <- matching_hmdbs.sample_cohort.1 %>% distinct()

matching_hmdbs.sample_cohort.mult_identifiers <- matching_hmdbs.sample_cohort.1 %>%
  group_by(HMDB) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_hmdbs.sample_cohort.2 <- matching_hmdbs.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      HMDB
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_hmdbs.sample_cohort.3 <- matching_hmdbs.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_hmdbs.sample_cohort.3 <- matching_hmdbs.sample_cohort.3 %>% distinct()

remaining_matching <- sample_cohort.clean.3 %>% filter(!(original_index %in% matching_hmdbs.sample_cohort.3$original_index))

# ==============================================================
# Matchup with PubChem (btw metabolite database and sample cohort)
# ==============================================================
remaining_matching.pubchem_expanded <- remaining_matching %>% tidyr::separate_rows(all_of(sample_cohort_pubchem_name), sep = regex_sep)
matching_pubchem.sample_cohort <- remaining_matching.pubchem_expanded %>% filter(!is.na(.data[[sample_cohort_pubchem_name]]) & .data[[sample_cohort_pubchem_name]] != "NA" & 
                                                                   (.data[[sample_cohort_pubchem_name]] %in% metabolite_database$input_pubchem))
pubchem_matches <- metabolite_database %>%
  mutate(match_pubchem = input_pubchem) %>%
  filter(!is.na(match_pubchem))

pubchem_matches <- pubchem_matches %>% select(
  match_pubchem,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_pubchem.sample_cohort.1 <- matching_pubchem.sample_cohort %>%
  left_join(
    pubchem_matches,
    by = setNames("match_pubchem", sample_cohort_pubchem_name),
    relationship = "many-to-many"
  )

matching_pubchem.sample_cohort.1 <- matching_pubchem.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_pubchem.sample_cohort.1 <- matching_pubchem.sample_cohort.1 %>% distinct()

matching_pubchem.sample_cohort.mult_identifiers <- matching_pubchem.sample_cohort.1 %>%
  group_by(.data[[sample_cohort_pubchem_name]]) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_pubchem.sample_cohort.2 <- matching_pubchem.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      .data[[sample_cohort_pubchem_name]]
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_pubchem.sample_cohort.3 <- matching_pubchem.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_pubchem.sample_cohort.3 <- matching_pubchem.sample_cohort.3 %>% distinct()

remaining_matching.1 <- remaining_matching %>% filter(!(original_index %in% matching_pubchem.sample_cohort.3$original_index))

# ==============================================================
# Matchup with Kegg (btw metabolite database and sample cohort)
# ==============================================================
remaining_matching.kegg_expanded <- remaining_matching.1 %>% tidyr::separate_rows(all_of(sample_cohort_kegg_name), sep = regex_sep)
matching_kegg.sample_cohort <- remaining_matching.kegg_expanded %>% filter(!is.na(.data[[sample_cohort_kegg_name]]) & .data[[sample_cohort_kegg_name]] != "NA" & 
                                                                  (.data[[sample_cohort_kegg_name]] %in% metabolite_database$input_kegg))
kegg_matches <- metabolite_database %>%
  mutate(match_kegg = input_kegg) %>%
  filter(!is.na(match_kegg))

kegg_matches <- kegg_matches %>% select(
  match_kegg,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_kegg.sample_cohort.1 <- matching_kegg.sample_cohort %>%
  left_join(
    kegg_matches,
    by = setNames("match_kegg", sample_cohort_kegg_name),
    relationship = "many-to-many"
  )

matching_kegg.sample_cohort.1 <- matching_kegg.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_kegg.sample_cohort.1 <- matching_kegg.sample_cohort.1 %>% distinct()

matching_kegg.sample_cohort.mult_identifiers <- matching_kegg.sample_cohort.1 %>%
  group_by(.data[[sample_cohort_kegg_name]]) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_kegg.sample_cohort.2 <- matching_kegg.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      .data[[sample_cohort_kegg_name]]
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_kegg.sample_cohort.3 <- matching_kegg.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_kegg.sample_cohort.3 <- matching_kegg.sample_cohort.3 %>% distinct()

remaining_matching.2 <- remaining_matching.1 %>% filter(!(original_index %in% matching_kegg.sample_cohort.3$original_index))

# ==============================================================
# Matchup with chemspider (btw metabolite database and sample cohort)
# ==============================================================
remaining_matching.chemspider_expanded <- remaining_matching.2 %>% tidyr::separate_rows(all_of(sample_cohort_chemspider_name), sep = regex_sep)
matching_chemspider.sample_cohort <- remaining_matching.chemspider_expanded %>% filter(!is.na(.data[[sample_cohort_chemspider_name]]) & .data[[sample_cohort_chemspider_name]] != "NA" & 
                                                                             (.data[[sample_cohort_chemspider_name]] %in% metabolite_database$input_chemspider))
chemspider_matches <- metabolite_database %>%
  mutate(
    match_chemspider = as.character(input_chemspider)
  ) %>%
  filter(!is.na(match_chemspider))

chemspider_matches <- chemspider_matches %>% select(
  match_chemspider,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_chemspider.sample_cohort.1 <- matching_chemspider.sample_cohort %>%
  left_join(
    chemspider_matches,
    by = setNames("match_chemspider", sample_cohort_chemspider_name),
    relationship = "many-to-many"
  )

matching_chemspider.sample_cohort.1 <- matching_chemspider.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_chemspider.sample_cohort.1 <- matching_chemspider.sample_cohort.1 %>% distinct()

matching_chemspider.sample_cohort.mult_identifiers <- matching_chemspider.sample_cohort.1 %>%
  group_by(.data[[sample_cohort_chemspider_name]]) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_chemspider.sample_cohort.2 <- matching_chemspider.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      .data[[sample_cohort_chemspider_name]]
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_chemspider.sample_cohort.3 <- matching_chemspider.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_chemspider.sample_cohort.3 <- matching_chemspider.sample_cohort.3 %>% distinct()

remaining_matching.3 <- remaining_matching.2 %>% filter(!(original_index %in% matching_chemspider.sample_cohort.3$original_index))

# ==============================================================
# Matchup with inchikey (btw metabolite database and sample cohort)
# ==============================================================
matching_inchikey.sample_cohort <- remaining_matching.3 %>% filter(!is.na(.data[[sample_cohort_inchikey_name]]) & .data[[sample_cohort_inchikey_name]] != "NA" & 
                                                                                         (.data[[sample_cohort_inchikey_name]] %in% metabolite_database$input_inchikey))
inchikey_matches <- metabolite_database %>%
  mutate(
    match_inchikey = as.character(input_inchikey)
  ) %>%
  filter(!is.na(match_inchikey))

inchikey_matches <- inchikey_matches %>% select(
  match_inchikey,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_inchikey.sample_cohort.1 <- matching_inchikey.sample_cohort %>%
  left_join(
    inchikey_matches,
    by = setNames("match_inchikey", sample_cohort_inchikey_name),
    relationship = "many-to-many"
  )

matching_inchikey.sample_cohort.1 <- matching_inchikey.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_inchikey.sample_cohort.1 <- matching_inchikey.sample_cohort.1 %>% distinct()

mamatching_inchikey.sample_cohort.mult_identifiers <- matching_inchikey.sample_cohort.1 %>%
  group_by(.data[[sample_cohort_inchikey_name]]) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_inchikey.sample_cohort.2 <- matching_inchikey.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      .data[[sample_cohort_inchikey_name]]
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_inchikey.sample_cohort.3 <- matching_inchikey.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_inchikey.sample_cohort.3 <- matching_inchikey.sample_cohort.3 %>% distinct()

remaining_matching.4 <- remaining_matching.3 %>% filter(!(original_index %in% matching_inchikey.sample_cohort.3$original_index))

# ========================================================================
# Matchup with Biochemical Name (btw metabolite database and sample cohort)
# ========================================================================
matching_name.sample_cohort <- remaining_matching.4 %>% filter(!is.na(.data[[sample_cohort_biochem_name]]) & .data[[sample_cohort_biochem_name]] != "NA" & .data[[sample_cohort_biochem_name]] != "")

name_matches <- metabolite_database.1 %>%
  mutate(match_name = input_metabolite_name_final) %>%
  bind_rows(
    metabolite_database %>%
      mutate(match_name = biochemical_final)
  ) %>%
  bind_rows(
    metabolite_database %>%
      mutate(match_name = MLR_Harmonized_Name_moore)
  ) %>%
  bind_rows(
    metabolite_database %>%
      mutate(match_name = MLR_Harmonized_Name_comets)
  ) %>%
  filter(!is.na(match_name))

name_matches <- name_matches %>% mutate(match_name=as.character(match_name))

name_matches <- name_matches %>% select(
  match_name,
  uid_01,
  metabolite_id,
  MLR_Harmonized_Name_moore,
  MLR_Harmonized_Name_comets
)

matching_name.sample_cohort.1 <- matching_name.sample_cohort %>%
  left_join(
    name_matches,
    by = setNames("match_name", sample_cohort_biochem_name),
    relationship = "many-to-many"
  )

matching_name.sample_cohort.1 <- matching_name.sample_cohort.1 %>%
  rename(moore_metabolite_id = metabolite_id)
matching_name.sample_cohort.1 <- matching_name.sample_cohort.1 %>% distinct()

matching_name.sample_cohort.mult_identifiers <- matching_name.sample_cohort.1 %>%
  group_by(.data[[sample_cohort_biochem_name]]) %>%
  summarise(
    n_uid = n_distinct(uid_01, na.rm = TRUE),
    n_moore = n_distinct(moore_metabolite_id, na.rm = TRUE),
    n_mlr_moore = n_distinct(MLR_Harmonized_Name_moore, na.rm = TRUE),
    n_mlr_comets = n_distinct(MLR_Harmonized_Name_comets, na.rm = TRUE)
  ) %>%
  filter(n_uid > 1 | n_moore > 1 | n_mlr_moore > 1 | n_mlr_comets > 1)

matching_name.sample_cohort.2 <- matching_name.sample_cohort.1 %>%
  mutate(
    final_metabolite_id = coalesce(
      as.character(uid_01),
      moore_metabolite_id,
      MLR_Harmonized_Name_moore,
      MLR_Harmonized_Name_comets,
      .data[[sample_cohort_biochem_name]]
    ),
    metabolite_id_origin = case_when(
      !is.na(uid_01) ~ "COMETS",
      !is.na(moore_metabolite_id) ~ "Moore",
      !is.na(MLR_Harmonized_Name_moore) |
        !is.na(MLR_Harmonized_Name_comets) ~ "MetLinkR",
      TRUE ~ "sample_cohort"
    )
  )

matching_name.sample_cohort.3 <- matching_name.sample_cohort.2 %>% select(-uid_01, -moore_metabolite_id, -MLR_Harmonized_Name_moore, -MLR_Harmonized_Name_comets)
matching_name.sample_cohort.3 <- matching_name.sample_cohort.3 %>% distinct()

remaining_matching.5 <- remaining_matching.4 %>% filter(!(original_index %in% matching_name.sample_cohort.3$original_index))

# ========================================================================
# Bind all outputs
# ========================================================================
remaining_matching.5 <- remaining_matching.5 %>% mutate(
  final_metabolite_id = coalesce(.data[[sample_cohort_hmdb_name]],
                                 .data[[sample_cohort_pubchem_name]],
                                 .data[[sample_cohort_kegg_name]],
                                 .data[[sample_cohort_biochem_name]]),
  metabolite_id_origin = "sample_cohort"
)

final_harmonization <- rbind(remaining_matching.5, matching_hmdbs.sample_cohort.3, matching_name.sample_cohort.3,
                             matching_pubchem.sample_cohort.3, matching_kegg.sample_cohort.3,
                             matching_chemspider.sample_cohort.3, matching_inchikey.sample_cohort.3)
final_harmonization_w_metabid <- final_harmonization %>% filter(
  metabolite_id_origin != "sample_cohort"
)
final_harmonization_wo_metabid <- final_harmonization %>% filter(
  metabolite_id_origin == "sample_cohort"
)

write.csv(final_harmonization, file.path(base_dir, "cohort_harmonization/output/alspac_sample_harmonization_wo_mlr.csv"))
# ========================================================================
# Install MetLinkR
# ========================================================================
# Locally install RaMP
if (!requireNamespace("devtools", quietly = TRUE)) {
  
  install.packages("devtools", repos = "https://cloud.r-project.org")
  
}


library(devtools)
#load_all("/Users/heysiri/MetLinkR")
#install_github("ncats/RAMP-DB")
devtools::install_github("hey0siri/MetLinkR")

# Load the package
library(RaMP)

# initializes the RaMP database object, downloading and caching the latest SQLite database
# if no version already exists in local cache.
rampDB <- RaMP()

# note that you can use the following method to check database versions hosted in your
# computer's local cache and databases that are available to download in our remote repository.
RaMP::listAvailableRaMPDbVersions()

# using that list of available RaMP DB versions, one can specify the database version to use
# if the selected version is not available on your computer, but is in our remote repository at GitHub,
# the SQLite DB file will be automatically downloaded into local file cache.
# RaMP is using the BiocFileCache package to manage a local file cache.
rampDB <- RaMP(version = "2.5.4")


# =============================================
# Build input dataframe + MetLinkR 
# =============================================
remaining_cohort_metabolite_file_path <- file.path(base_dir, "/output/remaining_cohort_metabolites.csv")
write.csv(final_harmonization_wo_metabid, remaining_cohort_metabolite_file_path)
metlinkr_input_file <- data.frame("FileNames" = c(remaining_cohort_metabolite_file_path),
                                  "ShortFileName" = c(sample_cohort_name),
                                  "HMDB" = c(sample_cohort_hmdb_name),
                                  "Metabolite_Name" = c(sample_cohort_biochem_name),
                                  "PubChem_CID" = c(sample_cohort_pubchem_name),
                                  "KEGG" = c(sample_cohort_kegg_name),
                                  "LIPIDMAPS" = c(NA),
                                  "chebi" = c(NA))

write.csv(metlinkr_input_file, file.path(base_dir, "metlinkr_input_files/metlinkr_input_file_for_merging_cohort.csv"))

# To separate the MetLinkR output that is generated from metabolite database generation
setwd(file.path(base_dir, "/cohort_harmonization"))
metLinkR::harmonizeInputSheets(file.path(base_dir, "metlinkr_input_files/metlinkr_input_file_for_merging_cohort.csv"))
setwd(file.path(base_dir))

# =====================================================
# Grabbing MetLinkR Standardized Names as metabolite IDs
# =====================================================
cohort_metlinkr_output_file <- readxl::read_xlsx(file.path(base_dir, "cohort_harmonization/metLinkR_output/mapping_library.xlsx"))

id_cols <- c(
  sample_cohort_hmdb_name,
  sample_cohort_biochem_name,
  sample_cohort_pubchem_name,
  sample_cohort_kegg_name,
  sample_cohort_chemspider_name,
  sample_cohort_inchikey_name
)

match_value_col = paste0("Input name (", sample_cohort_name, ")")

final_harmonization_wo_metabid <- final_harmonization_wo_metabid %>%
  mutate(across(all_of(id_cols), as.character))

metlinkr_lookup.cohort <- cohort_metlinkr_output_file %>%
  select(
    Harmonized_Name = `Harmonized name`, 
    match_value = .data[[match_value_col]]) %>%
  filter(match_value != "-") %>% 
  separate_rows(
    match_value,
    sep = "\\s*;\\s*"
  ) %>%
  mutate(
    match_value = stringr::str_trim(match_value)
  ) %>%
  distinct()


cohort_metlinkr.hmdb_merged <- final_harmonization_wo_metabid %>%
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_hmdb_name)
  ) %>%
  rename(HMDB_Harmonized_Name = Harmonized_Name)

cohort_metlinkr.pubchem_merged <- cohort_metlinkr.hmdb_merged %>% 
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_pubchem_name)
  ) %>%
  rename(Pubchem_Harmonized_Name = Harmonized_Name)

cohort_metlinkr.chemspider_merged <- cohort_metlinkr.pubchem_merged %>% 
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_chemspider_name)
  ) %>%
  rename(Chemspider_Harmonized_Name = Harmonized_Name)

cohort_metlinkr.kegg_merged <- cohort_metlinkr.chemspider_merged %>% 
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_kegg_name)
  ) %>%
  rename(Kegg_Harmonized_Name = Harmonized_Name)

cohort_metlinkr.inchikey_merged <- cohort_metlinkr.kegg_merged %>% 
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_inchikey_name)
  ) %>%
  rename(Inchikey_Harmonized_Name = Harmonized_Name)



cohort_metlinkr.biochemical_merged <- cohort_metlinkr.inchikey_merged %>%
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_biochem_name.orig)
  ) %>%
  rename(Biochemical_final_Harmonized_Name = Harmonized_Name) %>%
  left_join(
    metlinkr_lookup.cohort,
    by = setNames("match_value", sample_cohort_biochem_name)
  ) %>%
  rename(Biochemical_orig_Harmonized_Name = Harmonized_Name) %>%
  mutate(
    Biochemical_Harmonized_Name = coalesce(
      Biochemical_final_Harmonized_Name,
      Biochemical_orig_Harmonized_Name
    )
  ) %>%
  
  select(
    -Biochemical_final_Harmonized_Name,
    -Biochemical_orig_Harmonized_Name
  )

cohort_metlinkr.metlinkr <- cohort_metlinkr.biochemical_merged %>%
  mutate(
    MLR_Harmonized_Name = coalesce(
      HMDB_Harmonized_Name,
      Biochemical_Harmonized_Name,
      Inchikey_Harmonized_Name,
      Kegg_Harmonized_Name,
      Chemspider_Harmonized_Name,
      Pubchem_Harmonized_Name
    )
  ) %>%
  select(-HMDB_Harmonized_Name, -Biochemical_Harmonized_Name, -Inchikey_Harmonized_Name, -Kegg_Harmonized_Name, -Chemspider_Harmonized_Name, -Pubchem_Harmonized_Name) 

mlr_lookup <- cohort_metlinkr.metlinkr %>%
  filter(!is.na(MLR_Harmonized_Name)) %>%
  group_by(original_index) %>%
  summarise(
    MLR_Harmonized_Name = first(MLR_Harmonized_Name),
    .groups = "drop"
  )
cohort_metlinkr.metlinkr.1 <- cohort_metlinkr.metlinkr %>%
  select(-MLR_Harmonized_Name) %>%
  left_join(mlr_lookup, by = "original_index")

cohort_metlinkr.metlinkr.2 <- cohort_metlinkr.metlinkr.1 %>%
  mutate(
    final_metabolite_id = if_else(
      !is.na(MLR_Harmonized_Name),
      MLR_Harmonized_Name,
      final_metabolite_id
    ),
    metabolite_id_origin = if_else(
      !is.na(MLR_Harmonized_Name),
      "MetLinkR",
      metabolite_id_origin
    )
  )
cohort_metlinkr.metlinkr.2 <- cohort_metlinkr.metlinkr.2 %>% select(-"MLR_Harmonized_Name")
final_harmonization_w_mlr <- rbind(final_harmonization_w_metabid, cohort_metlinkr.metlinkr.2)
write.csv(final_harmonization_w_mlr, file.path(base_dir, "cohort_harmonization/output/alspac_sample_harmonization_w_mlr.csv"))

