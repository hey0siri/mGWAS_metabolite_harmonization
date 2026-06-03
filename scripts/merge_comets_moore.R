library(dplyr)
library(tidyverse)
library(purrr)
base_dir <- getwd()
source(file.path(base_dir, "scripts/helpers_metabolite_utils.R"))
# =========================
# Load preprocessed metabolite dataframes
# =========================
moore.preprocessed.orig <- read.csv(file.path(base_dir, "/output/moore_list_preprocessed_revised.csv"))
comets.preprocessed.orig <- read.csv(file.path(base_dir, "/output/comets_preprocessed_revised.csv"))
moore.metlinkr.combined <- read.csv(file.path(base_dir, "/output/moore_metlinkr_merged.csv"))
comets.metlinkr.combined <- read.csv(file.path(base_dir, "/output/comets_metlinkr_merged.csv"))


# =========================
# Merge based on hmdd -> metabolite id -> biochem name order
# =========================

id_cols <- c(
  "input_hmdb_id",
  "input_pubchem",
  "input_chemspider",
  "input_kegg",
  "input_chebi"
)

moore.metlinkr.combined <- moore.metlinkr.combined %>%
  mutate(across(all_of(id_cols), as.character))


moore.metlinkr.combined <- moore.metlinkr.combined %>%
  mutate(.row_id_moore = row_number())

comets.metlinkr.combined <- comets.metlinkr.combined %>%
  mutate(.row_id_comets = row_number())

moore.metlinkr.combined.hmdb_na <- moore.metlinkr.combined %>% filter(input_hmdb_id == "")
comets.metlinkr.combined.hmdb_na <- comets.metlinkr.combined %>% filter(hmdb_id == "")

moore.metlinkr.combined.hmdb <- moore.metlinkr.combined %>% filter(input_hmdb_id != "")
comets.metlinkr.combined.hmdb <- comets.metlinkr.combined %>% filter(hmdb_id != "")

join_hmdb <- full_join(
  moore.metlinkr.combined.hmdb,
  comets.metlinkr.combined.hmdb,
  by = c("input_hmdb_id" = "hmdb_id"),
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)

join_hmdb <- join_hmdb %>% distinct()
join_hmdb.1 <- join_hmdb %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))

moore_unmatched.hmdb <- anti_join(moore.metlinkr.combined.hmdb,
                                  comets.metlinkr.combined.hmdb,
                                    by = c("input_hmdb_id" = "hmdb_id"))

comets_unmatched.hmdb <- anti_join(comets.metlinkr.combined.hmdb,
                                   moore.metlinkr.combined.hmdb,
                                   by = c("hmdb_id" = "input_hmdb_id"))

moore_unmatched.1 <- rbind(moore.metlinkr.combined.hmdb_na, moore_unmatched.hmdb)
comets_unmatched.1 <- rbind(comets.metlinkr.combined.hmdb_na, comets_unmatched.hmdb)

join_uid <- full_join(
  moore_unmatched.1,
  comets_unmatched.1,
  by = c("metabolite_id" = "uid_01"),
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)

join_uid <- join_uid %>% distinct()
join_uid.1 <- join_uid %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))

moore_unmatched.2 <- anti_join(moore_unmatched.1,
                               comets_unmatched.1,
                               by = c("metabolite_id" = "uid_01"))

comets_unmatched.2 <- anti_join(comets_unmatched.1,
                                   moore_unmatched.1,
                                   by = c("uid_01" = "metabolite_id"))

join_name.1 <- full_join(
  moore_unmatched.2,
  comets_unmatched.2,
  by = c("input_metabolite_name_final" = "biochemical_final"),
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)

join_name.1 <- join_name.1 %>% distinct()
join_name.1.1 <- join_name.1 %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))

moore_unmatched.3 <- anti_join(moore_unmatched.2, comets_unmatched.2, by = c("input_metabolite_name_final" = "biochemical_final"))
comets_unmatched.3 <- anti_join(comets_unmatched.2,
                                moore_unmatched.2,
                                by = c("biochemical_final" = "input_metabolite_name_final"))

moore.metlinkr.combined.mlr_na <- moore_unmatched.3 %>% filter(MLR_Harmonized_Name == "")
comets.metlinkr.combined.mlr_na <- comets_unmatched.3 %>% filter(MLR_Harmonized_Name == "")

moore.metlinkr.combined.mlr <- moore_unmatched.3 %>% filter(MLR_Harmonized_Name != "")
comets.metlinkr.combined.mlr <- comets_unmatched.3 %>% filter(MLR_Harmonized_Name != "")

join.name.2 <- full_join(
  moore.metlinkr.combined.mlr,
  comets.metlinkr.combined.mlr,
  by = "MLR_Harmonized_Name",
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)
join_name.2 <- join.name.2 %>% distinct()
join_name.2.1 <- join_name.2 %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))


moore_unmatched.mlr <- anti_join(moore.metlinkr.combined.mlr, comets.metlinkr.combined.mlr, by="MLR_Harmonized_Name")
comets_unmatched.mlr <- anti_join(comets.metlinkr.combined.mlr, moore.metlinkr.combined.mlr, by="MLR_Harmonized_Name")

moore_unmatched.4 <- rbind(moore_unmatched.mlr, moore.metlinkr.combined.mlr_na)
comets_unmatched.4 <- rbind(comets_unmatched.mlr, comets.metlinkr.combined.mlr_na)

#SKIP
join.name.3 <- full_join(
  moore_unmatched.4,
  comets_unmatched.4,
  by = c("input_metabolite_name" = "biochemical"),
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)
join_name.3 <- join.name.3 %>% distinct()
join_name.3.1 <- join_name.3 %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))

#SKIP
join.uid_to_misc <- full_join(
  moore_unmatched.4,
  comets_unmatched.4,
  by = c("input_chebi" = "uid_01"),
  relationship="many-to-many",
  suffix = c("_moore", "_comets"),
  keep=TRUE
)
join.uid_to_misc <- join.uid_to_misc %>% distinct()
join_name.4.1 <- join.uid_to_misc %>% filter(!is.na(.row_id_moore) & !is.na(.row_id_comets))


all_merged <- rbind(join_hmdb.1, join_uid.1, join_name.1.1,join_name.2.1, join.name.3)

all_merged <- all_merged %>% select(-X_moore, -X_comets)


# =========================
# Collapse duplicate rows
# =========================
# [6] "input_hmdb_id"               "metabolite_id"               "input_pubchem"               "input_metabolite_name"       "input_chemspider"           
# [11] "input_kegg"                  "input_inchikey"              "input_chebi"                 "multiple_identifiers_moore"  "multiple_isomers_moore"     
# [16] "original_moore"              "MLR_Harmonized_Name_moore"   ".row_id_moore"               "uid_01"                      "hmdb_id"                    
# [21] "biochemical"                 "original_index_comets"       "original_comets"             "contains_original_comets"    "contains_derived_comets"    
# [26] "n_source_rows_comets"        "multiple_identifiers_comets" "multiple_isomers_comets"     "biochemical_final"           "MLR_Harmonized_Name_comets" 
# [31] ".row_id_comets"  

# group_cols <- c("input_metabolite_name_final", "biochemical_final")
# id_cols <- c("input_hmdb_id", "input_pubchem", "input_chemspider", "input_kegg", "input_inchikey", "input_chebi", "hmdb_id", "uid_01")
# 
# mult_rows <- all_merged %>%
#   dplyr::filter(
#     Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
#   )
# df_clean <- all_merged %>%
#   dplyr::filter(
#     !Reduce(`|`, lapply(across(all_of(id_cols)), has_delimiter))
#   )
# expanded.test_df <- collapse_if_consistent_with_provenance(df_clean, group_cols, id_cols)
# 
# 
# collapse_unique <- function(x) {
#   
#   vals <- x %>%
#     as.character() %>%
#     na.omit() %>%
#     unique()
#   
#   vals <- vals[vals != ""]
#   
#   if(length(vals) == 0) return(NA)
#   
#   paste(vals, collapse = "|")
# }
# 
# collapse_index <- function(x) {
#   
#   x %>%
#     as.character() %>%
#     str_split("\\|") %>%
#     unlist() %>%
#     unique() %>%
#     sort() %>%
#     paste(collapse="|")
# }
# 
# group_cols <- c(
#   "input_hmdb_id",
#   "hmdb_id",
#   "MLR_Harmonized_Name_moore",
#   "MLR_Harmonized_Name_comets"
# )
# 
# collapsed_merged <- all_merged %>%
#   
#   group_by(across(all_of(group_cols))) %>%
#   
#   summarise(
#     
#     across(
#       everything(),
#       collapse_index
#     ),
#     
#     n_merged_rows = n(),
#     
#     .groups = "drop"
#   )
# 
# write.csv(collapsed_merged, file.path(base_dir, "/output/draft_comets_moore_merged_collapsed.csv"))
write.csv(all_merged, file.path(base_dir, "/output/draft_comets_moore_merged.csv"))