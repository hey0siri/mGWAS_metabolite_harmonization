library(dplyr)
library(tidyverse)
library(purrr)

# =========================
# Load preprocessed metabolite dataframes
# =========================
moore.preprocessed <- read.csv("moore_list_preprocessed.csv")
comets.preprocessed <- read.csv("comets_preprocessed.csv")

# =========================
# Select relevant columns
# =========================
moore.preprocessed <- moore.preprocessed %>% select(
  input_pubchem_final,
  input_hmdb_id_final,
  input_kegg_final,
  input_chebi_final,
  input_chemspider_final,
  input_inchikey_final,
  input_metabolite_name_final,
  metabolite_id_final_2
) %>% distinct()

comets.preprocessed <- comets.preprocessed %>% select(
  hmdb_id_final,
  biochemical_final_final_final,
  uid_01_final_2
) %>% distinct()

# =========================
# Merge based on hmdd -> metabolite id -> biochem name order
# =========================
moore.preprocessed <- moore.preprocessed %>%
  mutate(.row_id_moore = row_number())

comets.preprocessed <- comets.preprocessed %>%
  mutate(.row_id_comets = row_number())

join_hmdb <- full_join(
  moore.preprocessed,
  comets.preprocessed %>% mutate(hmdb_id_final_comets = hmdb_id_final),
  by = c("input_hmdb_id_final" = "hmdb_id_final_comets"),
  suffix = c("_moore", "_comets")
)

join_hmdb <- join_hmdb %>%
  mutate(
    matched_hmdb = !is.na(.row_id_moore) & !is.na(.row_id_comets)
  )

moore_unmatched.hmdb <- anti_join(moore.preprocessed,
                                    comets.preprocessed,
                                    by = c("input_hmdb_id_final" = "hmdb_id_final"))

comets_unmatched.hmdb <- anti_join(comets.preprocessed,
                                   moore.preprocessed,
                                   by = c("hmdb_id_final" = "input_hmdb_id_final"))

join_uid <- full_join(
  moore_unmatched.hmdb,
  comets_unmatched.hmdb %>% mutate(uid_01_final_2_comets = uid_01_final_2),
  by = c("metabolite_id_final_2" = "uid_01_final_2_comets"),
  suffix = c("_moore", "_comets")
)

join_uid <- join_uid %>%
  mutate(
    matched_uid = !is.na(.row_id_moore) & !is.na(.row_id_comets)
  )

moore_unmatched.uid <- anti_join(moore_unmatched.hmdb,
                                  comets_unmatched.hmdb,
                                  by = c("metabolite_id_final_2" = "uid_01_final_2"))

comets_unmatched.uid <- anti_join(comets_unmatched.hmdb,
                                   moore_unmatched.hmdb,
                                   by = c("uid_01_final_2" = "metabolite_id_final_2"))


join_name <- full_join(
  moore_unmatched.uid,
  comets_unmatched.uid %>% mutate(biochemical_final_final_final_comets = biochemical_final_final_final),
  by = c(
    "input_metabolite_name_final" =
      "biochemical_final_final_final_comets"
  ),
  suffix = c("_moore", "_comets")
) %>%
  mutate(
    matched_name = !is.na(.row_id_moore) & !is.na(.row_id_comets)
  )

moore_unmatched.name <- anti_join(moore_unmatched.uid,
                                  comets_unmatched.uid,
                                  by = c("input_metabolite_name_final" =
      "biochemical_final_final_final"))

comets_unmatched.name <- anti_join(comets_unmatched.uid,
                                  moore_unmatched.uid,
                                  by = c("biochemical_final_final_final" =
                                           "input_metabolite_name_final"))

matched_hmdb <- join_hmdb %>% filter(matched_hmdb) %>% select(-matched_hmdb)
matched_uid <- join_uid %>% filter(matched_uid) %>% select(-matched_uid)
matched_name <- join_name %>% filter(matched_name) %>% select(-matched_name)

comets_moore_sequential <- bind_rows(
  matched_hmdb,
  matched_uid,
  matched_name
) %>% mutate(source = "COMETS_MOORE")


comets_moore_sequential.1 <- comets_moore_sequential %>%
  distinct(.row_id_moore, .row_id_comets, .keep_all = TRUE)

moore_unmatched.name <- moore_unmatched.name %>% mutate(source = "MOORE")
comets_unmatched.name <- comets_unmatched.name %>% mutate(source = "COMETS")

all_unmatched <- full_join(moore_unmatched.name, comets_unmatched.name, by = "source")

comets_moore.merged <- rbind(comets_moore_sequential, all_unmatched) 
comets_moore_merged <- comets_moore.merged %>% select(-.row_id_moore, -.row_id_comets) %>% distinct()
