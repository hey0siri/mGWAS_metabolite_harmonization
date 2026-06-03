library(dplyr)
library(tidyverse)
library(purrr)

base_dir <- getwd()

# =========================
# Load MetLinkR (from GitHub)
# =========================
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

# =========================
# Load data
# =========================
#merged_database <- read.csv("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/output/comets_moore.merged.csv")

moore_list.original <- read.csv("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/fromSteveMoore_metabolite_list.csv")
comets.original <- load(file.path("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/compileduids.RData"))

moore_list.preprocessed <- read.csv(file.path(base_dir, "output/moore_list_preprocessed_revised.csv"))
comets.preprocessed <- read.csv(file.path(base_dir, "output/comets_preprocessed_revised.csv"))

# =========================
# Build input dataframe 
# =========================
#running comets and moore jointly (as a test)
COMETS_file_name <- file.path(base_dir, "output/comets_preprocessed_revised.csv")
moore_file_name <- file.path(base_dir, "output/moore_list_preprocessed_revised.csv")
metlinkr_input_file <- data.frame("FileNames" = c(COMETS_file_name, moore_file_name),
                                  "ShortFileName" = c("COMETS", "Moore"),
                                  "HMDB" = c("hmdb_id", "input_hmdb_id"),
                                  "Metabolite_Name" = c("biochemical_final", "input_metabolite_name_final"),
                                  "PubChem_CID" = c(NA, "input_pubchem"),
                                  "KEGG" = c(NA, "input_kegg"),
                                  "LIPIDMAPS" = c(NA, NA),
                                  "chebi" = c(NA, "input_chebi"))

write.csv(metlinkr_input_file, file.path(base_dir, "metlinkr_input_files/metlinkr_input_file_for_merging.csv"))

#metLinkR::harmonizeInputSheets(file.path(base_dir, "metlinkr_input_files/metlinkr_input_file_for_merging.csv"))
mapping_library_path <- file.path(base_dir, "metLinkR_output/mapping_library.xlsx")


if (!file.exists(mapping_library_path)) {
  
  metLinkR::harmonizeInputSheets(
    
    file.path(base_dir, "metlinkr_input_files/metlinkr_input_file_for_merging.csv")
    
  )
  
} else {
  
  message("Using existing metLinkR_output/mapping_library.xlsx; skipping MetLinkR rerun.")
  
}

# ================================
# Output Sanity Check + COMETS Merging 
# ================================
metlinkr_output_file <- readxl::read_xlsx(file.path(base_dir, "metLinkR_output/mapping_library.xlsx"))
comets_columns <- c("biochemical_final", "hmdb_id")

metlinkr_lookup.comets <- metlinkr_output_file %>%
  select(
    Harmonized_Name = `Harmonized name`, 
    match_value = `Input name (COMETS)`) %>%
  filter(match_value != "-") %>% 
  separate_rows(
    match_value,
    sep = "\\s*;\\s*"
  ) %>%
  mutate(
    match_value = stringr::str_trim(match_value)
  ) %>%
  distinct()

comets.hmdb_merged <- comets.preprocessed %>%
  left_join(
    metlinkr_lookup.comets,
    by = c("hmdb_id" = "match_value")
  ) %>%
  rename(HMDB_Harmonized_Name = Harmonized_Name)


comets.biochemical_merged <- comets.hmdb_merged %>%
  left_join(
    metlinkr_lookup.comets,
    by = c("biochemical_final" = "match_value")
  ) %>%
  rename(Biochemical_final_Harmonized_Name = Harmonized_Name) %>%
  left_join(
    metlinkr_lookup.comets,
    by = c("biochemical" = "match_value")
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

comets.metlinkr <- comets.biochemical_merged %>%
  mutate(
    MLR_Harmonized_Name = coalesce(
      HMDB_Harmonized_Name,
      Biochemical_Harmonized_Name
    )
  ) %>%
  select(-HMDB_Harmonized_Name, -Biochemical_Harmonized_Name) 


# ================================
# Output Sanity Check + MetLinkR Merging 
# ================================
id_cols <- c(
  "input_hmdb_id",
  "input_pubchem",
  "input_chemspider",
  "input_kegg",
  "input_chebi"
)

moore_list.preprocessed <- moore_list.preprocessed %>%
  mutate(across(all_of(id_cols), as.character))

metlinkr_lookup.moore <- metlinkr_output_file %>%
  select(
    Harmonized_Name = `Harmonized name`, 
    match_value = `Input name (Moore)`) %>%
  filter(match_value != "-") %>% 
  separate_rows(
    match_value,
    sep = "\\s*;\\s*"
  ) %>%
  mutate(
    match_value = stringr::str_trim(match_value)
  ) %>%
  distinct()

moore.hmdb_merged <- moore_list.preprocessed %>%
  left_join(
    metlinkr_lookup.moore,
    by = c("input_hmdb_id" = "match_value")
  ) %>%
  rename(HMDB_Harmonized_Name = Harmonized_Name)

moore.pubchem_merged <- moore.hmdb_merged %>% 
  left_join(
    metlinkr_lookup.moore,
    by = c("input_pubchem" = "match_value")
  ) %>%
  rename(Pubchem_Harmonized_Name = Harmonized_Name)

moore.chemspider_merged <- moore.pubchem_merged %>% 
  left_join(
    metlinkr_lookup.moore,
    by = c("input_chemspider" = "match_value")
  ) %>%
  rename(Chemspider_Harmonized_Name = Harmonized_Name)

moore.kegg_merged <- moore.chemspider_merged %>% 
  left_join(
    metlinkr_lookup.moore,
    by = c("input_kegg" = "match_value")
  ) %>%
  rename(Kegg_Harmonized_Name = Harmonized_Name)

moore.chebi_merged <- moore.kegg_merged %>% 
  left_join(
    metlinkr_lookup.moore,
    by = c("input_chebi" = "match_value")
  ) %>%
  rename(Chebi_Harmonized_Name = Harmonized_Name)


moore.inchikey_merged <- moore.chebi_merged %>% 
  left_join(
    metlinkr_lookup.moore,
    by = c("input_inchikey" = "match_value")
  ) %>%
  rename(Inchikey_Harmonized_Name = Harmonized_Name)



moore.biochemical_merged <- moore.inchikey_merged %>%
  left_join(
    metlinkr_lookup.moore,
    by = c("input_metabolite_name_final" = "match_value")
  ) %>%
  rename(Biochemical_final_Harmonized_Name = Harmonized_Name) %>%
  left_join(
    metlinkr_lookup.moore,
    by = c("input_metabolite_name" = "match_value")
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

moore.metlinkr <- moore.biochemical_merged %>%
  mutate(
    MLR_Harmonized_Name = coalesce(
      HMDB_Harmonized_Name,
      Biochemical_Harmonized_Name,
      Inchikey_Harmonized_Name,
      Chebi_Harmonized_Name, 
      Kegg_Harmonized_Name,
      Chemspider_Harmonized_Name,
      Pubchem_Harmonized_Name
    )
  ) %>%
  select(-HMDB_Harmonized_Name, -Biochemical_Harmonized_Name, -Inchikey_Harmonized_Name, -Chebi_Harmonized_Name, -Kegg_Harmonized_Name, -Chemspider_Harmonized_Name, -Pubchem_Harmonized_Name) 

write.csv(comets.metlinkr, file.path(base_dir, "/output/comets_metlinkr_merged.csv"))
write.csv(moore.metlinkr, file.path(base_dir, "/output/moore_metlinkr_merged.csv"))
