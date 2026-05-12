library(dplyr)
library(tidyverse)
library(purrr)

# =========================
# Load MetLinkR (from GitHub)
# =========================
# Locally install RaMP
install.packages("devtools")
library(devtools)
load_all("/Users/heysiri/MetLinkR")
#install_github("ncats/RAMP-DB")

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
merged_database <- read.csv("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/output/comets_moore.merged.csv")

moore_list <- read.csv("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/fromSteveMoore_metabolite_list.csv")
comets_original <- load(file.path("/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/compileduids.RData"))

# =========================
# Build input dataframe 
# =========================
#for comets, moore separately
write.csv(comets_original, "/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/COMETS_metabolite_df.csv")
COMETS_file_name <- "/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/COMETS_metabolite_df.csv"
moore_file_name <- "/Users/heysiri/Documents/mGWAS_metabolite_harmonization/original_data/fromSteveMoore_metabolite_list.csv"
metlinkr_input_file <- data.frame("FileNames" = c(COMETS_file_name, moore_file_name),
                                  "ShortFileName" = c("COMETS", "Moore"),
                                  "HMDB" = c("hmdb_id", "input_hmdb_id"),
                                  "Metabolite_Name" = c("biochemical", "input_metabolite_name"),
                                  "PubChem_CID" = c(NA, "input_pubchem"),
                                  "KEGG" = c(NA, "input_kegg"),
                                  "LIPIDMAPS" = c(NA, NA),
                                  "chebi" = c(NA, "input_chebi"))

write.csv(metlinkr_input_file, "metlinkr_input_file_for_merging.csv")
metLinkR::harmonizeInputSheets("metlinkr_input_file_for_merging.csv")

#for merged database
merged_database$hmdb_merged <- with(
  merged_database,
  ifelse(
    is.na(input_hmdb_id_final) & is.na(hmdb_id_final),
    NA,
    ifelse(
      is.na(input_hmdb_id_final),
      hmdb_id_final,
      ifelse(
        is.na(hmdb_id_final),
        input_hmdb_id_final,
        ifelse(
          input_hmdb_id_final == hmdb_id_final,
          input_hmdb_id_final,
          paste(input_hmdb_id_final, hmdb_id_final, sep = ";")
        )
      )
    )
  )
)
merged_database <-merged_database %>% select (-input_hmdb_id_final, -hmdb_id_final)

merged_database$biochemical_merged <- with(
  merged_database,
  ifelse(
    is.na(biochemical_final_final_final) & is.na(input_metabolite_name_final),
    NA,
    ifelse(
      is.na(biochemical_final_final_final),
      input_metabolite_name_final,
      ifelse(
        is.na(input_metabolite_name_final),
        biochemical_final_final_final,
        ifelse(
          biochemical_final_final_final == input_metabolite_name_final,
          biochemical_final_final_final,
          paste(biochemical_final_final_final, input_metabolite_name_final, sep = ";")
        )
      )
    )
  )
)
merged_database <-merged_database %>% select (-biochemical_final_final_final, -input_metabolite_name_final)
merged_file_name <- "/Users/heysiri/Documents/mGWAS_metabolite_harmonization/output/comets_moore_merged_collapsed.csv"

write_csv(merged_database, merged_file_name)
merged_metlinkr_input_file <- data.frame("FileNames" = merged_file_name,
                                         "ShortFileName" = "moore_COMETS_merged",
                                         "HMDB" = "hmdb_merged",
                                         "Metabolite_Name" = "biochemical_merged",
                                         "PubChem_CID" = "input_pubchem_final",
                                         "KEGG" = "input_kegg_final",
                                         "LIPIDMAPS" = NA,
                                         "chebi" = "input_chebi_final")
write.csv(merged_metlinkr_input_file, "metlinkr_input_file_merged_for_merging.csv")
metLinkR::harmonizeInputSheets("metlinkr_input_file_merged_for_merging.csv")

