# Metabolite Harmonization for Metabolite GWAS

> **Status:** ⚠️ Work in progress
>
> This repository contains an evolving metabolite harmonization workflow for integrating metabolite annotations across multiple metabolomics resources. Scripts, matching logic, and outputs are under active development and should be reviewed before production use.

---

# Repository Purpose

This repository implements a metabolite harmonization pipeline designed to:

* Standardize metabolite identifiers and names from multiple reference resources.
* Resolve internal inconsistencies within source metabolite lists.
* Incorporate external harmonization from MetLinkR / RaMP.
* Generate a merged COMETS–Moore reference database for downstream metabolite GWAS harmonization.

The current workflow focuses on integrating:

* COMETS metabolite annotations
* Dr. Steven Moore metabolite annotations
* MetLinkR harmonization results

---

# Input Data

The `original_data/` directory contains the source files used by the pipeline:

```
original_data/
├── compileduids.RData
├── fromSteveMoore_metabolite_list.csv
└── rows_split.csv
```

## Source Descriptions

### COMETS

Primary source:

```
compileduids.RData
```

Contains COMETS metabolite identifiers, HMDB mappings, and biochemical names.

### Moore Metabolite List

```
fromSteveMoore_metabolite_list.csv
```

Contains metabolite identifiers collected and curated by Dr. Steven Moore.

### Supporting Files

```
rows_split.csv
```

File generated from Claude for the splitting of concatenated metabolite identifiers curated by Dr. Steven Moore.

---

# Workflow Overview

The harmonization pipeline proceeds through four major stages:

```
COMETS
        \
         \
          --> Preprocessing --> MetLinkR Harmonization --> COMETS-Moore Merge --> Reference Harmonization Database
         /                                     
Moore -- 
                                                 
```

Scripts should generally be executed in the following order.

---

# 1. COMETS Preprocessing

Script:

```
scripts/preprocess_comets.R
```

## Purpose

Standardize and expand COMETS metabolite annotations before harmonization.

## Operations Performed

### Identifier Expansion

Handles metabolites represented by multiple identifiers, including:

* `#`-concatenated identifiers
* biochemical names joined with `"or"`

Examples:

```
HMDB00001#HMDB00002
```

becomes

```
HMDB00001
HMDB00002
```

while preserving provenance.

### Variant Handling

Investigates and optionally collapses identifier variants such as:

```
_UID_B1
_UID_B2
```

when annotations are consistent.

### Name Standardization

Applies biochemical name normalization through helper utilities.

### HMDB Standardization

Normalizes HMDB identifiers into a consistent format.

## Output

```
output/comets_preprocessed_revised.csv
```

---

# 2. Moore Preprocessing

Script:

```
scripts/preprocess_moore.R
```

## Purpose

Standardize metabolite identifiers and names from the Moore metabolite list.

## Operations Performed

* Identifier cleaning
* HMDB normalization
* Metabolite name normalization
* Detection of multiple identifiers
* Detection of isomer annotations
* Provenance tracking

## Output

```
output/moore_list_preprocessed_revised.csv
```

---

# 3. MetLinkR Harmonization

Script:

```
scripts/moore_comets_metlinkr.R
```

## Purpose

Assign standardized metabolite names using MetLinkR.

## External Dependencies

* MetLinkR
* RaMP

The script generates a MetLinkR input specification and runs harmonization if a previous mapping library does not already exist. This script depends on a particular fork of the MetLinkR repo.

### Generated Input

```
metlinkr_input_files/
└── metlinkr_input_file_for_merging.csv
```

### MetLinkR Output

```
metLinkR_output/
└── mapping_library.xlsx
```

If `mapping_library.xlsx` already exists, the MetLinkR step is skipped and the existing file is reused.

---

## Harmonization Priority

### COMETS

Matches are attempted using:

1. HMDB ID
2. Normalized biochemical name
3. Original biochemical name

### Moore

Matches are attempted using:

1. HMDB ID
2. Metabolite name
3. InChIKey
4. ChEBI
5. KEGG
6. ChemSpider
7. PubChem CID

The first available harmonized name is retained as:

```
MLR_Harmonized_Name
```

---

## Outputs

```
output/comets_metlinkr_merged.csv
output/moore_metlinkr_merged.csv
```

---

# 4. COMETS–Moore Reference Database Construction

Script:

```
scripts/merge_comets_moore.R
```

## Purpose

Construct a merged reference harmonization database linking COMETS and Moore metabolites.

## Sequential Matching Strategy

Matches are attempted in the following order:

### Step 1: HMDB

```
input_hmdb_id ↔ hmdb_id
```

### Step 2: Metabolite UID

```
metabolite_id ↔ uid_01
```

### Step 3: Standardized Metabolite Name

```
input_metabolite_name_final
    ↔
biochemical_final
```

### Step 4: MetLinkR Harmonized Name

```
MLR_Harmonized_Name
```

Additional exploratory matching strategies are included but currently marked as experimental.

---

## Join Characteristics

* Many-to-many joins are permitted.
* Unmatched records are propagated to subsequent matching stages.
* Provenance columns are retained throughout processing.

---

## Output

```
output/draft_comets_moore_merged.csv
```

This file currently serves as the primary merged reference database.

---

# Helper Utilities

Script:

```
scripts/helpers_metabolite_utils.R
```

Provides reusable functions for:

* HMDB standardization
* Biochemical name normalization
* Delimiter expansion
* Variant detection
* Provenance tracking
* Identifier consistency checks
* Duplicate collapsing utilities

Many preprocessing steps depend on these helper functions.

---

# Harmonization Philosophy

The workflow prioritizes structured identifiers over free-text names whenever possible.

Approximate priority:

1. HMDB
2. Metabolite UID
3. PubChem CID
4. InChIKey
5. ChEBI
6. KEGG
7. ChemSpider
8. Harmonized MetLinkR name
9. Standardized metabolite name
10. Original metabolite name

Higher-confidence identifier matches are preferred before lower-confidence name-based matching.

---

# Key Outputs

## Preprocessed Data

```
output/comets_preprocessed_revised.csv
output/moore_list_preprocessed_revised.csv
```

## MetLinkR Harmonization

```
output/comets_metlinkr_merged.csv
output/moore_metlinkr_merged.csv
```

## Merged Reference Database

```
output/draft_comets_moore_merged.csv
```

---

# Known Limitations

* Many-to-many matches are intentionally retained and require interpretation.
* Some identifier expansion logic (e.g. underscore-delimited COMETS identifiers) remains under development.
* Several exploratory merge strategies are present but not yet finalized.
* Current ouput of merged Moore's and COMETS databases will result to many-to-many matches due to many duplicates. Will be revised with experimenting harmonization.

Results should always be reviewed for biological plausibility before downstream analysis.

---

*Last updated: June 2026*
