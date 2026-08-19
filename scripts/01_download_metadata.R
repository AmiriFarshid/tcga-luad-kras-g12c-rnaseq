# Download the public metadata used to define the TCGA-LUAD cohort.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(readr)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/gdc_maf", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("metadata", recursive = TRUE, showWarnings = FALSE)

message("Querying primary-tumor STAR count files...")
rna_query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access = "open"
)

rna_manifest <- getResults(rna_query) %>%
  select(any_of(c(
    "file_id", "file_name", "cases", "sample_type", "access",
    "workflow_type", "data_type", "experimental_strategy"
  )))

if (nrow(rna_manifest) == 0) stop("No RNA-seq files were returned.")

write_csv(rna_manifest, "metadata/gdc_star_counts_manifest.csv")
saveRDS(rna_query, "data/processed/rna_query.rds")

message("Querying the GDC masked ensemble mutation MAF...")
maf_query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

write_csv(
  getResults(maf_query),
  "metadata/gdc_masked_maf_manifest.csv"
)

GDCdownload(maf_query, method = "api", directory = "data/raw/gdc_maf")
maf <- GDCprepare(
  maf_query,
  directory = "data/raw/gdc_maf",
  summarizedExperiment = FALSE
)
saveRDS(as.data.frame(maf), "data/processed/luad_masked_maf.rds")

message("Downloading indexed clinical data...")
clinical <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type = "clinical",
  save.csv = FALSE
)
saveRDS(as.data.frame(clinical), "data/processed/luad_clinical.rds")

message(
  "Retrieved ", nrow(rna_manifest), " RNA-seq files, ",
  nrow(maf), " mutation records, and ", nrow(clinical),
  " clinical records."
)
