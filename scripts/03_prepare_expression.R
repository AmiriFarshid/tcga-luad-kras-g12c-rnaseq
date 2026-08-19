# Download STAR counts and create one participant-level count matrix.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
})

cohort <- read_csv("data/processed/primary_cohort.csv", show_col_types = FALSE)
manifest <- read_csv(
  "metadata/gdc_star_counts_manifest.csv",
  show_col_types = FALSE
) %>%
  mutate(participant_id = substr(cases, 1, 12)) %>%
  semi_join(cohort, by = "participant_id")

if (n_distinct(manifest$participant_id) != nrow(cohort)) {
  stop("Not every cohort participant has an RNA-seq file.")
}

message("Querying ", nrow(manifest), " STAR count files...")
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = unique(manifest$cases)
)

returned <- getResults(query)
if (!setequal(manifest$file_id, returned$file_id)) {
  stop("The GDC query does not match the saved cohort manifest.")
}

download_dir <- "data/raw/gdc_rnaseq"
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
GDCdownload(
  query,
  method = "api",
  directory = download_dir,
  files.per.chunk = 20
)

se <- GDCprepare(
  query,
  directory = download_dir,
  summarizedExperiment = TRUE
)

if (!"unstranded" %in% assayNames(se)) {
  stop("The GDC object does not contain unstranded raw counts.")
}

counts <- assay(se, "unstranded")
if (anyNA(counts) || any(counts < 0) || any(counts != floor(counts))) {
  stop("The count matrix is not a non-negative integer matrix.")
}
storage.mode(counts) <- "integer"

sample_table <- tibble(
  aliquot_barcode = colnames(counts),
  participant_id = substr(aliquot_barcode, 1, 12),
  library_size = as.numeric(colSums(counts))
) %>%
  group_by(participant_id) %>%
  arrange(desc(library_size), aliquot_barcode, .by_group = TRUE) %>%
  mutate(selected = row_number() == 1) %>%
  ungroup()

selected <- sample_table %>%
  filter(selected) %>%
  arrange(participant_id)

keep <- match(selected$aliquot_barcode, colnames(counts))
clinical_order <- match(selected$participant_id, cohort$participant_id)
if (anyNA(keep) || anyNA(clinical_order)) stop("Sample alignment failed.")

cohort <- cohort[clinical_order, ] %>%
  mutate(
    kras_group = factor(
      kras_group,
      levels = c("KRAS_wild_type", "KRAS_G12C")
    ),
    sex_at_birth = factor(sex_at_birth, levels = c("female", "male")),
    stage_group = factor(
      stage_group,
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
    ),
    smoking_group = factor(
      smoking_group,
      levels = c("Never smoker", "Ever smoker")
    ),
    library_size = selected$library_size,
    aliquot_barcode = selected$aliquot_barcode
  )
cohort <- as.data.frame(cohort)
rownames(cohort) <- cohort$aliquot_barcode

analysis_se <- SummarizedExperiment(
  assays = list(counts = counts[, keep, drop = FALSE]),
  rowData = rowData(se),
  colData = S4Vectors::DataFrame(cohort)
)
colnames(analysis_se) <- cohort$participant_id

if (anyDuplicated(colnames(analysis_se))) {
  stop("Participant identifiers are not unique.")
}

saveRDS(analysis_se, "data/processed/analysis_counts.rds", compress = "gzip")
write_csv(sample_table, "metadata/rna_aliquot_selection.csv")

message(
  "Saved ", nrow(analysis_se), " genes and ", ncol(analysis_se),
  " participants. Excluded ", sum(!sample_table$selected),
  " additional aliquots."
)
