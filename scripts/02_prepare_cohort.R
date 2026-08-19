# Define KRAS groups and prepare clinical covariates.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in% c(
    "", "not reported", "unknown", "not available", "--"
  )] <- NA_character_
  x
}

collapse_stage <- function(x) {
  x <- tolower(clean_value(x))
  case_when(
    grepl("^stage iv", x) ~ "Stage IV",
    grepl("^stage iii", x) ~ "Stage III",
    grepl("^stage ii", x) ~ "Stage II",
    grepl("^stage i", x) ~ "Stage I",
    TRUE ~ NA_character_
  )
}

collapse_smoking <- function(x) {
  x <- tolower(clean_value(x))
  case_when(
    grepl("never|non-smoker", x) ~ "Never smoker",
    grepl("smoker|current|former|reformed", x) ~ "Ever smoker",
    TRUE ~ NA_character_
  )
}

protein_altering <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
  "Translation_Start_Site", "Nonstop_Mutation"
)

rna_manifest <- read_csv(
  "metadata/gdc_star_counts_manifest.csv",
  show_col_types = FALSE
) %>%
  mutate(
    aliquot_barcode = cases,
    participant_id = substr(cases, 1, 12)
  )

maf_manifest <- read_csv(
  "metadata/gdc_masked_maf_manifest.csv",
  show_col_types = FALSE
) %>%
  mutate(participant_id = substr(cases, 1, 12))

maf <- readRDS("data/processed/luad_masked_maf.rds") %>%
  mutate(participant_id = substr(Tumor_Sample_Barcode, 1, 12))

kras_variants <- maf %>%
  filter(
    Hugo_Symbol == "KRAS",
    Variant_Classification %in% protein_altering
  )

g12c <- kras_variants %>%
  filter(HGVSp_Short == "p.G12C") %>%
  pull(participant_id) %>%
  unique()

kras_mutant <- unique(kras_variants$participant_id)
mutation_available <- unique(maf_manifest$participant_id)

cohort <- rna_manifest %>%
  count(participant_id, name = "rna_file_count") %>%
  mutate(
    mutation_data_available = participant_id %in% mutation_available,
    kras_group = case_when(
      participant_id %in% g12c ~ "KRAS_G12C",
      participant_id %in% kras_mutant ~ "KRAS_other",
      mutation_data_available ~ "KRAS_wild_type",
      TRUE ~ "mutation_data_unavailable"
    )
  )

clinical <- readRDS("data/processed/luad_clinical.rds")
if (anyDuplicated(clinical$submitter_id)) {
  stop("Clinical submitter_id values are not unique.")
}

primary <- cohort %>%
  filter(kras_group %in% c("KRAS_G12C", "KRAS_wild_type"))

idx <- match(primary$participant_id, clinical$submitter_id)
age_index <- suppressWarnings(as.numeric(clean_value(clinical$age_at_index[idx])))
age_diagnosis <- suppressWarnings(
  as.numeric(clean_value(clinical$age_at_diagnosis[idx])) / 365.25
)

primary <- primary %>%
  mutate(
    age_years = coalesce(age_index, age_diagnosis),
    sex_at_birth = tolower(clean_value(clinical$sex_at_birth[idx])),
    stage_group = collapse_stage(clinical$ajcc_pathologic_stage[idx]),
    smoking_group = collapse_smoking(clinical$tobacco_smoking_status[idx]),
    pack_years_smoked = suppressWarnings(
      as.numeric(clean_value(clinical$pack_years_smoked[idx]))
    )
  )

unknown_smoking <- clean_value(clinical$tobacco_smoking_status[idx])[
  is.na(primary$smoking_group)
]
unknown_smoking <- unique(unknown_smoking[!is.na(unknown_smoking)])
if (length(unknown_smoking)) {
  warning("Unmapped smoking labels: ", paste(unknown_smoking, collapse = "; "))
}

cohort_summary <- primary %>%
  count(kras_group, name = "participants")

baseline_summary <- primary %>%
  group_by(kras_group) %>%
  summarise(
    participants = n(),
    median_age = median(age_years, na.rm = TRUE),
    female = sum(sex_at_birth == "female", na.rm = TRUE),
    stage_i = sum(stage_group == "Stage I", na.rm = TRUE),
    ever_smoker = sum(smoking_group == "Ever smoker", na.rm = TRUE),
    never_smoker = sum(smoking_group == "Never smoker", na.rm = TRUE),
    smoking_missing = sum(is.na(smoking_group)),
    .groups = "drop"
  )

observed_drivers <- maf %>%
  filter(
    Hugo_Symbol %in% c("TP53", "STK11", "KEAP1"),
    Variant_Classification %in% protein_altering,
    participant_id %in% primary$participant_id
  ) %>%
  distinct(participant_id, gene = Hugo_Symbol) %>%
  mutate(mutated = TRUE)

driver_status <- expand_grid(
  participant_id = primary$participant_id,
  gene = c("TP53", "STK11", "KEAP1")
) %>%
  left_join(observed_drivers, by = c("participant_id", "gene")) %>%
  mutate(mutated = coalesce(mutated, FALSE)) %>%
  left_join(primary %>% select(participant_id, kras_group), by = "participant_id")

driver_columns <- driver_status %>%
  select(participant_id, gene, mutated) %>%
  pivot_wider(
    names_from = gene,
    values_from = mutated,
    names_prefix = "mutated_"
  )
primary <- primary %>% left_join(driver_columns, by = "participant_id")

co_mutation_summary <- driver_status %>%
  group_by(kras_group, gene) %>%
  summarise(
    mutated = sum(mutated),
    participants = n(),
    frequency = mutated / participants,
    .groups = "drop"
  )

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(primary, "data/processed/primary_cohort.csv")
write_csv(cohort_summary, "results/tables/cohort_summary.csv")
write_csv(baseline_summary, "results/tables/baseline_summary.csv")
write_csv(co_mutation_summary, "results/tables/co_mutation_summary.csv")

message("Primary cohort:")
print(cohort_summary)
message("\nCo-mutation summary:")
print(co_mutation_summary)
