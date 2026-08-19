# Fit the primary DESeq2 model and two sensitivity models.

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

fit_deseq2 <- function(se, samples, design_formula) {
  metadata <- droplevels(as.data.frame(colData(se)[samples, ]))
  metadata$age_per_10_years <-
    (metadata$age_years - mean(metadata$age_years)) / 10
  metadata$stage_group <- factor(
    make.names(as.character(metadata$stage_group)),
    levels = make.names(levels(metadata$stage_group))
  )
  metadata$smoking_group <- factor(
    make.names(as.character(metadata$smoking_group)),
    levels = make.names(levels(metadata$smoking_group))
  )

  counts <- assay(se, "counts")[, samples, drop = FALSE]
  keep <- rowSums(counts >= 10) >= 10

  if (qr(model.matrix(design_formula, metadata))$rank !=
      ncol(model.matrix(design_formula, metadata))) {
    stop("The DESeq2 design matrix is not full rank.")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts[keep, , drop = FALSE],
    colData = metadata,
    design = design_formula
  )
  # Use the explicit DESeq2 workflow to allow a larger iteration limit.
  # Cook's outliers remain filtered by results() rather than replaced.
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds, maxit = 1000)

  non_converged <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
  if (non_converged > 0) {
    stop(non_converged, " genes did not converge in the DESeq2 model.")
  }

  raw <- results(
    dds,
    contrast = c("kras_group", "KRAS_G12C", "KRAS_wild_type")
  )
  coefficient <- grep(
    "kras_group.*KRAS_G12C",
    resultsNames(dds),
    value = TRUE
  )
  if (length(coefficient) != 1) stop("KRAS coefficient was not identified.")

  shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm")
  genes <- as.data.frame(rowData(se))[match(rownames(dds), rownames(se)), ]

  table <- tibble(
    ensembl_gene_id = rownames(dds),
    gene_name = genes$gene_name,
    gene_type = genes$gene_type,
    base_mean = raw$baseMean,
    raw_log2_fold_change = raw$log2FoldChange,
    shrunken_log2_fold_change = shrunk$log2FoldChange,
    statistic = raw$stat,
    p_value = raw$pvalue,
    adjusted_p_value = raw$padj
  ) %>%
    arrange(adjusted_p_value)

  list(dds = dds, table = table)
}

se <- readRDS("data/processed/analysis_counts.rds")
metadata <- as.data.frame(colData(se))

primary_complete <- complete.cases(
  metadata[, c(
    "age_years", "sex_at_birth", "stage_group",
    "smoking_group", "kras_group"
  )]
)
primary <- fit_deseq2(
  se,
  primary_complete,
  ~ age_per_10_years + sex_at_birth + stage_group + smoking_group + kras_group
)

ever_smoker_complete <- !is.na(metadata$smoking_group) &
  metadata$smoking_group == "Ever smoker" &
  complete.cases(
    metadata[, c("age_years", "sex_at_birth", "stage_group", "kras_group")]
  )
ever_smoker <- fit_deseq2(
  se,
  ever_smoker_complete,
  ~ age_per_10_years + sex_at_birth + stage_group + kras_group
)

co_mutation_complete <- complete.cases(
  metadata[, c(
    "age_years", "sex_at_birth", "stage_group", "smoking_group",
    "mutated_TP53", "mutated_STK11", "mutated_KEAP1", "kras_group"
  )]
)
co_mutation <- fit_deseq2(
  se,
  co_mutation_complete,
  ~ age_per_10_years + sex_at_birth + stage_group + smoking_group +
    mutated_TP53 + mutated_STK11 + mutated_KEAP1 + kras_group
)

alpha <- 0.05
lfc_threshold <- 1
classify <- function(x) {
  x %>%
    mutate(
      selected = !is.na(adjusted_p_value) &
        !is.na(shrunken_log2_fold_change) &
        adjusted_p_value < alpha &
        abs(shrunken_log2_fold_change) >= lfc_threshold,
      direction = case_when(
        selected & shrunken_log2_fold_change > 0 ~ "Higher in G12C",
        selected & shrunken_log2_fold_change < 0 ~ "Lower in G12C",
        TRUE ~ "Not selected"
      )
    )
}

primary_table <- classify(primary$table)
ever_table <- classify(ever_smoker$table)
co_mutation_table <- classify(co_mutation$table)

compare_with_primary <- function(sensitivity_table, model_name) {
  comparison <- primary_table %>%
    select(
      ensembl_gene_id,
      primary_lfc = shrunken_log2_fold_change,
      primary_selected = selected,
      primary_direction = direction
    ) %>%
    inner_join(
      sensitivity_table %>%
        select(
          ensembl_gene_id,
          sensitivity_lfc = shrunken_log2_fold_change,
          sensitivity_selected = selected,
          sensitivity_direction = direction
        ),
      by = "ensembl_gene_id"
    )

  overlap <- comparison %>% filter(primary_selected, sensitivity_selected)
  tibble(
    sensitivity_model = model_name,
    spearman_lfc_correlation = cor(
      comparison$primary_lfc,
      comparison$sensitivity_lfc,
      method = "spearman",
      use = "complete.obs"
    ),
    selected_gene_overlap = nrow(overlap),
    direction_agreement_in_overlap = if (nrow(overlap)) {
      mean(overlap$primary_direction == overlap$sensitivity_direction)
    } else {
      NA_real_
    }
  )
}

sensitivity_summary <- tibble(
  model = c(
    "Primary adjusted",
    "Ever smokers only",
    "Adjusted for TP53, STK11, and KEAP1"
  ),
  participants = c(
    sum(primary_complete),
    sum(ever_smoker_complete),
    sum(co_mutation_complete)
  ),
  g12c_participants = c(
    sum(metadata$kras_group[primary_complete] == "KRAS_G12C"),
    sum(metadata$kras_group[ever_smoker_complete] == "KRAS_G12C"),
    sum(metadata$kras_group[co_mutation_complete] == "KRAS_G12C")
  ),
  selected_genes = c(
    sum(primary_table$selected),
    sum(ever_table$selected),
    sum(co_mutation_table$selected)
  ),
  higher_in_g12c = c(
    sum(primary_table$direction == "Higher in G12C"),
    sum(ever_table$direction == "Higher in G12C"),
    sum(co_mutation_table$direction == "Higher in G12C")
  ),
  lower_in_g12c = c(
    sum(primary_table$direction == "Lower in G12C"),
    sum(ever_table$direction == "Lower in G12C"),
    sum(co_mutation_table$direction == "Lower in G12C")
  )
)

concordance <- bind_rows(
  compare_with_primary(ever_table, "Ever smokers only"),
  compare_with_primary(
    co_mutation_table,
    "Adjusted for TP53, STK11, and KEAP1"
  )
)

volcano_data <- primary_table %>%
  filter(!is.na(adjusted_p_value), !is.na(shrunken_log2_fold_change)) %>%
  mutate(
    plot_padj = pmax(adjusted_p_value, .Machine$double.xmin),
    direction = factor(
      direction,
      levels = c("Not selected", "Higher in G12C", "Lower in G12C")
    )
  )

colors <- c(
  "Not selected" = "grey75",
  "Higher in G12C" = "#E45756",
  "Lower in G12C" = "#4C78A8"
)
volcano <- ggplot(
  volcano_data,
  aes(shrunken_log2_fold_change, -log10(plot_padj), color = direction)
) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(alpha), linetype = 2, color = "grey40") +
  scale_color_manual(values = colors) +
  labs(
    title = "KRAS G12C-associated differential expression",
    subtitle = "Covariate-adjusted DESeq2 model",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write_csv(primary_table, "data/processed/primary_de_results.csv")
write_csv(ever_table, "data/processed/ever_smoker_de_results.csv")
write_csv(
  co_mutation_table,
  "data/processed/co_mutation_adjusted_de_results.csv"
)
write_csv(
  primary_table %>% filter(selected) %>% slice_head(n = 100),
  "results/tables/differential_expression_top100.csv"
)
write_csv(sensitivity_summary, "results/tables/sensitivity_summary.csv")
write_csv(concordance, "results/tables/sensitivity_concordance.csv")
ggsave(
  "results/figures/volcano.png",
  volcano,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

message("Model summary:")
print(sensitivity_summary)
message("\nSensitivity concordance:")
print(concordance)
