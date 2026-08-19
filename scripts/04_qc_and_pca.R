# Calculate sample-level QC metrics and inspect global expression patterns.

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

robust_z <- function(x) {
  spread <- mad(x, na.rm = TRUE)
  if (is.na(spread) || spread == 0) return(rep(0, length(x)))
  (x - median(x, na.rm = TRUE)) / spread
}

se <- readRDS("data/processed/analysis_counts.rds")
counts <- assay(se, "counts")
gene_names <- as.character(rowData(se)$gene_name)
mitochondrial <- !is.na(gene_names) & grepl("^MT-", gene_names)

qc <- tibble(
  participant_id = colnames(se),
  kras_group = as.character(colData(se)$kras_group),
  library_size = as.numeric(colSums(counts)),
  detected_genes = as.numeric(colSums(counts >= 10)),
  mitochondrial_percent = as.numeric(
    100 * colSums(counts[mitochondrial, , drop = FALSE]) / colSums(counts)
  )
) %>%
  mutate(
    library_outlier = abs(robust_z(log10(library_size))) > 3,
    detected_gene_outlier = abs(robust_z(detected_genes)) > 3,
    mitochondrial_outlier = robust_z(mitochondrial_percent) > 3,
    qc_flag = library_outlier | detected_gene_outlier | mitochondrial_outlier
  )

colData(se)$qc_flag <- qc$qc_flag
saveRDS(se, "data/processed/analysis_counts.rds", compress = "gzip")

qc_summary <- qc %>%
  group_by(kras_group) %>%
  summarise(
    participants = n(),
    median_library_size = median(library_size),
    median_detected_genes = median(detected_genes),
    median_mitochondrial_percent = median(mitochondrial_percent),
    flagged_samples = sum(qc_flag),
    .groups = "drop"
  )

keep_genes <- rowSums(counts >= 10) >= 10
dds <- DESeqDataSetFromMatrix(
  countData = counts[keep_genes, , drop = FALSE],
  colData = as.data.frame(colData(se)),
  design = ~ 1
)
vst_matrix <- assay(vst(dds, blind = TRUE))
variances <- matrixStats::rowVars(vst_matrix, useNames = FALSE)
top <- order(variances, decreasing = TRUE)[seq_len(min(500, length(variances)))]
pca <- prcomp(t(vst_matrix[top, , drop = FALSE]))
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- as_tibble(pca$x[, 1:2], rownames = "participant_id") %>%
  left_join(
    as.data.frame(colData(se)) %>%
      select(participant_id, kras_group, smoking_group),
    by = "participant_id"
  ) %>%
  mutate(smoking_group = coalesce(as.character(smoking_group), "Missing"))

colors <- c(KRAS_wild_type = "#4C78A8", KRAS_G12C = "#E45756")
pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, color = kras_group, shape = smoking_group)
) +
  geom_point(alpha = 0.75, size = 2.2) +
  scale_color_manual(values = colors) +
  labs(
    title = "PCA of TCGA-LUAD RNA-seq samples",
    subtitle = "VST; 500 most variable genes",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    color = "KRAS group",
    shape = "Smoking status"
  ) +
  theme_minimal(base_size = 12)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
write_csv(qc_summary, "results/tables/qc_summary.csv")
write_csv(qc, "data/processed/sample_qc_metrics.csv")
ggsave(
  "results/figures/pca.png",
  pca_plot,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

message("Genes used for VST: ", sum(keep_genes))
message("QC-flagged samples: ", sum(qc$qc_flag))
print(qc_summary)
