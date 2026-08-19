# Run preranked GSEA with MSigDB Hallmark gene sets.

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

format_pathway <- function(x) {
  x <- tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))
  replacements <- c(
    "Tnfa" = "TNFα", "Nfkb" = "NF-κB", "Tgf Beta" = "TGF-β",
    "E2f" = "E2F", "G2m" = "G2/M", "Mtorc1" = "mTORC1",
    "Il2 Stat5" = "IL-2/STAT5", "Kras" = "KRAS", "Uv" = "UV",
    "P53" = "p53"
  )
  for (pattern in names(replacements)) {
    x <- gsub(pattern, replacements[[pattern]], x, fixed = TRUE)
  }
  x
}

de <- read_csv(
  "data/processed/primary_de_results.csv",
  show_col_types = FALSE
) %>%
  mutate(ensembl_gene = sub("\\..*$", "", ensembl_gene_id)) %>%
  filter(!is.na(p_value), !is.na(statistic), is.finite(statistic)) %>%
  group_by(ensembl_gene) %>%
  slice_max(abs(statistic), n = 1, with_ties = FALSE) %>%
  ungroup()

ranks <- de$statistic
names(ranks) <- de$ensembl_gene
ranks <- sort(ranks, decreasing = TRUE)

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)
pathways <- split(hallmark$ensembl_gene, hallmark$gs_name)
pathways <- lapply(pathways, unique)

set.seed(2026)
gsea <- fgseaMultilevel(
  pathways = pathways,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
) %>%
  as_tibble() %>%
  arrange(padj) %>%
  transmute(
    pathway,
    pathway_label = format_pathway(pathway),
    direction = if_else(NES > 0, "Higher in G12C", "Lower in G12C"),
    size,
    enrichment_score = ES,
    normalized_enrichment_score = NES,
    p_value = pval,
    adjusted_p_value = padj
  )

significant <- gsea %>% filter(adjusted_p_value < 0.05)
if (nrow(significant) == 0) stop("No significant Hallmark pathways were found.")
plot_data <- significant %>%
  slice_min(
    adjusted_p_value,
    n = min(20, nrow(significant)),
    with_ties = FALSE
  ) %>%
  arrange(normalized_enrichment_score) %>%
  mutate(
    pathway_label = factor(pathway_label, levels = pathway_label),
    significance = -log10(pmax(adjusted_p_value, .Machine$double.xmin))
  )

colors <- c("Higher in G12C" = "#E45756", "Lower in G12C" = "#4C78A8")
gsea_plot <- ggplot(
  plot_data,
  aes(
    normalized_enrichment_score,
    pathway_label,
    color = direction,
    size = significance
  )
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = colors) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "Hallmark pathway enrichment",
    subtitle = "Preranked GSEA using DESeq2 Wald statistics",
    x = "Normalized enrichment score",
    y = NULL,
    color = NULL,
    size = "-log10(padj)"
  ) +
  theme_minimal(base_size = 12)

write_csv(gsea, "results/tables/hallmark_gsea.csv")
ggsave(
  "results/figures/hallmark_gsea.png",
  gsea_plot,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

message(
  nrow(significant), " of ", nrow(gsea),
  " Hallmark pathways had adjusted p-value < 0.05."
)
