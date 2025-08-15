library(lionessR)
library(DESeq2)
library(SummarizedExperiment)

metadata <- vroom::vroom("data/metadata.csv")
exp <- vroom::vroom("results/count_table_hisat2.txt") |>
    dplyr::distinct(Geneid, .keep_all = T) |> 
    tibble::column_to_rownames("Geneid")

colnames(exp) <- gsub(".bam", "", colnames(exp))
exp <- exp[, colnames(exp) %in% metadata$sample]

rowData <- DataFrame(row.names = rownames(exp), gene = rownames(exp))
colData <- DataFrame(row.names = metadata$sample, sample = as.character(metadata$sample), condition = metadata$condition)

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp)), 
                           colData = colData, rowData = rowData)

dds <- DESeqDataSet(se, design = ~ condition)

vst_data <- vst(dds, blind = TRUE)

# Extract the normalized data matrix
normalized_counts <- assay(vst_data)

rownames(normalized_counts)[which.max(rowSums(normalized_counts))] # q3

cat("\nDimensions of normalized data:", dim(normalized_counts), "\n")

sample_ids <- metadata$sample

if (all(sample_ids %in% colnames(normalized_counts))) {
  expr_mat <- as.matrix(normalized_counts[, sample_ids, drop = FALSE])
} else if (all(sample_ids %in% rownames(normalized_counts))) {
  # transpose if samples are in rows
  expr_mat <- t(as.matrix(normalized_counts[sample_ids, , drop = FALSE]))
} else {
  stop("Sample names in metadata$condition were not found in rows or columns of normalized_counts.")
}

lioness_networks <- lioness(expr_mat, cor)
W  <- assay(se)

flight_idx <- metadata$condition == "flight"
ground_idx <- metadata$condition == "ground"

# Per-gene t-test
pvals <- apply(W, 1, function(x) {
  # variances within each group
  var_flight <- var(x[flight_idx], na.rm = TRUE)
  var_ground <- var(x[ground_idx], na.rm = TRUE)
  
  if (var_flight == 0 || var_ground == 0) {
    return(NA_real_)  # no variance in at least one group
  } else {
    return(t.test(x[flight_idx], x[ground_idx])$p.value)
  }
})

# Adjust for multiple testing
padj <- p.adjust(pvals, method = "BH")

# Build results table
res <- data.frame(
  gene = rownames(W),
  mean_flight = rowMeans(W[, flight_idx, drop = FALSE]),
  mean_ground = rowMeans(W[, ground_idx, drop = FALSE]),
  delta = rowMeans(W[, flight_idx, drop = FALSE]) -
          rowMeans(W[, ground_idx, drop = FALSE]),
  pval = pvals,
  padj = padj
)

sig_genes <- subset(res, !is.na(padj) & padj < 0.05)
sig_genes <- sig_genes[order(sig_genes$padj), ]
nrow(sig_genes) #q4

library(clusterProfiler)
library(org.Mm.eg.db) 

sig_gene_ids <- sig_genes$gene

ego <- enrichGO(
  gene          = sig_gene_ids,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego[which.min(ego$p.adjust), "ID"] #q5
