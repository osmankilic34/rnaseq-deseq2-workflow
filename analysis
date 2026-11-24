# =============================================================
# RNA-seq Differential Expression Workflow (DESeq2)
# Input : featureCounts count matrix
# Output: DE results, QC/figures and up/down-regulated gene tables
# =============================================================

# --- 1. LIBRARIES AND SETTINGS ---
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(vsn)

# Set working directory (folder containing the count file; outputs are written here)
setwd("C:/R")

# --- 2. INPUT SETTINGS (edit for each analysis) ---
dosya_adi  <- "counts_11_0vs24h.out"    # Count matrix file to read
analiz_adi <- "11_0vs24h_Analiz"        # Label used as prefix for output files

# Check that the input file exists
if (!file.exists(dosya_adi)) stop(paste("ERROR:", dosya_adi, "not found in the working directory!"))

# Read the count matrix
countData <- read.table(dosya_adi, header = TRUE, row.names = 1, sep = "\t")

# --- DATA CLEANING ---
# 1. Remove the 'LOC' prefix from gene IDs
rownames(countData) <- gsub("LOC", "", rownames(countData))

# 2. Remove file extensions from sample names
colnames(countData) <- gsub(".sorted.bam", "", colnames(countData))

# Convert to matrix
countData <- as.matrix(countData)

# --- 3. EXPERIMENTAL DESIGN ---
# Conditions must be listed in the same order as the count matrix columns
condition <- c("Control", "Control", "Exp", "Exp", "Exp")

# Safety check: stop if the number of conditions does not match the number of samples
if (length(condition) != ncol(countData)) stop("Number of conditions does not match number of samples!")
print(data.frame(Sample = colnames(countData), Condition = condition))  # check sample-condition mapping

colData <- data.frame(row.names = colnames(countData), condition = factor(condition))
colData$condition <- relevel(colData$condition, ref = "Control")

# --- 4. DIFFERENTIAL EXPRESSION ANALYSIS ---
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Exp", "Control"))

# --- 5. FIGURES ---

# A) P-value histogram
png(paste0(analiz_adi, "_1_Histogram.png"), width = 800, height = 600)
hist(res$pvalue, breaks = 50, col = "grey", main = paste("P-Value:", analiz_adi), xlab = "P", ylab = "Frekans")
dev.off()

# B) PCA (variance-stabilized counts)
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condition") + theme_bw() + ggtitle(paste("PCA:", analiz_adi))
ggsave(paste0(analiz_adi, "_2_PCA.png"), plot = pca_plot, width = 8, height = 6, dpi = 300)

# C) Volcano plot
res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  mutate(durum = case_when(
    log2FoldChange > 1 & padj < 0.05 ~ "ARTIS (UP)",
    log2FoldChange < -1 & padj < 0.05 ~ "AZALIS (DOWN)",
    TRUE ~ "ANLAMSIZ"
  ))

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = durum)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "blue", "red")) +
  theme_bw() +
  labs(title = paste("Volcano:", analiz_adi), x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
ggsave(paste0(analiz_adi, "_3_Volcano.png"), plot = volcano_plot, width = 8, height = 6, dpi = 300)

# D) Heatmap of the top 40 genes by adjusted p-value
top40_sirasi <- head(order(res$padj), 40)
heatmap_data <- assay(vsd)[top40_sirasi, ]
pheatmap(heatmap_data, cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = FALSE, show_colnames = TRUE, annotation_col = colData,
         scale = "row", fontsize_row = 10, fontsize_col = 12,
         main = paste("Heatmap Top 40:", analiz_adi),
         filename = paste0(analiz_adi, "_4_Heatmap.png"), width = 8, height = 10)

# --- 6. RESULT TABLES ---
res_export <- res_df %>%
  arrange(padj) %>%
  dplyr::select(GeneID, baseMean, log2FoldChange, lfcSE, pvalue, padj)

# All genes
write.csv2(res_export, paste0(analiz_adi, "_Tablo1_Tum_Genler.csv"), row.names = FALSE)

# Up-regulated genes (log2FC > 1, padj < 0.05)
up_genes <- res_export %>% filter(log2FoldChange > 1 & padj < 0.05)
write.csv2(up_genes, paste0(analiz_adi, "_Tablo2_UP.csv"), row.names = FALSE)

# Down-regulated genes (log2FC < -1, padj < 0.05)
down_genes <- res_export %>% filter(log2FoldChange < -1 & padj < 0.05)
write.csv2(down_genes, paste0(analiz_adi, "_Tablo3_DOWN.csv"), row.names = FALSE)

# Summary
print(paste(analiz_adi, "analysis completed!"))
print(paste("Number of UP genes:", nrow(up_genes)))
print(paste("Number of DOWN genes:", nrow(down_genes)))
