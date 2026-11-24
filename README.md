# rnaseq-deseq2-workflow
Automated R workflow for RNA-seq differential expression: DESeq2 analysis, QC and publication-ready plots (PCA, volcano, heatmap), and up/down-regulated gene tables
# RNA-seq Differential Expression Workflow (DESeq2)

An automated R workflow that takes a gene count matrix and produces differential expression results, QC and publication-ready figures, and up/down-regulated gene tables in a single run.

Developed for my M.Sc. thesis on biotic and abiotic stress-responsive genes in barley (*Hordeum vulgare*). Upstream processing (SRA Toolkit, FastQC, fastp, MultiQC, HISAT2, SAMtools, featureCounts) was performed on the Linux command line; selected differentially expressed genes were validated by RT-qPCR.

## What it does

1. Reads a featureCounts count matrix and cleans gene and sample names
2. Runs differential expression analysis with DESeq2 (Exp vs Control)
3. Generates figures:
   - P-value histogram
   - PCA plot (VST-transformed counts)
   - Volcano plot
   - Heatmap of the top 40 genes by adjusted p-value
4. Exports results tables:
   - All genes, sorted by adjusted p-value
   - Up-regulated genes (log2FC > 1, padj < 0.05)
   - Down-regulated genes (log2FC < -1, padj < 0.05)

## Requirements

R (≥ 4.0) with the following packages:

- Bioconductor: `DESeq2`, `vsn`
- CRAN: `tidyverse`, `pheatmap`, `ggrepel`

## Configuration

Edit these values at the top of the script for each analysis:

| Setting | Description |
| --- | --- |
| `setwd()` | Folder containing the count file; outputs are written here |
| `dosya_adi` | Count matrix file name (featureCounts output) |
| `analiz_adi` | Label used as a prefix for all output files |
| `condition` | Sample conditions, **in the same order as the count matrix columns** |

The script stops if the number of conditions does not match the number of samples, and prints the sample–condition mapping so it can be checked before analysis.

## Usage

Rscript deseq2_analysis.R

or open the script in RStudio and click **Source**.

## Output

All files are prefixed with `analiz_adi`:

- `_1_Histogram.png`, `_2_PCA.png`, `_3_Volcano.png`, `_4_Heatmap.png`
- `_Tablo1_Tum_Genler.csv` (all genes), `_Tablo2_UP.csv`, `_Tablo3_DOWN.csv`

## Planned improvements

- Port the full workflow, including upstream steps, to Nextflow
- Read sample conditions from a metadata file instead of editing the script
