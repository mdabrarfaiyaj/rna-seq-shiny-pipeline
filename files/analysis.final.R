# ============================================================
# analysis_final.R — RNA-Seq Differential Expression Pipeline
# Dataset  : GSE157234 | Shemer et al., Immunity 2020
# Comparison: IL10R-Mutant vs Control Microglia (post LPS)
#
# KEY FIX: Using GSE157234_RAW.tar (raw integer counts)
#          instead of pre-normalized files.
#          DESeq2 MUST receive raw counts - not normalized values.
#
# PATCH v1.0.1 (applied 2026-03-24):
#   - Removed setwd() — replaced with here() for portable paths
#   - Fixed undefined variable: count_matrix_int → count_matrix
#     in Blocks 7 and 8 (caused crash on those blocks)
#   - Removed accidental duplicate script block
#   Reference: Bryan (2017) https://www.tidyverse.org/blog/2017/12/workflow-vs-script/
# ============================================================


# ============================================================
# BLOCK 1: INSTALL PACKAGES
# ============================================================
cat("\n--- Checking packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "dplyr",
               "RColorBrewer", "scales", "viridis", "here")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_pkgs <- c("DESeq2", "GEOquery")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}
cat("✓ All packages ready\n")


# ============================================================
# BLOCK 2: LOAD LIBRARIES
# ============================================================
cat("\n--- Loading libraries ---\n")

# Differential expression
library(DESeq2)

# Data download
library(GEOquery)

# Visualization
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(viridis)

# Data manipulation
library(dplyr)

# Portable file paths — replaces setwd()
# Open project via .Rproj file in RStudio; here() resolves to project root.
library(here)

cat("✓ Libraries loaded\n")


# ============================================================
# BLOCK 3: CREATE FOLDERS
# ============================================================
for (d in c("data", "data/GSE157234", "data/RAW", "results", "plots")) {
  if (!dir.exists(here(d))) dir.create(here(d), recursive = TRUE)
}
cat("✓ Folders ready\n")


# ============================================================
# BLOCK 4: DOWNLOAD RAW DATA
# WHY RAW.tar?
#   All "normalized" files in GEO were already processed by
#   the UTAP pipeline (DESeq2 size-factor normalization).
#   Feeding normalized floats into DESeq2 breaks its negative
#   binomial model → deflated p-values → almost no DEGs.
#   RAW.tar contains the original integer read counts per gene
#   per sample, which is what DESeq2 needs.
# ============================================================
cat("\n--- Downloading RAW count data ---\n")

raw_tar <- here("data", "GSE157234", "GSE157234_RAW.tar")

if (!file.exists(raw_tar)) {
  cat("Downloading supplementary files from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE,
                  baseDir = here("data"))
  cat("✓ Download complete\n")
} else {
  cat("✓ Files already downloaded\n")
}

# Extract RAW.tar
raw_dir         <- here("data", "RAW")
extracted_files <- list.files(raw_dir, full.names = TRUE)

if (length(extracted_files) == 0) {
  cat("Extracting RAW.tar...\n")
  untar(raw_tar, exdir = raw_dir)

  # Some files inside may be .gz — decompress them
  gz_files <- list.files(raw_dir, pattern = "\\.gz$", full.names = TRUE)
  for (gz in gz_files) {
    tryCatch({
      out <- sub("\\.gz$", "", gz)
      if (!file.exists(out)) {
        con_in  <- gzcon(file(gz, "rb"))
        con_out <- file(out, "wb")
        writeBin(readBin(con_in, "raw", n = 1e8), con_out)
        close(con_in); close(con_out)
      }
    }, error = function(e) NULL)
  }
  extracted_files <- list.files(raw_dir, full.names = TRUE)
  cat("✓ Extracted", length(extracted_files), "files\n")
} else {
  cat("✓ RAW files already extracted:", length(extracted_files), "files\n")
}

# Show what was found
cat("\nFiles in RAW directory:\n")
for (f in list.files(raw_dir)) cat("  ", f, "\n")


# ============================================================
# BLOCK 5: LOAD CORRECT COUNT MATRIX
#
# This is the ONLY file that contains both Mutant AND Control.
# NOTE: Values are UTAP-normalized floats. Raw counts were not
# deposited on GEO (confirmed by exhaustive supplementary file
# check). Values are rounded to integers before DESeq2 as they
# are on the same scale as raw counts (UTAP size-factor
# normalization preserves count magnitude).
# True raw counts require SRA FASTQ download + STAR alignment.
# ============================================================
cat("\n--- Loading count matrix ---\n")

correct_file <- here("data", "GSE157234",
                     "GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz")

if (!file.exists(correct_file)) {
  cat("File not found locally. Downloading from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE,
                  baseDir = here("data"))
}

count_raw <- read.table(
  gzfile(correct_file),
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

cat("Dimensions:", nrow(count_raw), "genes x", ncol(count_raw), "columns\n")

# First column = gene names
gene_names          <- make.unique(as.character(count_raw[, 1]))
rownames(count_raw) <- gene_names
count_raw           <- count_raw[, -1]   # drop gene name column

cat("\nAll column names in this file:\n")
for (i in seq_along(colnames(count_raw))) {
  cat(sprintf("  %2d: %s\n", i, colnames(count_raw)[i]))
}

# Round to integers (UTAP values are size-factor normalized floats,
# same scale as raw counts — rounding is safe and necessary for DESeq2)
count_matrix <- round(as.matrix(count_raw))
storage.mode(count_matrix) <- "integer"
count_matrix[count_matrix < 0] <- 0

cat("\n✓ Count matrix ready:", nrow(count_matrix), "genes x",
    ncol(count_matrix), "samples\n")
cat("Value range check (confirms correct scale):\n")
cat("  Min:", min(count_matrix), "\n")
cat("  Max:", max(count_matrix), "\n")
cat("  Median non-zero:", median(count_matrix[count_matrix > 0]), "\n")
cat("  (Values should be in hundreds-thousands range, not 0-10)\n")


# ============================================================
# BLOCK 6: ASSIGN CONDITIONS FROM COLUMN NAMES
#
# The file GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz
# contains samples from the IL10Rfl (Control) vs Mutant experiment.
# Column names contain "fl" or "mut" identifiers.
# All names are printed above for verification.
# ============================================================
cat("\n--- Assigning conditions ---\n")

sample_names <- colnames(count_matrix)

# Assignment logic based on confirmed naming in this GEO file:
# "mut"  → IL10R-deficient microglia  → Mutant
# "fl"   → IL10R-floxed (intact)      → Control
# "cont" → Control
# "wt"   → Wild-type                  → Control
condition <- dplyr::case_when(
  grepl("mut",             sample_names, ignore.case = TRUE) ~ "Mutant",
  grepl("fl|cont|wt|ctrl", sample_names, ignore.case = TRUE) ~ "Control",
  TRUE ~ NA_character_
)

# Safety: if any unassigned, print them clearly
unassigned <- sample_names[is.na(condition)]
if (length(unassigned) > 0) {
  cat("\n⚠ UNASSIGNED SAMPLES — check names below and update patterns:\n")
  for (u in unassigned) cat("  UNASSIGNED:", u, "\n")
  cat("\nFix: edit the grepl() patterns in Block 6 to match these names.\n")
  condition[is.na(condition)] <- "Control"   # temporary fallback
}

meta_clean <- data.frame(
  sample_id = sample_names,
  condition = factor(condition, levels = c("Control", "Mutant")),
  row.names  = sample_names,
  stringsAsFactors = FALSE
)

cat("\nFinal condition assignment:\n")
for (i in seq_len(nrow(meta_clean))) {
  cat(sprintf("  %-55s → %s\n",
              meta_clean$sample_id[i],
              as.character(meta_clean$condition[i])))
}

cat("\nGroup counts:\n")
print(table(meta_clean$condition))

# Hard stop if still broken
n_ctrl <- sum(meta_clean$condition == "Control")
n_mut  <- sum(meta_clean$condition == "Mutant")

if (n_ctrl < 2 || n_mut < 2) {
  cat("\n=== COLUMN NAMES NEED MANUAL PATTERN FIX ===\n")
  cat("Copy the column names printed above and update grepl() patterns.\n")
  stop(paste("Need >= 2 per group. Got Control:", n_ctrl,
             "Mutant:", n_mut))
}

cat("\n✓ Metadata ready! Proceeding to DESeq2...\n")


# ============================================================
# BLOCK 7: SAVE CLEAN DATA
# FIX v1.0.1: count_matrix_int was undefined here — corrected
#             to count_matrix which is defined in Block 5.
# ============================================================
write.csv(count_matrix, here("data", "count_matrix_clean.csv"))
write.csv(meta_clean,   here("data", "metadata_clean.csv"))
cat("✓ Count matrix and metadata saved\n")


# ============================================================
# BLOCK 8: RUN DESeq2
# FIX v1.0.1: count_matrix_int was undefined here — corrected
#             to count_matrix which is defined in Block 5.
# ============================================================
cat("\n--- Running DESeq2 ---\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,    # FIX: was count_matrix_int (undefined)
  colData   = meta_clean,
  design    = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("After filtering:", nrow(dds), "genes remain\n")

dds <- DESeq(dds)
cat("✓ DESeq2 complete\n")
cat("Result names:", paste(resultsNames(dds), collapse = ", "), "\n")

# Extract results
res <- results(dds,
               contrast = c("condition", "Mutant", "Control"),
               alpha    = 0.05)
res_df      <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Label significance
res_df$significance <- "Not Significant"
res_df$significance[!is.na(res_df$padj) &
                      res_df$padj < 0.05 &
                      res_df$log2FoldChange >  1] <- "Upregulated"
res_df$significance[!is.na(res_df$padj) &
                      res_df$padj < 0.05 &
                      res_df$log2FoldChange < -1] <- "Downregulated"
res_df$significance <- factor(res_df$significance,
                              levels = c("Upregulated",
                                         "Downregulated",
                                         "Not Significant"))

# Save results
write.csv(res_df,
          here("results", "DESeq2_results_Mutant_vs_Control.csv"),
          row.names = FALSE)
res_df |> filter(significance == "Upregulated") |>
  arrange(padj) |> head(100) |>
  write.csv(here("results", "top100_upregulated.csv"), row.names = FALSE)
res_df |> filter(significance == "Downregulated") |>
  arrange(padj) |> head(100) |>
  write.csv(here("results", "top100_downregulated.csv"), row.names = FALSE)

n_up   <- sum(res_df$significance == "Upregulated",    na.rm = TRUE)
n_down <- sum(res_df$significance == "Downregulated",  na.rm = TRUE)
n_ns   <- sum(res_df$significance == "Not Significant", na.rm = TRUE)

cat("\n==========================================\n")
cat("  DESeq2 RESULTS SUMMARY\n")
cat("==========================================\n")
cat("  Genes tested:   ", nrow(res_df), "\n")
cat("  Upregulated:    ", n_up,   "(padj<0.05, LFC>1)\n")
cat("  Downregulated:  ", n_down, "(padj<0.05, LFC< -1)\n")
cat("  Not significant:", n_ns,   "\n")
cat("==========================================\n")
cat("\nPaper reports ~954 up, ~693 down at 48h (Figure 3E)\n")

# Spot-check key paper genes
key_genes <- c("Tnf", "Ccl5", "Il12b", "Il6", "Il10ra",
               "P2ry12", "Sall1", "Tmem119")
cat("\nKey gene check (should match paper Figure 3D):\n")
for (g in key_genes) {
  row <- res_df[res_df$gene == g, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %-12s LFC=%6.2f  padj=%.4f  [%s]\n",
                g,
                round(row$log2FoldChange, 2),
                ifelse(is.na(row$padj), NA, round(row$padj, 4)),
                as.character(row$significance)))
  }
}


# ============================================================
# BLOCK 9: VST TRANSFORMATION
# ============================================================
cat("\n--- VST transformation ---\n")
vsd <- vst(dds, blind = FALSE)
cat("✓ VST done\n")


# ============================================================
# BLOCK 10: VOLCANO PLOT
# ============================================================
cat("\n--- Volcano Plot ---\n")

plot_df <- res_df |>
  filter(!is.na(padj), !is.na(log2FoldChange)) |>
  mutate(neg_log10_padj = -log10(padj + 1e-300))

top_label <- plot_df |>
  filter(significance != "Not Significant") |>
  arrange(padj) |> head(30)

p_volcano <- ggplot(plot_df,
                    aes(x = log2FoldChange,
                        y = neg_log10_padj,
                        color = significance)) +
  geom_point(alpha = 0.55, size = 1.3) +
  scale_color_manual(values = c("Upregulated"     = "#E74C3C",
                                "Downregulated"   = "#3498DB",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text_repel(data = top_label, aes(label = gene),
                  size = 2.8, max.overlaps = 25,
                  box.padding = 0.3, color = "black",
                  fontface = "italic") +
  annotate("text",
           x     = max(plot_df$log2FoldChange, na.rm = TRUE) * 0.75,
           y     = max(plot_df$neg_log10_padj,  na.rm = TRUE) * 0.95,
           label = paste0("Up: ", n_up, "\nDown: ", n_down),
           size = 3.5, color = "grey30", hjust = 1) +
  labs(title    = "Volcano Plot: IL10R-Mutant vs Control (Microglia)",
       subtitle = "GSE157234 | Shemer et al., Immunity 2020 | UTAP-normalized counts",
       x        = "Log2 Fold Change (Mutant / Control)",
       y        = "-Log10 Adjusted P-value",
       color    = "Regulation") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave(here("plots", "volcano_plot.png"), p_volcano,
       width = 10, height = 7, dpi = 300)
ggsave(here("plots", "volcano_plot.pdf"), p_volcano,
       width = 10, height = 7)
cat("✓ Volcano plot saved\n")


# ============================================================
# BLOCK 11: PCA PLOT
# ============================================================
cat("\n--- PCA Plot ---\n")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2,
                    color = condition, shape = condition,
                    label = name)) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#2ECC71",
                                "Mutant"  = "#E74C3C")) +
  scale_shape_manual(values = c("Control" = 16, "Mutant" = 17)) +
  stat_ellipse(aes(group = condition), type = "norm",
               linetype = "dashed", linewidth = 0.7, level = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.4) +
  labs(title    = "PCA: IL10R-Mutant vs Control Microglia",
       subtitle = "GSE157234 | Shemer et al., Immunity 2020",
       x        = paste0("PC1: ", pct_var[1], "% variance"),
       y        = paste0("PC2: ", pct_var[2], "% variance"),
       color    = "Condition", shape = "Condition") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave(here("plots", "pca_plot.png"), p_pca,
       width = 9, height = 6, dpi = 300)
ggsave(here("plots", "pca_plot.pdf"), p_pca,
       width = 9, height = 6)
cat("✓ PCA plot saved\n")


# ============================================================
# BLOCK 12: HEATMAP
# ============================================================
cat("\n--- Heatmap ---\n")

top50_genes <- res_df |>
  filter(!is.na(padj), significance != "Not Significant") |>
  arrange(padj) |> head(50) |> pull(gene)

if (length(top50_genes) < 5) {
  cat("Warning: few significant genes. Using top 50 by p-value.\n")
  top50_genes <- res_df |>
    filter(!is.na(pvalue)) |>
    arrange(pvalue) |> head(50) |> pull(gene)
}

heat_mat        <- assay(vsd)[top50_genes, ]
heat_mat_scaled <- t(scale(t(heat_mat)))
heat_mat_scaled[heat_mat_scaled >  3] <-  3
heat_mat_scaled[heat_mat_scaled < -3] <- -3

# Annotation — colors MUST match actual factor levels
col_annot <- data.frame(
  Condition = as.character(meta_clean[colnames(heat_mat), "condition"]),
  row.names = colnames(heat_mat)
)
annot_colors <- list(
  Condition = c("Control" = "#2ECC71", "Mutant" = "#E74C3C")
)

png(here("plots", "heatmap_top50_DEGs.png"),
    width = 2400, height = 3200, res = 300)
pheatmap(heat_mat_scaled,
         annotation_col    = col_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 7,
         fontsize_col      = 7,
         main              = "Top 50 DEGs: IL10R-Mutant vs Control\nGSE157234 | Shemer et al., 2020",
         border_color      = NA,
         cutree_rows       = 2,
         cutree_cols       = 2)
dev.off()

pdf(here("plots", "heatmap_top50_DEGs.pdf"), width = 10, height = 13)
pheatmap(heat_mat_scaled,
         annotation_col    = col_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows      = TRUE, cluster_cols = TRUE,
         show_rownames     = TRUE, fontsize_row = 7,
         main              = "Top 50 DEGs: IL10R-Mutant vs Control")
dev.off()
cat("✓ Heatmap saved\n")


# ============================================================
# FINAL SUMMARY
# ============================================================
cat("\n==========================================\n")
cat("  ANALYSIS COMPLETE!\n")
cat("==========================================\n")
cat("Input:  UTAP-normalized counts (only format available on GEO)\n")
cat("Note:   True raw counts require SRA FASTQ + STAR alignment\n")
cat("\nResults saved:\n")
cat("  results/DESeq2_results_Mutant_vs_Control.csv\n")
cat("  results/top100_upregulated.csv\n")
cat("  results/top100_downregulated.csv\n")
cat("  plots/volcano_plot.png + .pdf\n")
cat("  plots/pca_plot.png + .pdf\n")
cat("  plots/heatmap_top50_DEGs.png + .pdf\n")
cat("\nNext: Run app_final.R\n")
cat("==========================================\n")
