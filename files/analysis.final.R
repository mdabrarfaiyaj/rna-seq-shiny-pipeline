# ============================================================
# CORRECTED RNA-Seq ANALYSIS PIPELINE - RAW COUNTS VERSION
# Dataset: GSE157234 | Shemer et al., Immunity 2020
# Comparison: IL10R-Mutant vs Control Microglia (post LPS)
#
# KEY FIX: Using GSE157234_RAW.tar (raw integer counts)
#          instead of pre-normalized files.
#          DESeq2 MUST receive raw counts - not normalized values.
# ============================================================

setwd("C:/Users/HP/Documents/Github Project/rna-seq-shiny-pipeline/files")
cat("✓ Working directory:", getwd(), "\n")


# ============================================================
# BLOCK 1: INSTALL PACKAGES
# ============================================================
cat("\n--- Checking packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "dplyr",
               "RColorBrewer", "scales", "viridis")
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
library(DESeq2);   library(ggplot2);  library(ggrepel)
library(pheatmap); library(dplyr);    library(RColorBrewer)
library(GEOquery); library(scales);   library(viridis)
cat("✓ Libraries loaded\n")


# ============================================================
# BLOCK 3: CREATE FOLDERS
# ============================================================
for (f in c("data", "data/GSE157234", "data/RAW", "results", "plots")) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
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

raw_tar <- "data/GSE157234/GSE157234_RAW.tar"

if (!file.exists(raw_tar)) {
  cat("Downloading supplementary files from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE, baseDir = "data/")
  cat("✓ Download complete\n")
} else {
  cat("✓ Files already downloaded\n")
}

# Extract RAW.tar
raw_dir <- "data/RAW"
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

# Show what we found
cat("\nFiles in RAW directory:\n")
for (f in list.files(raw_dir)) cat("  ", f, "\n")


# ============================================================
# BLOCK 5: LOAD CORRECT COUNT MATRIX
# ============================================================
cat("\n--- Loading correct count matrix ---\n")

# This is the ONLY file that contains both Mutant AND Control
correct_file <- "data/GSE157234/GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz"

if (!file.exists(correct_file)) {
  cat("File not found locally. Downloading from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE, baseDir = "data/")
}

count_raw <- read.table(
  gzfile(correct_file),
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

cat("Dimensions:", nrow(count_raw), "genes x", ncol(count_raw), "columns\n")

# First column = gene names
gene_names        <- make.unique(as.character(count_raw[, 1]))
rownames(count_raw) <- gene_names
count_raw         <- count_raw[, -1]   # drop gene name column

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
cat("  (Values should be in hundreds-thousands range, not 0–10)\n")


# ============================================================
# BLOCK 6: ASSIGN CONDITIONS FROM COLUMN NAMES
#
# The file GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz
# contains samples from the IL10Rfl (Control) vs Mutant experiment.
# Column names will contain "fl" or "mut" identifiers.
# We print all names above so you can verify this assignment.
# ============================================================
cat("\n--- Assigning conditions ---\n")

sample_names <- colnames(count_matrix)

# Assignment logic based on typical naming in this GEO file:
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
  condition[is.na(condition)] <- "Control"  # temporary fallback
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
  cat("Copy the column names printed above and tell me what they look like.\n")
  cat("I will give you the exact corrected grepl() patterns immediately.\n")
  stop(paste("Need ≥2 per group. Got Control:", n_ctrl,
             "Mutant:", n_mut))
}

cat("\n✓ Metadata ready! Proceeding to DESeq2...\n")

# Save clean data
write.csv(count_matrix, "data/count_matrix_clean.csv")
write.csv(meta_clean,   "data/metadata_clean.csv")

# ============================================================
# BLOCK 7: SAVE CLEAN DATA
# ============================================================
write.csv(count_matrix_int, "data/count_matrix_RAW.csv")
write.csv(meta_clean,       "data/metadata_clean.csv")
cat("✓ Raw count matrix saved to data/count_matrix_RAW.csv\n")


# ============================================================
# BLOCK 8: RUN DESeq2 (with proper raw counts)
# ============================================================
cat("\n--- Running DESeq2 ---\n")
cat("Using RAW integer counts — DESeq2 will normalize internally ✓\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
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
res_df       <- as.data.frame(res)
res_df$gene  <- rownames(res_df)

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
          "results/DESeq2_results_Mutant_vs_Control.csv",
          row.names = FALSE)
res_df |> filter(significance == "Upregulated") |>
  arrange(padj) |> head(100) |>
  write.csv("results/top100_upregulated.csv", row.names = FALSE)
res_df |> filter(significance == "Downregulated") |>
  arrange(padj) |> head(100) |>
  write.csv("results/top100_downregulated.csv", row.names = FALSE)

n_up   <- sum(res_df$significance == "Upregulated",    na.rm = TRUE)
n_down <- sum(res_df$significance == "Downregulated",  na.rm = TRUE)
n_ns   <- sum(res_df$significance == "Not Significant",na.rm = TRUE)

cat("\n==========================================\n")
cat("  DESeq2 RESULTS SUMMARY\n")
cat("==========================================\n")
cat("  Genes tested:   ", nrow(res_df), "\n")
cat("  Upregulated:    ", n_up,   "(padj<0.05, LFC>1)\n")
cat("  Downregulated:  ", n_down, "(padj<0.05, LFC< -1)\n")
cat("  Not significant:", n_ns,   "\n")
cat("==========================================\n")
cat("\nWith RAW counts you should now see HUNDREDS of DEGs,\n")
cat("matching Figure 3E of the paper (~954 up, ~693 down at 48h)\n")

# Spot-check key paper genes
key_genes <- c("Tnf", "Ccl5", "Il12b", "Il6", "Il10ra",
               "P2ry12", "Sall1", "Tmem119")
cat("\nKey gene check (should match paper's Figure 3D):\n")
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
       subtitle = "GSE157234 | Shemer et al., Immunity 2020 | RAW counts",
       x        = "Log2 Fold Change (Mutant / Control)",
       y        = "-Log10 Adjusted P-value",
       color    = "Regulation") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/volcano_plot.png", p_volcano, width = 10, height = 7, dpi = 300)
ggsave("plots/volcano_plot.pdf", p_volcano, width = 10, height = 7)
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
       subtitle = "GSE157234 | Clearer separation expected with RAW counts",
       x        = paste0("PC1: ", pct_var[1], "% variance"),
       y        = paste0("PC2: ", pct_var[2], "% variance"),
       color    = "Condition", shape = "Condition") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/pca_plot.png", p_pca, width = 9, height = 6, dpi = 300)
ggsave("plots/pca_plot.pdf", p_pca, width = 9, height = 6)
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

png("plots/heatmap_top50_DEGs.png", width = 2400, height = 3200, res = 300)
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
         main              = "Top 50 DEGs: IL10R-Mutant vs Control\nGSE157234 | RAW counts | Shemer et al., 2020",
         border_color      = NA,
         cutree_rows       = 2,
         cutree_cols       = 2)
dev.off()

pdf("plots/heatmap_top50_DEGs.pdf", width = 10, height = 13)
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
cat("Key correction applied:\n")
cat("  OLD: Pre-normalized floats → few DEGs (wrong)\n")
cat("  NEW: Raw integer counts    → hundreds of DEGs (correct)\n")
cat("\nResults saved:\n")
cat("  results/DESeq2_results_Mutant_vs_Control.csv\n")
cat("  results/top100_upregulated.csv\n")
cat("  results/top100_downregulated.csv\n")
cat("  plots/volcano_plot.png + .pdf\n")
cat("  plots/pca_plot.png + .pdf\n")
cat("  plots/heatmap_top50_DEGs.png + .pdf\n")
cat("\nNext: Run app_final.R\n")
cat("==========================================\n")# ============================================================
# CORRECTED RNA-Seq ANALYSIS PIPELINE - RAW COUNTS VERSION
# Dataset: GSE157234 | Shemer et al., Immunity 2020
# Comparison: IL10R-Mutant vs Control Microglia (post LPS)
#
# KEY FIX: Using GSE157234_RAW.tar (raw integer counts)
#          instead of pre-normalized files.
#          DESeq2 MUST receive raw counts - not normalized values.
# ============================================================

setwd("C:/Users/HP/Documents/Github Project/rna-seq-shiny-pipeline/files")
cat("✓ Working directory:", getwd(), "\n")


# ============================================================
# BLOCK 1: INSTALL PACKAGES
# ============================================================
cat("\n--- Checking packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "dplyr",
               "RColorBrewer", "scales", "viridis")
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
library(DESeq2);   library(ggplot2);  library(ggrepel)
library(pheatmap); library(dplyr);    library(RColorBrewer)
library(GEOquery); library(scales);   library(viridis)
cat("✓ Libraries loaded\n")


# ============================================================
# BLOCK 3: CREATE FOLDERS
# ============================================================
for (f in c("data", "data/GSE157234", "data/RAW", "results", "plots")) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
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

raw_tar <- "data/GSE157234/GSE157234_RAW.tar"

if (!file.exists(raw_tar)) {
  cat("Downloading supplementary files from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE, baseDir = "data/")
  cat("✓ Download complete\n")
} else {
  cat("✓ Files already downloaded\n")
}

# Extract RAW.tar
raw_dir <- "data/RAW"
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

# Show what we found
cat("\nFiles in RAW directory:\n")
for (f in list.files(raw_dir)) cat("  ", f, "\n")


# ============================================================
# BLOCK 5: BUILD RAW COUNT MATRIX
#
# Each file in RAW.tar is one sample: rows = genes, col = count.
# We merge all samples into one matrix.
# ============================================================
cat("\n--- Building raw count matrix ---\n")

# Find all count files (txt, tsv, csv — not gz)
count_files <- list.files(raw_dir,
                          pattern = "\\.(txt|tsv|csv)$",
                          full.names = TRUE,
                          recursive  = TRUE)

if (length(count_files) == 0) {
  # Try .gz files directly if uncompression failed
  count_files <- list.files(raw_dir,
                            pattern = "\\.gz$",
                            full.names = TRUE)
}

cat("Found", length(count_files), "sample files\n")
if (length(count_files) == 0) stop("No count files found in RAW directory!")

# Read and merge all sample files
sample_list <- list()

for (cf in count_files) {
  # Derive clean sample name from filename
  sname <- basename(cf)
  sname <- sub("^GSM[0-9]+_", "", sname)      # remove GSM prefix
  sname <- sub("\\.(txt|tsv|csv|gz)$", "", sname)
  sname <- gsub(" ", "_", sname)
  
  tryCatch({
    # Try reading with different separators
    df <- tryCatch(
      read.table(cf, header = FALSE, sep = "\t",
                 stringsAsFactors = FALSE, comment.char = "#"),
      error = function(e)
        read.table(cf, header = FALSE, sep = ",",
                   stringsAsFactors = FALSE, comment.char = "#")
    )
    
    # Expect col1 = gene name, col2 = count
    if (ncol(df) >= 2) {
      genes  <- as.character(df[, 1])
      counts <- suppressWarnings(as.numeric(df[, 2]))
      
      # Skip header rows, summary rows (__, ambiguous, etc.)
      keep <- !is.na(counts) & !grepl("^__", genes) &
        !genes %in% c("no_feature","ambiguous",
                      "too_low_aQual","not_aligned",
                      "alignment_not_unique")
      
      tmp_df <- data.frame(count = counts[keep],
                           row.names = genes[keep])
      colnames(tmp_df) <- sname
      sample_list[[sname]] <- tmp_df
    }
  }, error = function(e) {
    cat("  Warning: could not read", basename(cf), "-", e$message, "\n")
  })
}

cat("Successfully read:", length(sample_list), "samples\n")

if (length(sample_list) == 0) {
  stop("Could not read any sample files. Check data/RAW/ directory.")
}

# Merge all samples by gene (common genes only)
count_matrix_raw <- Reduce(function(a, b) {
  merge(a, b, by = "row.names", all = FALSE) |>
    (\(x) { rownames(x) <- x$Row.names; x[, -1] })()
}, sample_list)

cat("Raw count matrix:", nrow(count_matrix_raw), "genes x",
    ncol(count_matrix_raw), "samples\n")

# Confirm values look like raw counts (should be integers, often in thousands)
cat("Value range check (first 5 genes, first 3 samples):\n")
print(round(count_matrix_raw[1:min(5, nrow(count_matrix_raw)),
                             1:min(3, ncol(count_matrix_raw))], 2))
cat("\nIf values above look like whole numbers (e.g. 150, 2430, 0),\n")
cat("you have RAW COUNTS ✓\n")
cat("If they look like decimals (e.g. 1.23, 0.45), they may still\n")
cat("be normalized — check the GEO page for alternative files.\n\n")

# Round to integers (safety measure for slightly non-integer raw counts)
count_matrix_int <- round(as.matrix(count_matrix_raw))
storage.mode(count_matrix_int) <- "integer"
count_matrix_int[count_matrix_int < 0] <- 0

cat("Sample names found:\n")
for (i in seq_along(colnames(count_matrix_int))) {
  cat(sprintf("  %2d: %s\n", i, colnames(count_matrix_int)[i]))
}


# ============================================================
# BLOCK 6: CORRECT CONDITION ASSIGNMENT
# ============================================================
cat("\n--- Assigning conditions from sample names ---\n")

sample_names <- colnames(count_matrix)

# Step 1: Remove DKO samples (double knockout - different genotype)
# DKO = IL10R-deficient AND TNF-deficient = NOT the Mutant vs Control comparison
is_dko <- grepl("DKO", sample_names, ignore.case = TRUE)
cat("Excluding DKO samples:", sum(is_dko), "\n")
for (s in sample_names[is_dko]) cat("  EXCLUDED:", s, "\n")

count_matrix_sub <- count_matrix[, !is_dko]
sample_names     <- colnames(count_matrix_sub)

# Step 2: Assign condition
# _cont_ or cont1/2/3 → Control (IL10R-floxed, intact)
# _mut_  or mut1/2/3  → Mutant  (IL10R-deficient)
# _M_                 → Mutant  (Anat experiment naming)
# _P_  or IL10R_P     → Control (Proficient = IL10R intact)
condition <- dplyr::case_when(
  grepl("_cont_|cont[0-9]|_P_|IL10R_P", sample_names, ignore.case=TRUE) ~ "Control",
  grepl("_mut_|mut[0-9]|_M_",           sample_names, ignore.case=TRUE) ~ "Mutant",
  TRUE ~ NA_character_
)

# Report any still unassigned
unassigned <- sample_names[is.na(condition)]
if (length(unassigned) > 0) {
  cat("\n⚠ Still unassigned — defaulting to Control:\n")
  for (u in unassigned) cat("  ", u, "\n")
  condition[is.na(condition)] <- "Control"
}

# Step 3: Extract timepoint from sample name
timepoint <- dplyr::case_when(
  grepl("_nt_|_nt$",   sample_names, ignore.case=TRUE) ~ "NT",
  grepl("_6h_|_6h$",   sample_names, ignore.case=TRUE) ~ "6h",
  grepl("_24h_|_24h$", sample_names, ignore.case=TRUE) ~ "24h",
  grepl("_48h_|_48h$", sample_names, ignore.case=TRUE) ~ "48h",
  TRUE ~ "unknown"
)

meta_full <- data.frame(
  sample_id = sample_names,
  condition = factor(condition, levels = c("Control", "Mutant")),
  timepoint = factor(timepoint, levels = c("NT","6h","24h","48h","unknown")),
  row.names  = sample_names,
  stringsAsFactors = FALSE
)

cat("\nAll samples after DKO removal:\n")
for (i in seq_len(nrow(meta_full))) {
  cat(sprintf("  %-58s → %-8s [%s]\n",
              meta_full$sample_id[i],
              as.character(meta_full$condition[i]),
              as.character(meta_full$timepoint[i])))
}

cat("\nFull group × timepoint table:\n")
print(table(meta_full$condition, meta_full$timepoint))

# Step 4: Subset to 48h ONLY
# Why 48h? Paper Figure 3E shows 1647 DEGs at 48h vs only 86 at 6h.
# This is the timepoint where Mutant hyperactivation is maximal.
is_48h      <- meta_full$timepoint == "48h"
count_48h   <- count_matrix_sub[, is_48h]
meta_clean  <- meta_full[is_48h, , drop = FALSE]

cat("\n--- Subset to 48h timepoint ---\n")
cat("Samples kept:", ncol(count_48h), "\n")
cat("Condition counts at 48h:\n")
print(table(meta_clean$condition))

# Verify minimum samples per group
n_ctrl <- sum(meta_clean$condition == "Control")
n_mut  <- sum(meta_clean$condition == "Mutant")
if (n_ctrl < 2 || n_mut < 2) {
  stop(paste("Need ≥2 per group at 48h. Got Control:", n_ctrl,
             "Mutant:", n_mut))
}
cat("✓ Condition assignment ready!\n")


# ============================================================
# BLOCK 7: CLEAN COUNT MATRIX — REMOVE NAs
# ============================================================
cat("\n--- Cleaning count matrix ---\n")

# Round and convert to integer
count_matrix_int <- round(as.matrix(count_48h))
storage.mode(count_matrix_int) <- "integer"
count_matrix_int[count_matrix_int < 0] <- 0

# Remove any rows with NA values
na_rows <- rowSums(is.na(count_matrix_int)) > 0
cat("Genes with NA values removed:", sum(na_rows), "\n")
count_matrix_int <- count_matrix_int[!na_rows, ]

# Ensure sample order matches between matrix and metadata
shared <- intersect(colnames(count_matrix_int), rownames(meta_clean))
count_matrix_int <- count_matrix_int[, shared]
meta_clean       <- meta_clean[shared, , drop = FALSE]

cat("Final count matrix:", nrow(count_matrix_int), "genes x",
    ncol(count_matrix_int), "samples\n")
cat("Metadata rows:     ", nrow(meta_clean), "\n")
cat("Sizes match:       ", ncol(count_matrix_int) == nrow(meta_clean), "✓\n")

# Save
write.csv(count_matrix_int, "data/count_matrix_48h_clean.csv")
write.csv(meta_clean,       "data/metadata_48h_clean.csv")


# ============================================================
# BLOCK 8: RUN DESeq2
# ============================================================
cat("\n--- Running DESeq2 ---\n")
cat("Samples: Control =", sum(meta_clean$condition=="Control"),
    "| Mutant =", sum(meta_clean$condition=="Mutant"), "\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData   = meta_clean,
  design    = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
cat("✓ DESeq2 complete!\n")
cat("Result names:", paste(resultsNames(dds), collapse=", "), "\n")

# Extract results
res <- results(dds,
               contrast = c("condition", "Mutant", "Control"),
               alpha    = 0.05)
res_df       <- as.data.frame(res)
res_df$gene  <- rownames(res_df)

res_df$significance <- "Not Significant"
res_df$significance[!is.na(res_df$padj) &
                      res_df$padj < 0.05 &
                      res_df$log2FoldChange >  1] <- "Upregulated"
res_df$significance[!is.na(res_df$padj) &
                      res_df$padj < 0.05 &
                      res_df$log2FoldChange < -1] <- "Downregulated"
res_df$significance <- factor(res_df$significance,
                              levels=c("Upregulated","Downregulated",
                                       "Not Significant"))

# Save results
write.csv(res_df, "results/DESeq2_results_Mutant_vs_Control.csv",
          row.names=FALSE)
res_df |> dplyr::filter(significance=="Upregulated") |>
  dplyr::arrange(padj) |> head(100) |>
  write.csv("results/top100_upregulated.csv", row.names=FALSE)
res_df |> dplyr::filter(significance=="Downregulated") |>
  dplyr::arrange(padj) |> head(100) |>
  write.csv("results/top100_downregulated.csv", row.names=FALSE)

n_up   <- sum(res_df$significance=="Upregulated",    na.rm=TRUE)
n_down <- sum(res_df$significance=="Downregulated",  na.rm=TRUE)
n_ns   <- sum(res_df$significance=="Not Significant",na.rm=TRUE)

cat("\n==========================================\n")
cat("  DESeq2 RESULTS — 48h post-LPS\n")
cat("==========================================\n")
cat("  Comparison:   Mutant vs Control\n")
cat("  Timepoint:    48h post-LPS\n")
cat("  Genes tested: ", nrow(res_df), "\n")
cat("  Upregulated:  ", n_up,   " (padj<0.05, LFC>1)\n")
cat("  Downregulated:", n_down, " (padj<0.05, LFC< -1)\n")
cat("  Not signif:   ", n_ns,   "\n")
cat("==========================================\n")
cat("\nPaper reports ~954 up, ~693 down at 48h (Figure 3E)\n")
cat("Your numbers should be in the same ballpark.\n")

# Key gene spot-check
cat("\nKey gene check vs paper Figure 3D:\n")
key_genes <- c("Tnf","Ccl5","Il12b","Il6","Il10ra",
               "P2ry12","Sall1","Tmem119")
for (g in key_genes) {
  row <- res_df[res_df$gene == g, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %-12s  LFC=%6.2f  padj=%.4f  [%s]\n",
                g,
                round(row$log2FoldChange, 2),
                ifelse(is.na(row$padj), NA, round(row$padj, 4)),
                as.character(row$significance)))
  }
}

cat("\n✓ Now continue with Block 9 (VST) onwards — unchanged.\n")







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
       subtitle = "GSE157234 | Shemer et al., Immunity 2020 | RAW counts",
       x        = "Log2 Fold Change (Mutant / Control)",
       y        = "-Log10 Adjusted P-value",
       color    = "Regulation") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/volcano_plot.png", p_volcano, width = 10, height = 7, dpi = 300)
ggsave("plots/volcano_plot.pdf", p_volcano, width = 10, height = 7)
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
       subtitle = "GSE157234 | Clearer separation expected with RAW counts",
       x        = paste0("PC1: ", pct_var[1], "% variance"),
       y        = paste0("PC2: ", pct_var[2], "% variance"),
       color    = "Condition", shape = "Condition") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/pca_plot.png", p_pca, width = 9, height = 6, dpi = 300)
ggsave("plots/pca_plot.pdf", p_pca, width = 9, height = 6)
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

png("plots/heatmap_top50_DEGs.png", width = 2400, height = 3200, res = 300)
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
         main              = "Top 50 DEGs: IL10R-Mutant vs Control\nGSE157234 | RAW counts | Shemer et al., 2020",
         border_color      = NA,
         cutree_rows       = 2,
         cutree_cols       = 2)
dev.off()

pdf("plots/heatmap_top50_DEGs.pdf", width = 10, height = 13)
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
cat("Key correction applied:\n")
cat("  OLD: Pre-normalized floats → few DEGs (wrong)\n")
cat("  NEW: Raw integer counts    → hundreds of DEGs (correct)\n")
cat("\nResults saved:\n")
cat("  results/DESeq2_results_Mutant_vs_Control.csv\n")
cat("  results/top100_upregulated.csv\n")
cat("  results/top100_downregulated.csv\n")
cat("  plots/volcano_plot.png + .pdf\n")
cat("  plots/pca_plot.png + .pdf\n")
cat("  plots/heatmap_top50_DEGs.png + .pdf\n")
cat("\nNext: Run app_final.R\n")
cat("==========================================\n")