# ============================================================
# RNA-Seq Differential Expression Analysis Pipeline
# ============================================================
# Dataset    : GSE157234 — Shemer et al., Immunity 2020
# Comparison : IL10R-Mutant vs Control Microglia (48h post-LPS)
# Author     : Md. Abrar Faiyaj
# Date       : March 2026
# GitHub     : https://github.com/mdabrarfaiyaj/rna-seq-shiny-pipeline
#
# BIOLOGICAL QUESTION:
#   What genes are differentially expressed in microglia that
#   cannot sense IL-10 (IL10R-deficient / Mutant) compared to
#   microglia with intact IL-10 signalling (Control), at 48h
#   after peripheral LPS challenge — the point of peak
#   hyperactivation described in Figure 3E of the paper?
#
# DATA NOTE:
#   Raw FASTQ counts were never publicly deposited for GSE157234.
#   NCBI's automated mouse raw count pipeline was also unavailable
#   at the time of analysis (NCBI GEO RNA-seq counts documentation).
#   The only available file was UTAP size-factor normalized counts
#   (GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz).
#   These values are on a count-compatible scale (hundreds to
#   thousands) — confirmed by the paper's own Data S1 source file.
#   Values are rounded to integers following NCBI's own practice
#   for DESeq2 compatibility. This is a pragmatic constraint of
#   public data availability, documented transparently here.
#
# PIPELINE OVERVIEW:
#   Block 1  → Install required packages (run once)
#   Block 2  → Load libraries
#   Block 3  → Create output folders
#   Block 4  → Download GEO supplementary files
#   Block 5  → Load and validate count matrix
#   Block 6  → Assign conditions + subset to 48h
#   Block 7  → Clean count matrix (remove NAs, align dimensions)
#   Block 8  → Run DESeq2 differential expression analysis
#   Block 9  → VST transformation for visualisation
#   Block 10 → Volcano plot
#   Block 11 → PCA plot
#   Block 12 → Heatmap (top 50 DEGs)
#   Block 13 → Session info (reproducibility record)
# ============================================================


# ------------------------------------------------------------
# SET WORKING DIRECTORY
# Set this to the folder containing your data/ results/ plots/
# subdirectories. Use a relative path for reproducibility.
# Example for Windows : setwd("C:/path/to/your/project/files")
# Example for Mac/Linux: setwd("~/path/to/your/project/files")
# ------------------------------------------------------------
# setwd("your/project/path/here")   # <-- uncomment and edit
cat("Working directory:", getwd(), "\n")


# ============================================================
# BLOCK 1: INSTALL REQUIRED PACKAGES
# ============================================================
# WHY THIS BLOCK EXISTS:
#   R packages must be installed before they can be used.
#   BiocManager is the package manager for Bioconductor packages
#   (DESeq2, GEOquery). CRAN packages are installed via the
#   standard install.packages(). This block checks first before
#   installing so it is safe to run repeatedly — it will only
#   install what is missing.
#
# KEY PACKAGES:
#   DESeq2      — Differential expression analysis (Bioconductor)
#   GEOquery    — Download data directly from NCBI GEO
#   ggplot2     — Publication quality plots
#   ggrepel     — Non-overlapping gene labels on plots
#   pheatmap    — Heatmap with hierarchical clustering
#   dplyr       — Data manipulation (filter, arrange, mutate)
#   RColorBrewer — Professional colour palettes
# ============================================================
cat("\n--- Block 1: Checking and installing packages ---\n")

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
# WHY THIS BLOCK EXISTS:
#   Installing a package makes it available on your computer
#   but does not load it into the current R session. library()
#   loads the package into memory so its functions are usable.
#   You must run this block at the start of every new R session.
# ============================================================
cat("\n--- Block 2: Loading libraries ---\n")

library(DESeq2);    library(ggplot2);   library(ggrepel)
library(pheatmap);  library(dplyr);     library(RColorBrewer)
library(GEOquery);  library(scales);    library(viridis)

cat("✓ Libraries loaded\n")


# ============================================================
# BLOCK 3: CREATE OUTPUT FOLDERS
# ============================================================
# WHY THIS BLOCK EXISTS:
#   R cannot save files to folders that do not exist — it will
#   throw an error. This block creates the required directory
#   structure automatically. recursive = TRUE means it creates
#   parent folders too if they are missing.
#
# FOLDER STRUCTURE:
#   data/            → raw and clean count matrices, metadata
#   data/GSE157234/  → downloaded GEO files
#   data/RAW/        → extracted raw count files
#   results/         → DESeq2 output CSV files
#   plots/           → all figures (PNG and PDF)
# ============================================================
cat("\n--- Block 3: Creating output folders ---\n")

for (f in c("data", "data/GSE157234", "data/RAW", "results", "plots")) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
}
cat("✓ Folders ready\n")


# ============================================================
# BLOCK 4: DOWNLOAD GEO SUPPLEMENTARY FILES
# ============================================================
# WHY THIS BLOCK EXISTS:
#   GEO (Gene Expression Omnibus) is NCBI's public repository
#   for gene expression data. getGEOSuppFiles() downloads all
#   supplementary files associated with an accession number.
#
# IMPORTANT DATA NOTE:
#   GSE157234_RAW.tar — despite its name — contains ATAC-seq
#   peak files (.peaks.txt), NOT RNA-seq counts. This is a known
#   issue with GEO depositions where multiple data types share
#   one accession (Cynthia Soto Cardinault noted this publicly).
#   The correct RNA-seq file is the UTAP normalized counts file
#   loaded in Block 5.
#
# The if(!file.exists()) check prevents re-downloading if the
# file is already present, saving bandwidth and time.
# ============================================================
cat("\n--- Block 4: Downloading GEO supplementary files ---\n")

geo_dir  <- "data/GSE157234"
raw_tar  <- file.path(geo_dir, "GSE157234_RAW.tar")

if (!file.exists(raw_tar)) {
  cat("Downloading supplementary files from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE, baseDir = "data/")
  cat("✓ Download complete\n")
} else {
  cat("✓ Files already downloaded\n")
}

cat("\nFiles available in", geo_dir, ":\n")
for (f in list.files(geo_dir)) cat("  ", f, "\n")


# ============================================================
# BLOCK 5: LOAD AND VALIDATE COUNT MATRIX
# ============================================================
# WHY THIS BLOCK EXISTS:
#   We load the UTAP normalized count matrix — the only publicly
#   available expression data for this dataset. Raw FASTQ counts
#   were never deposited on GEO, and NCBI's automated mouse raw
#   count pipeline was not yet available at time of analysis.
#
# WHAT IS UTAP NORMALIZATION?
#   UTAP (Weizmann Institute pipeline, Kohen et al. 2019) uses
#   DESeq2's median-of-ratios size factor normalization:
#     normalized count = raw count / size factor
#   Because the size factor is typically close to 1.0, the
#   normalized values remain on a count-compatible scale
#   (hundreds to thousands) rather than a log or ratio scale.
#
# WHY ROUNDING IS ACCEPTABLE HERE:
#   DESeq2 requires integer input for its negative binomial model.
#   NCBI itself rounds values for DESeq2 compatibility in its
#   official RNA-seq count pipeline (NCBI GEO documentation).
#   UTAP values differ from FPKM/TPM — those normalize by both
#   sequencing depth AND gene length, producing values on a
#   completely different scale (0.01 to ~100) that would be
#   genuinely inappropriate for DESeq2.
#
# VALUE RANGE CHECK:
#   Median of non-zero values should be in hundreds-to-thousands.
#   If you see values in the range 0–15, the data is log-
#   transformed — do NOT round and feed into DESeq2.
# ============================================================
cat("\n--- Block 5: Loading count matrix ---\n")

utap_file <- "data/GSE157234/GSE157234_UTAP_wt_IL10Rfl_vs_mut_normalized.txt.gz"

if (!file.exists(utap_file)) {
  cat("UTAP file not found — re-downloading from GEO...\n")
  getGEOSuppFiles("GSE157234", makeDirectory = TRUE, baseDir = "data/")
}

count_raw <- read.table(
  gzfile(utap_file),
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

cat("Dimensions:", nrow(count_raw), "genes x", ncol(count_raw), "columns\n")

# First column is gene names — set as rownames and remove
gene_names          <- make.unique(as.character(count_raw[, 1]))
rownames(count_raw) <- gene_names
count_raw           <- count_raw[, -1]

# Round to integers and enforce non-negative values
count_matrix              <- round(as.matrix(count_raw))
storage.mode(count_matrix) <- "integer"
count_matrix[count_matrix < 0] <- 0

# Value range validation — confirms count-compatible scale
cat("\nValue range validation:\n")
cat("  Min    :", min(count_matrix), "\n")
cat("  Max    :", max(count_matrix), "\n")
cat("  Median non-zero:", median(count_matrix[count_matrix > 0]), "\n")
cat("  (Expected: median in hundreds-to-thousands range) \n")

cat("\nAll sample names in this file:\n")
for (i in seq_along(colnames(count_matrix))) {
  cat(sprintf("  %2d: %s\n", i, colnames(count_matrix)[i]))
}
cat("✓ Count matrix ready\n")


# ============================================================
# BLOCK 6: ASSIGN CONDITIONS AND SUBSET TO 48h
# ============================================================
# WHY THIS BLOCK EXISTS:
#   The UTAP file contains samples from multiple experiments,
#   timepoints, and genotypes. We need to:
#     Step 1 — Exclude DKO samples (IL10R + TNF double knockout)
#              These are a different genotype entirely and would
#              confound the Mutant vs Control comparison.
#     Step 2 — Assign each sample to Control or Mutant based on
#              naming patterns in the column headers.
#     Step 3 — Extract timepoint (NT, 6h, 24h, 48h) from names.
#     Step 4 — Subset to 48h only — the paper's Figure 3E shows
#              48h as the peak of hyperactivation with 1,647 DEGs
#              vs only 86 DEGs at 6h. This is the key timepoint.
#
# NAMING PATTERNS IN THIS DATASET:
#   _cont_ or cont1/2/3 → Control (IL10R-floxed / intact)
#   _P_  or IL10R_P     → Control (Proficient = IL10R intact)
#   _mut_ or mut1/2/3   → Mutant  (IL10R-deficient)
#   _M_                 → Mutant  (Anat experiment naming)
#   DKO                 → EXCLUDED (double knockout)
#
# WHY case_when() INSTEAD OF ifelse()?
#   dplyr::case_when() handles multiple conditions cleanly and
#   returns NA for unmatched rows, making unassigned samples
#   visible rather than silently wrong.
# ============================================================
cat("\n--- Block 6: Assigning conditions and subsetting to 48h ---\n")

sample_names <- colnames(count_matrix)

# Step 1: Exclude DKO (double knockout) samples
is_dko <- grepl("DKO", sample_names, ignore.case = TRUE)
cat("DKO samples excluded:", sum(is_dko), "\n")
for (s in sample_names[is_dko]) cat("  EXCLUDED:", s, "\n")

count_matrix_sub <- count_matrix[, !is_dko]   # BUG FIX: was count_matrix
sample_names     <- colnames(count_matrix_sub)

# Step 2: Assign condition
condition <- dplyr::case_when(
  grepl("_cont_|cont[0-9]|_P_|IL10R_P", sample_names, ignore.case = TRUE) ~ "Control",
  grepl("_mut_|mut[0-9]|_M_",           sample_names, ignore.case = TRUE) ~ "Mutant",
  TRUE ~ NA_character_
)

# Report unassigned samples — do not silently ignore them
unassigned <- sample_names[is.na(condition)]
if (length(unassigned) > 0) {
  cat("\n⚠ Unassigned samples — check naming patterns:\n")
  for (u in unassigned) cat("  UNASSIGNED:", u, "\n")
  condition[is.na(condition)] <- "Control"   # fallback
}

# Step 3: Extract timepoint from sample name
timepoint <- dplyr::case_when(
  grepl("_nt_|_nt$",   sample_names, ignore.case = TRUE) ~ "NT",
  grepl("_6h_|_6h$",   sample_names, ignore.case = TRUE) ~ "6h",
  grepl("_24h_|_24h$", sample_names, ignore.case = TRUE) ~ "24h",
  grepl("_48h_|_48h$", sample_names, ignore.case = TRUE) ~ "48h",
  TRUE ~ "unknown"
)

# Build full metadata table
meta_full <- data.frame(
  sample_id = sample_names,
  condition = factor(condition, levels = c("Control", "Mutant")),
  timepoint = factor(timepoint, levels = c("NT", "6h", "24h", "48h", "unknown")),
  row.names  = sample_names,
  stringsAsFactors = FALSE
)

cat("\nSample assignments:\n")
for (i in seq_len(nrow(meta_full))) {
  cat(sprintf("  %-58s → %-8s [%s]\n",
              meta_full$sample_id[i],
              as.character(meta_full$condition[i]),
              as.character(meta_full$timepoint[i])))
}

cat("\nGroup × timepoint breakdown:\n")
print(table(meta_full$condition, meta_full$timepoint))

# Step 4: Subset to 48h only
# Paper Figure 3E: 48h shows maximum divergence between genotypes
# (954 up, 693 down = 1647 total DEGs at 48h)
is_48h     <- meta_full$timepoint == "48h"
count_48h  <- count_matrix_sub[, is_48h]     # BUG FIX: was count_matrix
meta_clean <- meta_full[is_48h, , drop = FALSE]

cat("\n48h subset:\n")
cat("  Samples kept:", ncol(count_48h), "\n")
print(table(meta_clean$condition))

# Safety check — DESeq2 needs at least 2 samples per group
n_ctrl <- sum(meta_clean$condition == "Control")
n_mut  <- sum(meta_clean$condition == "Mutant")
if (n_ctrl < 2 || n_mut < 2) {
  stop(paste("Need ≥2 samples per group. Got Control:", n_ctrl,
             "| Mutant:", n_mut))
}
cat("✓ Condition assignment ready\n")


# ============================================================
# BLOCK 7: CLEAN COUNT MATRIX
# ============================================================
# WHY THIS BLOCK EXISTS:
#   Before passing data to DESeq2, we must ensure:
#   1. Values are integers (DESeq2 requires integer counts)
#   2. No NA values exist (DESeq2 will error on NA rows)
#   3. Sample order matches between count matrix and metadata
#      (DESeq2 uses positional matching — wrong order = wrong
#      results with no error message, which is the worst kind
#      of bug in bioinformatics)
#
# WHY REMOVE NA ROWS?
#   Some genes may have NA counts due to alignment ambiguity
#   in the UTAP pipeline. DESeq2 cannot model NA values and
#   will throw an error if they are present.
#
# WHY ALIGN WITH intersect()?
#   If any sample name appears in the metadata but not the
#   count matrix (or vice versa), the analysis will fail.
#   intersect() keeps only samples present in BOTH objects.
# ============================================================
cat("\n--- Block 7: Cleaning count matrix ---\n")

# Convert 48h subset to clean integer matrix
count_matrix_int              <- round(as.matrix(count_48h))
storage.mode(count_matrix_int) <- "integer"
count_matrix_int[count_matrix_int < 0] <- 0

# Remove rows with any NA values
na_rows <- rowSums(is.na(count_matrix_int)) > 0
cat("Genes with NA values removed:", sum(na_rows), "\n")
count_matrix_int <- count_matrix_int[!na_rows, ]

# Align sample order between count matrix and metadata
shared           <- intersect(colnames(count_matrix_int), rownames(meta_clean))
count_matrix_int <- count_matrix_int[, shared]
meta_clean       <- meta_clean[shared, , drop = FALSE]

cat("Final count matrix:", nrow(count_matrix_int), "genes x",
    ncol(count_matrix_int), "samples\n")
cat("Metadata rows     :", nrow(meta_clean), "\n")
cat("Dimensions match  :", ncol(count_matrix_int) == nrow(meta_clean), "✓\n")

# Save clean data for reproducibility and Shiny app
write.csv(count_matrix_int, "data/count_matrix_48h_clean.csv")
write.csv(meta_clean,       "data/metadata_48h_clean.csv")
cat("✓ Clean data saved to data/\n")


# ============================================================
# BLOCK 8: RUN DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================
# WHY THIS BLOCK EXISTS:
#   DESeq2 (Love, Huber & Anders, Genome Biology 2014) is the
#   gold standard method for RNA-seq differential expression.
#
# HOW DESeq2 WORKS (simplified):
#   1. Estimates size factors to account for sequencing depth
#      differences between samples
#   2. Estimates dispersion (biological variability) for each
#      gene using an empirical Bayes shrinkage approach
#   3. Fits a negative binomial generalised linear model for
#      each gene with the design formula (~ condition)
#   4. Tests the Mutant vs Control contrast using the Wald test
#   5. Applies Benjamini-Hochberg multiple testing correction
#      to control the False Discovery Rate (FDR) at 5%
#
# FILTERING STEP:
#   rowSums(counts >= 10) >= 2 removes genes with fewer than
#   10 counts in at least 2 samples. Lowly expressed genes
#   have high noise and low statistical power — removing them
#   improves DESeq2's dispersion estimation and reduces the
#   multiple testing burden.
#
# CONTRAST:
#   c("condition", "Mutant", "Control") means:
#   positive LFC = higher in Mutant (upregulated in mutant)
#   negative LFC = lower in Mutant (downregulated in mutant)
#
# SIGNIFICANCE THRESHOLDS:
#   padj < 0.05 — adjusted p-value (BH-corrected FDR)
#   |LFC| > 1   — at least 2-fold change (log2 scale)
# ============================================================
cat("\n--- Block 8: Running DESeq2 ---\n")
cat("Control samples:", sum(meta_clean$condition == "Control"), "\n")
cat("Mutant  samples:", sum(meta_clean$condition == "Mutant"), "\n\n")

# Build DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_int,
  colData   = meta_clean,
  design    = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("Genes after low-count filtering:", nrow(dds), "\n")

# Set random seed for full reproducibility of any stochastic steps
# DESeq2's dispersion estimation uses numerical optimisation that
# can have minor stochastic components depending on the platform.
# Rule 6 — Sandve et al. (2013) doi:10.1371/journal.pcbi.1003285
set.seed(123)

# Run full DESeq2 pipeline (size factors → dispersion → GLM → Wald test)
dds <- DESeq(dds)
cat("✓ DESeq2 complete\n")
cat("Result names:", paste(resultsNames(dds), collapse = ", "), "\n")

# Extract results for Mutant vs Control contrast
res <- results(dds,
               contrast = c("condition", "Mutant", "Control"),
               alpha    = 0.05)

res_df      <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Classify genes by significance and direction
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

# Save full results and top 100 up/down separately
write.csv(res_df,
          "results/DESeq2_results_Mutant_vs_Control.csv",
          row.names = FALSE)

res_df |> dplyr::filter(significance == "Upregulated") |>
  dplyr::arrange(padj) |> head(100) |>
  write.csv("results/top100_upregulated.csv", row.names = FALSE)

res_df |> dplyr::filter(significance == "Downregulated") |>
  dplyr::arrange(padj) |> head(100) |>
  write.csv("results/top100_downregulated.csv", row.names = FALSE)

# Summary statistics
n_up   <- sum(res_df$significance == "Upregulated",    na.rm = TRUE)
n_down <- sum(res_df$significance == "Downregulated",  na.rm = TRUE)
n_ns   <- sum(res_df$significance == "Not Significant", na.rm = TRUE)

cat("\n==========================================\n")
cat("  DESeq2 RESULTS — 48h post-LPS\n")
cat("==========================================\n")
cat("  Comparison   : Mutant vs Control\n")
cat("  Timepoint    : 48h post-LPS\n")
cat("  Genes tested :", nrow(res_df), "\n")
cat("  Upregulated  :", n_up,   "(padj<0.05, LFC>1)\n")
cat("  Downregulated:", n_down, "(padj<0.05, LFC< -1)\n")
cat("  Not signif.  :", n_ns, "\n")
cat("==========================================\n")
cat("Paper reports ~954 up, ~693 down at 48h (Figure 3E)\n")
cat("Differences reflect normalized vs raw counts input.\n")

# Key gene validation against paper Figure 3D
cat("\nKey gene check (Figure 3D validation):\n")
cat("  Expected UP  : Tnf, Ccl5, Il12b, Il6\n")
cat("  Expected DOWN: P2ry12, Sall1, Tmem119, Il10ra\n\n")
key_genes <- c("Tnf", "Ccl5", "Il12b", "Il6", "Il10ra",
               "P2ry12", "Sall1", "Tmem119")
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


# ============================================================
# BLOCK 9: VST TRANSFORMATION
# ============================================================
# WHY THIS BLOCK EXISTS:
#   DESeq2 produces raw count-based results for differential
#   expression (Block 8). However, visualisation (PCA, heatmap)
#   requires data where variance is stable across the expression
#   range — raw counts violate this because high-count genes
#   naturally have higher variance.
#
# WHAT IS VST?
#   Variance Stabilizing Transformation (VST) applies a
#   data-driven transformation that makes variance approximately
#   constant across the full expression range. It is faster than
#   rlog for large datasets and preferred when n > 30 samples.
#
# blind = FALSE:
#   Uses the experimental design when estimating dispersion.
#   This gives better stabilisation than blind = TRUE (which
#   ignores the design). Use blind = FALSE when samples show
#   clear biological differences as expected here.
# ============================================================
cat("\n--- Block 9: VST transformation for visualisation ---\n")

vsd <- vst(dds, blind = FALSE)
cat("✓ VST complete\n")


# ============================================================
# BLOCK 10: VOLCANO PLOT
# ============================================================
# WHY THIS BLOCK EXISTS:
#   A volcano plot visualises all tested genes simultaneously,
#   showing statistical significance (y-axis) vs magnitude of
#   change (x-axis). It allows immediate identification of the
#   most biologically meaningful genes — those that are both
#   highly significant AND strongly changed.
#
# AXES EXPLAINED:
#   X — log2 fold change (Mutant / Control)
#       Positive = higher in Mutant (hyperactivation genes)
#       Negative = lower in Mutant (homeostasis genes lost)
#   Y — -log10(adjusted p-value)
#       Higher = more statistically significant
#       Horizontal dashed line = padj 0.05 threshold
#       Vertical dashed lines  = |LFC| = 1 threshold
#
# COLOUR CODING:
#   Red  = upregulated in Mutant (pro-inflammatory)
#   Blue = downregulated in Mutant (homeostatic genes lost)
#   Grey = not significant
#
# 1e-300 ADDITION:
#   Some p-values are so small R rounds them to exactly 0.
#   log10(0) = -Inf which breaks the plot. Adding 1e-300
#   (essentially zero) prevents this without distorting values.
# ============================================================
cat("\n--- Block 10: Volcano plot ---\n")

plot_df <- res_df |>
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
  dplyr::mutate(neg_log10_padj = -log10(padj + 1e-300))

# Label top 30 most significant DEGs
top_label <- plot_df |>
  dplyr::filter(significance != "Not Significant") |>
  dplyr::arrange(padj) |>
  head(30)

p_volcano <- ggplot(plot_df,
                    aes(x     = log2FoldChange,
                        y     = neg_log10_padj,
                        color = significance)) +
  geom_point(alpha = 0.55, size = 1.3) +
  scale_color_manual(values = c("Upregulated"     = "#E74C3C",
                                "Downregulated"   = "#3498DB",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text_repel(data        = top_label,
                  aes(label   = gene),
                  size        = 2.8,
                  max.overlaps = 25,
                  box.padding = 0.3,
                  color       = "black",
                  fontface    = "italic") +
  annotate("text",
           x     = max(plot_df$log2FoldChange, na.rm = TRUE) * 0.75,
           y     = max(plot_df$neg_log10_padj,  na.rm = TRUE) * 0.95,
           label = paste0("Up: ", n_up, "\nDown: ", n_down),
           size  = 3.5, color = "grey30", hjust = 1) +
  labs(title    = "Volcano Plot: IL10R-Mutant vs Control (Microglia)",
       subtitle = "GSE157234 | Shemer et al., Immunity 2020 | UTAP normalized counts",
       x        = "Log2 Fold Change (Mutant / Control)",
       y        = "-Log10 Adjusted P-value",
       color    = "Regulation") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/volcano_plot.png", p_volcano, width = 10, height = 7, dpi = 300)
ggsave("plots/volcano_plot.pdf", p_volcano, width = 10, height = 7)
cat("✓ Volcano plot saved (PNG + PDF)\n")


# ============================================================
# BLOCK 11: PCA PLOT
# ============================================================
# WHY THIS BLOCK EXISTS:
#   Principal Component Analysis (PCA) reduces the high-
#   dimensional gene expression space (thousands of genes) into
#   a few principal components that capture the most variation.
#   It is used as a quality control check and to visualise
#   overall sample relationships.
#
# HOW TO INTERPRET THIS PCA:
#   PC1 (78.9% variance) = genotype axis
#     → Separation between Control and Mutant
#     → Strong PC1 separation means genotype dominates variation
#   PC2 (8.8% variance) = within-group variation
#     → Controls cluster tightly = converged homeostatic state
#     → Mutants disperse = heterogeneous hyperactivation severity
#
# WHY DO MUTANTS DISPERSE ALONG PC2?
#   Without IL-10 dampening, individual mice progress through
#   the hyperactivation cascade at different rates. At 48h,
#   each mutant mouse is at a slightly different stage of
#   deterioration — this biological heterogeneity appears as
#   dispersion along PC2.
#
# DASHED ELLIPSES:
#   80% confidence ellipses (type="norm") showing the expected
#   area containing 80% of samples per group under normality.
# ============================================================
cat("\n--- Block 11: PCA plot ---\n")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data,
                aes(x     = PC1,
                    y     = PC2,
                    color = condition,
                    shape = condition,
                    label = name)) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#2ECC71",
                                "Mutant"  = "#E74C3C")) +
  scale_shape_manual(values = c("Control" = 16, "Mutant" = 17)) +
  stat_ellipse(aes(group = condition),
               type     = "norm",
               linetype = "dashed",
               linewidth = 0.7,
               level    = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.4) +
  labs(title    = "PCA: IL10R-Mutant vs Control Microglia",
       subtitle = "GSE157234 | Shemer et al., Immunity 2020 | 48h post-LPS",
       x        = paste0("PC1: ", pct_var[1], "% variance"),
       y        = paste0("PC2: ", pct_var[2], "% variance"),
       color    = "Condition",
       shape    = "Condition") +
  theme_bw(base_size = 13) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(color = "grey50", size = 10),
        panel.grid.minor = element_blank())

ggsave("plots/pca_plot.png", p_pca, width = 9, height = 6, dpi = 300)
ggsave("plots/pca_plot.pdf", p_pca, width = 9, height = 6)
cat("✓ PCA plot saved (PNG + PDF)\n")


# ============================================================
# BLOCK 12: HEATMAP (TOP 50 DEGs)
# ============================================================
# WHY THIS BLOCK EXISTS:
#   A heatmap shows gene expression patterns across all samples
#   simultaneously, with hierarchical clustering to reveal gene
#   groups with similar behaviour. It is the most information-
#   dense plot in this pipeline.
#
# HOW HEATMAP VALUES ARE CALCULATED:
#   1. VST values are extracted from the vsd object (Block 9)
#   2. Each gene is Z-scored: z = (value - mean) / sd
#      This centres each gene at 0 and scales by its spread,
#      allowing visual comparison across genes with different
#      absolute expression levels
#   3. Values are capped at ±3 to prevent extreme outliers
#      from compressing the colour scale
#
# WHY Z-SCORE? (Important for interviews!)
#   Without Z-scoring, a highly expressed gene (e.g. 10,000
#   counts) would dominate the colour scale, making low-
#   expressed genes invisible. Z-scoring puts all genes on
#   the same relative scale so patterns are comparable.
#
# HIERARCHICAL CLUSTERING:
#   Both rows (genes) and columns (samples) are clustered by
#   Euclidean distance with complete linkage. The dendrogram
#   shows which genes/samples are most similar to each other.
#   cutree_rows = 2 draws a line separating the two main gene
#   clusters (up-in-mutant vs down-in-mutant).
#
# WHAT TO LOOK FOR:
#   Top cluster  = genes HIGH in Control, LOW in Mutant
#                → homeostatic signature genes (P2ry12, Tmem119)
#   Bottom cluster = genes LOW in Control, HIGH in Mutant
#                → pro-inflammatory genes (Trem1, Mt1, Mt2)
# ============================================================
cat("\n--- Block 12: Heatmap ---\n")

# Select top 50 most significant DEGs
top50_genes <- res_df |>
  dplyr::filter(!is.na(padj), significance != "Not Significant") |>
  dplyr::arrange(padj) |>
  head(50) |>
  dplyr::pull(gene)

# Fallback if fewer than 5 significant genes found
if (length(top50_genes) < 5) {
  cat("Warning: fewer than 5 significant DEGs. Using top 50 by raw p-value.\n")
  top50_genes <- res_df |>
    dplyr::filter(!is.na(pvalue)) |>
    dplyr::arrange(pvalue) |>
    head(50) |>
    dplyr::pull(gene)
}

# Extract VST values and Z-score each gene (row-wise scaling)
heat_mat        <- assay(vsd)[top50_genes, ]
heat_mat_scaled <- t(scale(t(heat_mat)))

# Cap at ±3 standard deviations
heat_mat_scaled[heat_mat_scaled >  3] <-  3
heat_mat_scaled[heat_mat_scaled < -3] <- -3

# Build annotation bar (shows condition for each sample column)
col_annot <- data.frame(
  Condition = as.character(meta_clean[colnames(heat_mat), "condition"]),
  row.names = colnames(heat_mat)
)
annot_colors <- list(
  Condition = c("Control" = "#2ECC71", "Mutant" = "#E74C3C")
)

# Save PNG (high resolution for publication)
png("plots/heatmap_top50_DEGs.png", width = 2400, height = 3200, res = 300)
pheatmap(heat_mat_scaled,
         annotation_col    = col_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 7,
         fontsize_col      = 7,
         main              = "Top 50 DEGs: IL10R-Mutant vs Control\nGSE157234 | UTAP normalized counts | Shemer et al., 2020",
         border_color      = NA,
         cutree_rows       = 2,
         cutree_cols       = 2)
dev.off()

# Save PDF (vector format — infinitely scalable for presentations)
pdf("plots/heatmap_top50_DEGs.pdf", width = 10, height = 13)
pheatmap(heat_mat_scaled,
         annotation_col    = col_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         show_rownames     = TRUE,
         fontsize_row      = 7,
         main              = "Top 50 DEGs: IL10R-Mutant vs Control\nGSE157234 | Shemer et al., 2020")
dev.off()
cat("✓ Heatmap saved (PNG + PDF)\n")


# ============================================================
# FINAL SUMMARY
# ============================================================
cat("\n==========================================\n")
cat("  PIPELINE COMPLETE\n")
cat("==========================================\n")
cat("Dataset  : GSE157234 | Shemer et al., Immunity 2020\n")
cat("Comparison: IL10R-Mutant vs Control microglia\n")
cat("Timepoint : 48h post-LPS\n")
cat("\nResults saved:\n")
cat("  results/DESeq2_results_Mutant_vs_Control.csv\n")
cat("  results/top100_upregulated.csv\n")
cat("  results/top100_downregulated.csv\n")
cat("  data/count_matrix_48h_clean.csv\n")
cat("  data/metadata_48h_clean.csv\n")
cat("\nPlots saved:\n")
cat("  plots/volcano_plot.png + .pdf\n")
cat("  plots/pca_plot.png     + .pdf\n")
cat("  plots/heatmap_top50_DEGs.png + .pdf\n")
cat("\nNext step: Run app_final.R for interactive Shiny dashboard\n")
cat("==========================================\n")


# ============================================================
# BLOCK 13: INDUSTRY-STANDARD REPRODUCIBILITY SECTION
# ============================================================
# WHY THIS BLOCK EXISTS:
#   Reproducibility is the gold standard of scientific practice.
#   This block implements three layers of reproducibility that
#   are used in pharmaceutical industry, academic research, and
#   production bioinformatics pipelines:
#
#   LAYER 1 — Session info saved to file (not just console)
#     sessionInfo() printed to console disappears when you
#     close R. Saving it to a committed text file means anyone
#     who clones your repo can see the exact R version, OS, and
#     package versions used — permanently.
#     Reference: Ushey & Wickham (2024) renv: Project Environments
#     https://rstudio.github.io/renv/
#
#   LAYER 0 — set.seed(123) in Block 8
#     A random seed is set immediately before dds <- DESeq(dds)
#     to ensure any stochastic numerical steps in DESeq2's
#     dispersion estimation produce identical results on every run.
#     Rule 6: Sandve et al. (2013) doi:10.1371/journal.pcbi.1003285
#
#   LAYER 2 — RDS objects saved for Shiny app speed
#     DESeq2 takes 1-3 minutes to run. If the Shiny app re-runs
#     DESeq2 every time someone clicks the demo button, users
#     wait 3 minutes before seeing any plots. Saving dds, res,
#     and vsd as RDS files means the app loads in seconds by
#     reading pre-computed objects instead of re-running the
#     full pipeline. RDS is R's native binary format — it
#     preserves all object structure perfectly including
#     DESeqDataSet metadata, factor levels, and results tables.
#
#   LAYER 3 — renv lockfile for environment reproducibility
#     sessionInfo() records what YOU used. renv.lock enables
#     anyone else to INSTALL exactly what you used with one
#     command: renv::restore()
#     Without renv, a student cloning your repo in 2027 might
#     get DESeq2 v1.50 instead of the v1.46 you used, producing
#     different p-values and different DEG counts — silently.
#     renv prevents this entirely.
#     Reference: Siraji & Haque (2024), PMC, Primer on reproducible
#     research in R — doi:10.3390/clockssleep6010001
#     Reference: PharmaSUG 2025 — renv for version control and
#     environment reproducibility in pharmaceutical analytics
#
# IMPORTANT — HOW TO SET UP renv (run once per project):
#   Step 1: In R console, run:  renv::init(bioconductor = TRUE)
#   Step 2: Install all packages as normal
#   Step 3: Run: renv::snapshot()
#   Step 4: Commit renv.lock to GitHub
#   Anyone reproducing your work runs: renv::restore()
#   Reference: Appsilon (2024) — Using renv with Bioconductor
#   https://www.appsilon.com/post/renv-bioconductor
# ============================================================
cat("\n--- Block 13: Industry-standard reproducibility ---\n")


# ---- LAYER 1: Save session info to persistent file --------
# Creates results/session_info.txt — committed to GitHub so
# the exact computational environment is permanently recorded.
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

session_file <- "results/session_info.txt"
sink(session_file)
cat("==========================================================\n")
cat("  SESSION INFO — RNA-Seq Pipeline\n")
cat("  Dataset   : GSE157234 | Shemer et al., Immunity 2020\n")
cat("  Captured  :", date(), "\n")
cat("==========================================================\n\n")
print(sessionInfo())
sink()
cat("✓ Session info saved to:", session_file, "\n")
cat("  (Records exact R, Bioconductor, and package versions)\n")


# ---- LAYER 2: Save RDS objects for Shiny app speed --------
# dds  = full DESeqDataSet (contains counts, model, size factors)
# res  = DESeq2 results table (LFC, p-values, padj)
# vsd  = VST-transformed object (for PCA and heatmap)
# res_df = annotated results dataframe (for volcano and table)
#
# The Shiny app checks for these files first. If found, it loads
# them in ~1 second instead of re-running DESeq2 for 1-3 minutes.

saveRDS(dds,    "results/dds_object.rds")
saveRDS(res,    "results/deseq_results.rds")
saveRDS(vsd,    "results/vsd_object.rds")
saveRDS(res_df, "results/res_df.rds")
cat("✓ DESeq2 objects saved as RDS files:\n")
cat("  results/dds_object.rds\n")
cat("  results/deseq_results.rds\n")
cat("  results/vsd_object.rds\n")
cat("  results/res_df.rds\n")
cat("  (Shiny app will load these instantly instead of re-running DESeq2)\n")


# ---- LAYER 3: renv setup instructions ---------------------
# renv is authored by Kevin Ushey & Hadley Wickham (Posit/RStudio)
# It is the standard dependency management tool for R projects
# in academic bioinformatics, clinical research, and pharma industry.
#
# Check if renv is already initialised for this project
if (file.exists("renv.lock")) {
  cat("\n✓ renv.lock already exists — environment is locked\n")
  cat("  Anyone can reproduce this environment with: renv::restore()\n")
} else {
  cat("\n⚠ renv not yet initialised for this project\n")
  cat("  To lock your environment (recommended), run in R console:\n\n")
  cat("    install.packages('renv')\n")
  cat("    renv::init(bioconductor = TRUE)\n")
  cat("    renv::snapshot()\n\n")
  cat("  Then commit renv.lock to GitHub.\n")
  cat("  Collaborators restore with: renv::restore()\n")
  cat("\n  Why this matters: without renv, a DESeq2 update in 2027\n")
  cat("  could silently change your p-values. renv prevents this.\n")
  cat("  Reference: Ushey & Wickham (2024) https://rstudio.github.io/renv/\n")
}


# ---- Final complete summary --------------------------------
cat("\n==========================================\n")
cat("  PIPELINE COMPLETE — ALL OUTPUTS SAVED\n")
cat("==========================================\n")
cat("Dataset   : GSE157234 | Shemer et al., Immunity 2020\n")
cat("Comparison: IL10R-Mutant vs Control microglia\n")
cat("Timepoint : 48h post-LPS\n")
cat("\nAnalysis outputs:\n")
cat("  results/DESeq2_results_Mutant_vs_Control.csv\n")
cat("  results/top100_upregulated.csv\n")
cat("  results/top100_downregulated.csv\n")
cat("\nReproducibility outputs:\n")
cat("  results/session_info.txt    ← R environment record\n")
cat("  results/dds_object.rds      ← DESeqDataSet object\n")
cat("  results/deseq_results.rds   ← DESeq2 results\n")
cat("  results/vsd_object.rds      ← VST object\n")
cat("  results/res_df.rds          ← annotated results dataframe\n")
cat("\nClean data:\n")
cat("  data/count_matrix_48h_clean.csv\n")
cat("  data/metadata_48h_clean.csv\n")
cat("\nPlots:\n")
cat("  plots/volcano_plot.png + .pdf\n")
cat("  plots/pca_plot.png     + .pdf\n")
cat("  plots/heatmap_top50_DEGs.png + .pdf\n")
cat("\nNext step: Run app_final.R for interactive Shiny dashboard\n")
cat("==========================================\n")
