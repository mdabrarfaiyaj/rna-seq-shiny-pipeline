# Contributing to RNA-seq Shiny Pipeline

Thank you for your interest in contributing to this project! This is an educational project that replicates and visualizes RNA-seq analysis from published datasets (e.g., Shemer et al. 2020, GSE157234) using an interactive Shiny dashboard.

This document outlines how to get involved, whether you're reporting bugs, suggesting improvements, or contributing code.

## Getting Started

### Development Setup

To set up the project locally:

1. **Clone the repository**
   ```bash
   git clone https://github.com/mdabrarfaiyaj/rna-seq-shiny-pipeline.git
   cd rna-seq-shiny-pipeline
   ```

2. **Install required R packages**
   ```r
   if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
   if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
   if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
   
   # See DESCRIPTION or README.md for full dependencies
   ```

3. **Run the Shiny app locally**
   ```r
   shiny::runApp("app.R")
   ```

4. **Test your changes**
   - Verify the dashboard loads without errors
   - Test all interactive elements (filters, plots, downloads)
   - Check that biological interpretations are accurate

## Reporting Bugs

### Before you submit a bug report

- **Check the existing issues** — someone may have already reported it
- **Check the README** — your question might be answered there
- **Try with fresh data** — test with the included example dataset first

### How to submit a bug report

1. Go to **Issues** and click **New Issue**
2. Select the **Bug report** template
3. Fill in all sections clearly:
   - What you were trying to do
   - What actually happened
   - Your R/package versions
   - Steps to reproduce (as specific as possible)

**Example:**
- R version: 4.3.0
- DESeq2 version: 1.40.0
- Browser: Chrome on Ubuntu 22.04
- Steps: Loaded GSE157234 data → clicked volcano plot tab → error appears in console

## Suggesting Improvements

### Before you submit an improvement suggestion

- **Check the README** — it might already be on the roadmap
- **Search issues** — someone might have suggested it already
- **Think biologically** — will this improve interpretation of results?

### How to submit a feature request

1. Go to **Issues** and click **New Issue**
2. Select the **Feature request** template
3. Explain:
   - What you'd like to add
   - Why it's useful for RNA-seq analysis
   - Example use case (optional)

**Example:**
- Request: Add pathway enrichment results visualization
- Why: Users could see KEGG or GO term enrichment in context of DE genes
- Use case: Linking top DE genes to biological processes

## Contributing Code

### Before you start coding

1. **Open an issue first** — discuss your idea with the maintainer
2. **Wait for approval** — this prevents wasted effort on incompatible changes
3. **Understand the project structure:**
   - `app.R` — main Shiny app
   - `server.R` (if separate) — server logic
   - `ui.R` (if separate) — user interface
   - `scripts/` — data processing and helper functions
   - `data/` — example datasets

### How to contribute code

1. **Fork the repository** on GitHub
2. **Create a feature branch** (not on `main`):
   ```bash
   git checkout -b fix/volcano-plot-labels
   ```
3. **Make your changes** — keep commits focused and well-documented
4. **Test thoroughly:**
   - Run the app locally
   - Test with different data inputs
   - Verify plots/tables render correctly
5. **Write clear commit messages:**
   ```
   Fix: correct NA handling in log2 fold-change calculation
   
   Previously, rows with NA values in baseMean would cause 
   the volcano plot to fail. Now filtering NAs before plotting.
   Closes #23
   ```
6. **Push to your fork** and **create a Pull Request**
7. **Fill in the PR template** with:
   - What your changes do
   - How to test them
   - Any breaking changes

### Code Style

We follow these conventions:

- **R files:** tidyverse style (see `style_guide.md` if present)
- **Function names:** `snake_case` for functions, `camelCase` for variables
- **Comments:** Explain *why*, not *what*. Good: `# Filter by adjusted p-value cutoff (0.05)`. Avoid: `# remove rows`
- **Biological accuracy:** Always verify statistical interpretations with the literature

### Example contribution: Adding a new filter

```r
# Good: explains biological rationale
# Filter by log2 fold-change threshold to focus on 
# genes with meaningful biological effect size
filtered_genes <- data[abs(data$log2FoldChange) > lfc_cutoff, ]

# Bad: unclear
genes_filt <- data[abs(data$lfc) > x, ]
```

## Commit Message Guidelines

Follow conventional commits (simplified version):

- `fix:` — bug fixes
- `feat:` — new features
- `docs:` — documentation only
- `style:` — formatting (no code changes)
- `refactor:` — code restructuring (no behavior change)

**Format:**
```
type: brief description (50 chars max)

Optional: longer explanation of why this change matters.
Include issue number if applicable: Closes #42
```

## Pull Request Process

1. Update the README if you're adding features
2. Update the CHANGELOG if one exists
3. Ensure the Shiny app runs without warnings/errors
4. Your PR will be reviewed within 1-2 weeks
5. Address review feedback — maintainers may suggest changes
6. Once approved, your PR will be merged

## Questions or Discussions?

- **General questions** — Start a Discussion (if enabled)
- **Urgent issues** — Email the maintainer
- **Bioinformatics questions** — Ask in the issue description with context

## Recognition

All contributors will be acknowledged in:
- The GitHub contributors page
- README.md (if significant contribution)

## Code of Conduct

This project adheres to the Contributor Covenant Code of Conduct. By participating, you agree to uphold this code. See `CODE_OF_CONDUCT.md` for details.

---

**Thank you for making this project better!** 🎉

Your contributions, whether code, bug reports, or feedback, help improve reproducibility and accessibility of RNA-seq analysis. Open science thrives when researchers can explore, understand, and build upon each other's work.
