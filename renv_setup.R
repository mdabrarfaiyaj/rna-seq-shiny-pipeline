# ============================================================
# REPRODUCIBLE ENVIRONMENT SETUP
# RNA-Seq DESeq2 Pipeline | GSE157234 | Shemer et al., 2020
# ============================================================
# PURPOSE:
#   This script sets up renv for environment reproducibility
#   and regenerates manifest.json for Posit Cloud deployment.
#
# WHAT IS renv AND WHY DOES IT MATTER?
#   renv (authored by Kevin Ushey & Hadley Wickham, Posit) is
#   the industry-standard dependency management tool for R.
#   Source: Ushey & Wickham (2024) https://rstudio.github.io/renv/
#
#   Without renv:
#     A researcher cloning your repo in 2027 might install
#     DESeq2 v1.52 instead of v1.46. DESeq2 updates can change
#     dispersion estimation, resulting in different p-values
#     and different DEG lists — silently, with no error.
#
#   With renv:
#     renv.lock records the EXACT version of every package.
#     Anyone runs renv::restore() and gets identical results.
#     This is required by Nature Methods, Genome Biology, and
#     most top bioinformatics journals for code submissions.
#
#   For Bioconductor specifically:
#     renv::init(bioconductor = TRUE) is required because
#     Bioconductor packages (DESeq2, GEOquery) are versioned
#     differently from CRAN packages. Each Bioconductor release
#     ties to a specific R version — renv captures this.
#     Source: Appsilon (2024) https://www.appsilon.com/post/renv-bioconductor
#
# HOW TO USE THIS SCRIPT:
#   RUN ONCE when setting up the project for the first time.
#   After that, renv manages itself automatically.
#
# ============================================================


# ============================================================
# SECTION 1: FIRST-TIME renv SETUP
# ============================================================
# Run this section ONCE to initialise renv for your project.
# Comment it out after first run — do not run every time.
# ============================================================

# Step 1: Install renv if not present
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
  cat("✓ renv installed\n")
} else {
  cat("✓ renv already installed\n")
}

# Step 2: Confirm working directory is the project folder, not home
# renv refuses to initialise in the home directory as a safety guard.
# Before running this script, set your working directory to the project root:
#   In RStudio: Session → Set Working Directory → To Source File Location
# The check below will stop you with a clear message if this is wrong.

home_dir <- path.expand("~")
current_dir <- getwd()

if (current_dir == home_dir) {
  stop(
    "\n\n",
    "=========================================================\n",
    "  ERROR: Working directory is your home folder.\n",
    "  renv will not initialise here — it would affect all\n",
    "  your R projects, not just this one.\n\n",
    "  Fix: In RStudio, go to:\n",
    "  Session → Set Working Directory → To Source File Location\n",
    "  Then run this script again.\n",
    "=========================================================\n"
  )
}

cat("\n✓ Working directory confirmed:", current_dir, "\n")

# Step 3: Check if renv is already initialised
if (file.exists("renv.lock")) {
  
  cat("\n✓ renv.lock already exists\n")
  cat("  Your environment is already locked.\n")
  cat("  To update after adding new packages, run:\n")
  cat("    renv::snapshot()\n\n")
  
} else {
  
  cat("\nInitialising renv for this project...\n")
  cat("bioconductor = TRUE ensures DESeq2 and GEOquery\n")
  cat("versions are locked alongside CRAN packages.\n\n")
  
  # This creates:
  #   renv/            — project-local package library
  #   renv.lock        — lockfile recording exact versions
  #   .Rprofile        — auto-activates renv for this project
  renv::init(bioconductor = TRUE)
  
  cat("\n✓ renv initialised\n")
  cat("  Now install your packages as normal, then run:\n")
  cat("    renv::snapshot()\n")
  cat("  to lock the current versions into renv.lock\n")
}


# ============================================================
# SECTION 2: SNAPSHOT CURRENT ENVIRONMENT
# ============================================================
# Run this after installing or updating any package.
# It updates renv.lock with the current package versions.
# Commit renv.lock to GitHub after every snapshot.
# ============================================================

cat("\n--- Snapshotting current package environment ---\n")
cat("This records exact versions of all packages into renv.lock\n\n")

# Uncomment the line below and run when ready to snapshot:
# renv::snapshot()

cat("⚠ renv::snapshot() is commented out — uncomment to run\n")
cat("  Run it after all packages are installed and working.\n")


# ============================================================
# SECTION 3: VERIFY ENVIRONMENT STATUS
# ============================================================
# Run this anytime to check if your installed packages match
# what is recorded in renv.lock
# ============================================================

cat("\n--- Checking environment status ---\n")

if (file.exists("renv.lock")) {
  cat("Checking if installed packages match renv.lock...\n\n")
  renv::status()
} else {
  cat("renv.lock not yet created — run Section 1 first\n")
}


# ============================================================
# SECTION 4: REGENERATE manifest.json FOR POSIT CLOUD
# ============================================================
# manifest.json tells Posit Connect Cloud exactly which
# packages and versions your Shiny app needs to run.
# It must be regenerated whenever you:
#   - Update any package
#   - Add new packages to the app
#   - Change R version
#   - Update the app code significantly
#
# RStudio regenerates this automatically when you click the
# blue Publish button. But if you need to do it manually
# (e.g. for CI/CD pipelines), use the code below.
# ============================================================

cat("\n--- Regenerating manifest.json for Posit Cloud ---\n")

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

library(rsconnect)

# Set working directory to deploy/ folder
# where app.R and data/ live
deploy_dir <- "deploy"   # relative to project root

if (dir.exists(deploy_dir)) {
  
  original_wd <- getwd()
  setwd(deploy_dir)
  
  cat("Generating manifest.json in", deploy_dir, "...\n")
  rsconnect::writeManifest()
  
  setwd(original_wd)
  cat("✓ manifest.json regenerated in", deploy_dir, "\n")
  cat("  This file is ready for Posit Cloud deployment.\n")
  cat("  Commit it to GitHub before republishing.\n")
  
} else {
  cat("⚠ deploy/ folder not found\n")
  cat("  Make sure you are running this from the project root:\n")
  cat("  ", getwd(), "\n")
}


# ============================================================
# SECTION 5: COMMIT CHECKLIST BEFORE REPUBLISHING
# ============================================================
cat("\n==========================================\n")
cat("  PRE-PUBLISH CHECKLIST\n")
cat("==========================================\n")
cat("Before clicking Publish in RStudio:\n\n")
cat("  ☐ Run analysis_final.R completely\n")
cat("    (generates RDS files + fresh CSVs)\n\n")
cat("  ☐ Copy fresh data files to deploy/:\n")
cat("    cp data/count_matrix_48h_clean.csv deploy/data/\n")
cat("    cp data/metadata_48h_clean.csv     deploy/data/\n")
cat("    cp results/res_df.rds              deploy/results/\n\n")
cat("  ☐ Run Section 4 above to regenerate manifest.json\n\n")
cat("  ☐ Test app locally with deploy/app.R open in RStudio\n\n")
cat("  ☐ Click blue Publish button → select existing app\n\n")
cat("  ☐ Verify live URL after deployment completes:\n")
cat("    https://019cd22f-689b-acc2-4b56-472725ef4a7b.share.connect.posit.cloud\n\n")
cat("  ☐ Commit everything to GitHub:\n")
cat("    git add .\n")
cat("    git commit -m 'your message'\n")
cat("    git push\n")
cat("==========================================\n")
