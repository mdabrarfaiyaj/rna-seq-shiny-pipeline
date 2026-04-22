# ============================================================
# app_v2.R — RNA-Seq DESeq2 Interactive Shiny Dashboard
# VERSION 2 — True Raw Counts from Galaxy featureCounts
# ============================================================
# Dataset    : GSE157234 — Shemer et al., Immunity 2020
# Comparison : IL10R-Mutant vs Control Microglia (48h post-LPS)
# Author     : Md. Abrar Faiyaj
# GitHub     : https://github.com/mdabrarfaiyaj/rna-seq-shiny-pipeline
#
# WHAT CHANGED FROM V1:
#   V1 used UTAP-normalized counts from GEO (pre-normalized,
#   caused double-normalization in DESeq2).
#   V2 uses true raw integer counts from Galaxy featureCounts
#   (STAR → featureCounts → DESeq2 — correct pipeline).
#
# PREREQUISITE:
#   Run analysis_v2.R first. It saves RDS objects to results/v2/
#   and clean data to data/v2/. This app loads those files.
#
# HOW TO RUN:
#   1. Run analysis_v2.R completely
#   2. Open app_v2.R in RStudio and click Run App
# ============================================================

shiny_pkgs <- c("shiny","shinydashboard","DT","plotly",
                "ggrepel","pheatmap","dplyr","RColorBrewer",
                "ggplot2","DESeq2","apeglm")
for (pkg in shiny_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}

library(shiny);        library(shinydashboard)
library(DESeq2);       library(ggplot2);  library(ggrepel)
library(pheatmap);     library(dplyr);    library(DT)
library(RColorBrewer); library(plotly);library(apeglm)


# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("dna"), " RNA-Seq Explorer V2"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Home",          tabName="home",    icon=icon("house")),
      menuItem("Upload Data",   tabName="upload",  icon=icon("upload")),
      menuItem("Volcano Plot",  tabName="volcano", icon=icon("volcano")),
      menuItem("PCA Plot",      tabName="pca",     icon=icon("circle-dot")),
      menuItem("Heatmap",       tabName="heatmap", icon=icon("table-cells")),
      menuItem("Results Table", tabName="table",   icon=icon("table")),
      menuItem("About",         tabName="about",   icon=icon("circle-info"))
    ),
    hr(),
    div(style="padding:10px;",
        h5(icon("sliders")," Parameters",
           style="color:white;font-weight:bold;"),
        sliderInput("padj_cutoff","Adj. P-value cutoff:",
                    value=0.05, min=0.001, max=0.1,  step=0.001),
        sliderInput("lfc_cutoff", "Log2 FC cutoff:",
                    value=1.0,  min=0.5,   max=3.0,  step=0.1),
        sliderInput("top_n_genes","Top genes to label:",
                    value=20,   min=5,     max=50,   step=5),
        sliderInput("top_n_heat", "Genes in heatmap:",
                    value=50,   min=10,    max=100,  step=10)
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper{background-color:#f5f7fa;}
      .box{border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.08);}
      .stat-card{background:white;border-radius:8px;padding:20px;
                 text-align:center;box-shadow:0 2px 8px rgba(0,0,0,0.08);}
      .v2-badge{background:#27ae60;color:white;padding:3px 8px;
                border-radius:4px;font-size:11px;font-weight:bold;}
    "))),

    tabItems(

      # HOME
      tabItem("home",
              fluidRow(
                box(width=12, status="primary",
                    h2(icon("dna")," RNA-Seq Differential Expression Dashboard",
                       style="color:#2C3E50;font-weight:bold;"),
                    span(class="v2-badge","V2 — True Raw Counts"),
                    br(),br(),
                    p("Professional RNA-Seq analysis using ",strong("DESeq2"),
                      " — the gold standard for differential expression."),
                    p("Dataset: ",strong("Shemer et al., Immunity 2020"),
                      " | Comparison: ",
                      strong("IL10R-Mutant vs Control microglia at 48h post-LPS.")),
                    div(style="background:#eaf4fb;padding:12px;border-radius:6px;
                               border-left:4px solid #3498db;margin:10px 0;",
                        icon("circle-info",style="color:#3498db;"),
                        strong(" V2 improvement: "),
                        "This version uses true raw integer counts from Galaxy ",
                        "featureCounts ( HISAT2 → featureCounts → DESeq2). ",
                        "V1 used UTAP-normalized counts which caused double-normalization. ",
                        "V2 is methodologically correct."
                    ),
                    hr(),
                    fluidRow(
                      column(4,div(class="stat-card",
                                   icon("flask",style="font-size:30px;color:#3498DB;"),
                                   br(),h4("GSE157234"),p("Real published dataset"))),
                      column(4,div(class="stat-card",
                                   icon("chart-line",style="font-size:30px;color:#E74C3C;"),
                                   br(),h4("DESeq2 + apeglm"),p("Modern correct pipeline"))),
                      column(4,div(class="stat-card",
                                   icon("upload",style="font-size:30px;color:#2ECC71;"),
                                   br(),h4("Your Data"),p("Upload your own counts")))
                    ),
                    hr(),
                    h4("How to use:"),
                    tags$ol(
                      tags$li("Run ",strong("analysis_v2.R")," in RStudio first"),
                      tags$li("Click ",strong("Upload Data")," → Load Demo (V2)"),
                      tags$li("Explore Volcano, PCA, Heatmap tabs"),
                      tags$li("Download results from Results Table")
                    )
                )
              ),
              fluidRow(uiOutput("home_stats"))
      ),

      # UPLOAD
      tabItem("upload",
              fluidRow(
                box(width=6, title="Upload Count Matrix", status="primary",
                    p("CSV: rows=genes, columns=samples, values=",
                      strong("raw integer counts")),
                    fileInput("count_file","Choose Count Matrix (.csv)",
                              accept=c(".csv",".txt")),
                    tags$small(style="color:grey;",
                               "Gene names in first column, sample IDs as header.",
                               " Must be raw counts — not normalized.")
                ),
                box(width=6, title="Upload Metadata", status="primary",
                    p("CSV: rows=samples. Must have a ",strong("condition"),
                      " column with exactly 2 groups."),
                    fileInput("meta_file","Choose Metadata (.csv)",accept=".csv")
                )
              ),
              fluidRow(
                box(width=12,
                    title="OR Load V2 Demo Dataset (GSE157234 — True Raw Counts)",
                    status="success",
                    p("Reads pre-computed files saved by ",
                      strong("analysis_v2.R"),":"),
                    tags$ul(
                      tags$li(code("results/v2/dds_object_v2.rds")),
                      tags$li(code("results/v2/vsd_object_v2.rds")),
                      tags$li(code("results/v2/res_df_v2.rds"))
                    ),
                    p(icon("triangle-exclamation",style="color:orange;"),
                      " Run ",strong("analysis_v2.R")," first if these don't exist."),
                    actionButton("load_demo","Load V2 Demo (True Raw Counts)",
                                 icon=icon("play"),class="btn-success btn-lg"),
                    br(),br(),
                    uiOutput("data_status")
                )
              ),
              fluidRow(
                box(width=12, title="Run Analysis (for custom uploads)", status="warning",
                    p("Only needed if you uploaded your own count matrix above.",
                      " Demo data loads with pre-computed results automatically."),
                    actionButton("run_analysis","Run DESeq2 Analysis",
                                 icon=icon("play"),class="btn-warning btn-lg"),
                    br(),br(),
                    verbatimTextOutput("analysis_log")
                )
              )
      ),

      # VOLCANO
      tabItem("volcano",
              fluidRow(
                box(width=12,
                    title="Volcano Plot: IL10R-Mutant vs Control (48h post-LPS)",
                    status="danger",
                    p("Red=upregulated in Mutant (Tnf, Ccl5, Il6, Il12b, Il1b…). ",
                      "Blue=downregulated in Mutant (P2ry12, Sall1, Tmem119…). ",
                      "Adjust thresholds in sidebar. apeglm LFC shrinkage applied."),
                    fluidRow(
                      column(3,downloadButton("dl_volcano_png","Download PNG",
                                              class="btn-sm")),
                      column(3,downloadButton("dl_volcano_pdf","Download PDF",
                                              class="btn-sm"))
                    ),
                    br(),
                    plotlyOutput("volcano_plot",height="600px")
                )
              ),
              fluidRow(
                valueBoxOutput("n_up",  width=4),
                valueBoxOutput("n_down",width=4),
                valueBoxOutput("n_ns",  width=4)
              )
      ),

      # PCA
      tabItem("pca",
              fluidRow(
                box(width=12, title="PCA: Sample Clustering", status="primary",
                    p("Mutant (red triangles) should separate clearly from ",
                      "Control (blue circles) along PC1 — confirms strong ",
                      "hyperactivation effect at 48h post-LPS. "),
                  
                    fluidRow(
                      column(3,downloadButton("dl_pca_png","Download PNG",
                                              class="btn-sm")),
                      column(3,downloadButton("dl_pca_pdf","Download PDF",
                                              class="btn-sm"))
                    ),
                    br(),
                    plotlyOutput("pca_plot",height="550px")
                )
              )
      ),

      # HEATMAP
      tabItem("heatmap",
              fluidRow(
                box(width=12, title="Heatmap: Top DEGs", status="warning",
                    p("Z-scored expression. Red=high, Blue=low. ",
                      "Expect Tnf, Ccl5, Il6 high in Mutant; ",
                      "P2ry12, Sall1, Tmem119 low in Mutant."),
                    fluidRow(
                      column(3,downloadButton("dl_heat_png","Download PNG",
                                              class="btn-sm")),
                      column(3,downloadButton("dl_heat_pdf","Download PDF",
                                              class="btn-sm"))
                    ),
                    br(),
                    plotOutput("heatmap_plot",height="700px")
                )
              )
      ),

      # TABLE
      tabItem("table",
              fluidRow(
                box(width=12, title="Differential Expression Results",
                    status="success",
                    fluidRow(
                      column(4,
                             selectInput("table_filter","Show:",
                                         choices=c("All Significant"  ="sig",
                                                   "Upregulated Only"  ="up",
                                                   "Downregulated Only"="down",
                                                   "All Genes"         ="all"))
                      ),
                      column(4,br(),
                             downloadButton("dl_results_csv","Download CSV",
                                            class="btn-success"))
                    ),
                    br(),
                    DTOutput("results_table")
                )
              )
      ),

      # ABOUT
      tabItem("about",
              fluidRow(
                box(width=8, title="About This App", status="info",
                    h4("RNA-Seq Differential Expression Dashboard — V2"),
                    p("Built with R + Shiny + DESeq2"),
                    hr(),
                    h4("Dataset"),
                    p(strong("GEO:"),"GSE157234"),
                    p(strong("Paper:"),"Shemer et al., Immunity 53, 1033–1049, 2020"),
                    p(strong("Comparison:"),"IL10R-Mutant vs Control microglia"),
                    p(strong("Timepoint:"),"48h post-LPS (peak hyperactivation)"),
                    p(strong("Key finding:"),
                      "Without IL-10 signalling, microglia hyperactivate and produce",
                      "toxic TNF, causing neuronal damage (Figure 3, paper)."),
                    hr(),
                    h4("V2 Pipeline"),
                    tags$ul(
                      tags$li("Input: True raw integer counts from Galaxy featureCounts"),
                      tags$li("Alignment: STAR → mm10 (NCBI RefSeq)"),
                      tags$li("Quantification: featureCounts (all exons)"),
                      tags$li("Samples: 6 Control + 3 Mutant (48h post-LPS)"),
                      tags$li("DESeq2 design: ~ condition"),
                      tags$li("Filter: genes ≥10 counts in ≥2 samples"),
                      tags$li("LFC shrinkage: apeglm (modern — replaces deprecated betaPrior)"),
                      tags$li("Significance: padj<0.05 AND |log2FC|>1"),
                      tags$li("VST transformation for PCA and heatmap")
                    ),
                    hr(),
                    h4("V1 vs V2 comparison"),
                    tags$table(
                      style="width:100%;border-collapse:collapse;",
                      tags$tr(style="background:#eee;",
                              tags$th(style="padding:6px;",""),
                              tags$th(style="padding:6px;","V1"),
                              tags$th(style="padding:6px;","V2")),
                      tags$tr(tags$td(style="padding:6px;","Input"),
                              tags$td(style="padding:6px;color:orange;","UTAP-normalized"),
                              tags$td(style="padding:6px;color:green;","True raw counts")),
                      tags$tr(style="background:#f9f9f9;",
                              tags$td(style="padding:6px;","Normalization"),
                              tags$td(style="padding:6px;color:orange;","Double-normalized"),
                              tags$td(style="padding:6px;color:green;","Once — by DESeq2")),
                      tags$tr(tags$td(style="padding:6px;","LFC method"),
                              tags$td(style="padding:6px;","Standard"),
                              tags$td(style="padding:6px;color:green;","apeglm shrinkage")),
                      tags$tr(style="background:#f9f9f9;",
                              tags$td(style="padding:6px;","DEGs up/down"),
                              tags$td(style="padding:6px;","669 / 894"),
                              tags$td(style="padding:6px;","621 / 976")),
                      tags$tr(tags$td(style="padding:6px;","Paper (Fig 3E)"),
                              tags$td(style="padding:6px;","954 / 693"),
                              tags$td(style="padding:6px;","954 / 693"))
                    ),
                    br(),
                    p(style="color:grey;font-size:12px;",
                      "V2 count differences from paper are due to annotation differences: ",
                      "paper used Gencode vM10 with MARS-seq 3'UTR counting; V2 uses ",
                      "NCBI RefSeq with standard exon counting. Biological conclusions ",
                      "are consistent — key markers confirmed in correct direction.")
                ),
                box(width=4, title="Developer", status="success",
                    h4("Md. Abrar Faiyaj"),
                    p("MSc Biotechnology | Junior Research Collaborator"),
                    p("BRAC University, Dhaka, Bangladesh"),
                    p(icon("github"),
                      a("github.com/mdabrarfaiyaj",
                        href="https://github.com/mdabrarfaiyaj",
                        target="_blank")),
                    hr(),
                    p(strong("This pipeline demonstrates:")),
                    tags$ul(
                      tags$li("Bulk RNA-seq: SRA → STAR → featureCounts → DESeq2"),
                      tags$li("Interactive Shiny dashboard development"),
                      tags$li("Reproducible bioinformatics pipeline design"),
                      tags$li("Honest limitation documentation"),
                      tags$li("V1 → V2 methodological improvement")
                    ),
                    hr(),
                    p(strong("Zenodo DOI:")),
                    p(a("10.5281/zenodo.19138922",
                        href="https://doi.org/10.5281/zenodo.19138922",
                        target="_blank"))
                )
              )
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  rv <- reactiveValues(
    counts=NULL, metadata=NULL, dds=NULL,
    res_df=NULL, vsd=NULL, ready=FALSE,
    log_text=paste0(
      "Waiting...\n\n",
      "Step 1: Run analysis_v2.R in RStudio\n",
      "Step 2: Click 'Load V2 Demo Dataset' above\n",
      "Step 3: Explore the plot tabs\n\n",
      "For custom data: upload count matrix + metadata,\n",
      "then click Run DESeq2 Analysis."
    )
  )

  # ---- Load V2 Demo: RDS priority loading -------------------
  # WHY RDS PRIORITY?
  #   RDS files from analysis_v2.R load in ~1 second.
  #   Falls back to CSV if RDS files are missing.
  # V2 RDS paths: results/v2/ (separate from V1 results/)
  observeEvent(input$load_demo, {
    withProgress(message="Loading V2 demo data...", {
      tryCatch({

        rds_dds    <- "results/v2/dds_object_v2.rds"
        rds_vsd    <- "results/v2/vsd_object_v2.rds"
        rds_res_df <- "results/v2/res_df_v2.rds"
        count_file <- "data/v2/count_matrix_raw_v2.csv"
        meta_file  <- "data/v2/metadata_v2.csv"

        # ---- PATH 1: RDS files exist — load instantly --------
        if (all(file.exists(rds_dds, rds_vsd, rds_res_df))) {

          incProgress(0.3, "Loading pre-computed V2 DESeq2 objects...")
          dds_loaded    <- readRDS(rds_dds)
          vsd_loaded    <- readRDS(rds_vsd)
          res_df_loaded <- readRDS(rds_res_df)

          incProgress(0.7, "Restoring reactive values...")
          rv$counts   <- counts(dds_loaded)
          rv$metadata <- as.data.frame(colData(dds_loaded))
          rv$dds      <- dds_loaded
          rv$vsd      <- vsd_loaded
          rv$res_df   <- res_df_loaded
          rv$ready    <- TRUE

          n_up   <- sum(res_df_loaded$significance == "Up in Mutant",   na.rm=TRUE)
          n_down <- sum(res_df_loaded$significance == "Down in Mutant",  na.rm=TRUE)
          n_ctrl <- sum(rv$metadata$condition == "Control")
          n_mut  <- sum(rv$metadata$condition == "Mutant")

          rv$log_text <- paste0(
            "✓ V2 results loaded instantly (true raw counts)!\n",
            "  Source: RDS objects from analysis_v2.R\n\n",
            "  Input:        True raw counts (Galaxy featureCounts)\n",
            "  Genes tested: ", nrow(res_df_loaded), "\n",
            "  Control:      ", n_ctrl, " samples (48h post-LPS)\n",
            "  Mutant:       ", n_mut,  " samples (48h post-LPS)\n",
            "  Up in Mutant: ", n_up,   " (padj<0.05, |LFC|>1)\n",
            "  Down in Mutant:", n_down, " (padj<0.05, |LFC|>1)\n",
            "  LFC method:   apeglm shrinkage\n\n",
            "→ Plots are ready! Click any tab above."
          )
          showNotification("✓ V2 results loaded (true raw counts)!",
                           type="message")

          # ---- PATH 2: Fallback — CSV loading ------------------
        } else {

          if (!file.exists(count_file))
            stop(paste("File not found:", count_file,
                       "\nPlease run analysis_v2.R first."))
          if (!file.exists(meta_file))
            stop(paste("File not found:", meta_file,
                       "\nPlease run analysis_v2.R first."))

          incProgress(0.3, "Reading V2 count matrix from CSV...")
          mat <- read.csv(count_file, row.names=1, check.names=FALSE)
          mat <- as.matrix(mat)
          storage.mode(mat) <- "integer"
          mat[mat < 0] <- 0L

          incProgress(0.6, "Reading metadata from CSV...")
          meta <- read.csv(meta_file, row.names=1, stringsAsFactors=FALSE)
          meta$condition <- factor(meta$condition,
                                   levels=c("Control","Mutant"))

          shared <- intersect(colnames(mat), rownames(meta))
          mat    <- mat[, shared]
          meta   <- meta[shared, , drop=FALSE]

          rv$counts   <- mat
          rv$metadata <- meta
          rv$ready    <- FALSE

          n_ctrl <- sum(meta$condition=="Control")
          n_mut  <- sum(meta$condition=="Mutant")

          rv$log_text <- paste0(
            "✓ V2 demo data loaded from CSV (fallback mode)\n",
            "  Genes:   ", nrow(mat),"\n",
            "  Samples: ", ncol(mat),"\n",
            "  Control: ", n_ctrl," samples (48h post-LPS)\n",
            "  Mutant:  ", n_mut, " samples (48h post-LPS)\n\n",
            "⚠ RDS files not found — click 'Run DESeq2 Analysis' to proceed.\n",
            "  Tip: Run Block 13 of analysis_v2.R to generate RDS files."
          )
          showNotification("V2 demo loaded (CSV fallback — click Run DESeq2)",
                           type="warning")
        }

        incProgress(1.0)

      }, error=function(e) {
        rv$log_text <- paste0("ERROR: ", e$message,
                              "\n\nFix: Run analysis_v2.R first.")
        showNotification(paste("Error:", e$message),
                         type="error", duration=15)
      })
    })
  })

  # ---- Upload count file ----
  observeEvent(input$count_file, {
    req(input$count_file)
    tryCatch({
      df  <- read.csv(input$count_file$datapath,
                      row.names=1, check.names=FALSE)
      mat <- as.matrix(df)
      storage.mode(mat) <- "integer"
      rv$counts   <- mat
      rv$log_text <- paste0("✓ Count matrix uploaded!\n",
                            "  Genes: ",nrow(mat),
                            "  Samples: ",ncol(mat),
                            "\n\nUpload metadata then click Run.")
      showNotification("Count matrix loaded!", type="message")
    }, error=function(e) {
      showNotification(paste("Error:",e$message), type="error")
    })
  })

  # ---- Upload metadata ----
  observeEvent(input$meta_file, {
    req(input$meta_file)
    tryCatch({
      meta <- read.csv(input$meta_file$datapath, row.names=1)
      if (!"condition" %in% colnames(meta))
        stop("Metadata must have a 'condition' column!")
      meta$condition <- factor(meta$condition)
      rv$metadata <- meta
      showNotification("Metadata loaded!", type="message")
    }, error=function(e) {
      showNotification(paste("Error:",e$message), type="error")
    })
  })

  # ---- Run DESeq2 (for custom uploads) ----
  observeEvent(input$run_analysis, {
    req(rv$counts, rv$metadata)
    withProgress(message="Running DESeq2...", value=0, {
      tryCatch({
        shared     <- intersect(colnames(rv$counts), rownames(rv$metadata))
        counts_sub <- rv$counts[, shared]
        meta_sub   <- rv$metadata[shared, , drop=FALSE]

        na_rows    <- rowSums(is.na(counts_sub)) > 0
        counts_sub <- counts_sub[!na_rows, ]

        incProgress(0.1, "Creating DESeq2 object...")
        dds  <- DESeqDataSetFromMatrix(countData=counts_sub,
                                       colData=meta_sub,
                                       design=~condition)
        keep <- rowSums(counts(dds) >= 10) >= 2
        dds  <- dds[keep, ]

        incProgress(0.3, "Running DESeq2 (1-3 min)...")
        set.seed(123)
        dds <- DESeq(dds)

        incProgress(0.6, "Extracting results with apeglm shrinkage...")
        lvls <- levels(meta_sub$condition)

        res <- results(dds,
                       contrast = c("condition", lvls[2], lvls[1]),
                       alpha    = 0.05)

        res_shrunk <- lfcShrink(dds,
                                coef = paste0("condition_", lvls[2], "_vs_", lvls[1]),
                                type = "apeglm")

        res_df      <- as.data.frame(res_shrunk)
        res_df$gene <- rownames(res_df)

        res_df$significance <- dplyr::case_when(
          !is.na(res_df$padj) &
            res_df$padj < input$padj_cutoff &
            res_df$log2FoldChange > input$lfc_cutoff  ~ "Up in Mutant",
          !is.na(res_df$padj) &
            res_df$padj < input$padj_cutoff &
            res_df$log2FoldChange < -input$lfc_cutoff ~ "Down in Mutant",
          TRUE ~ "Not Significant"
        )

        incProgress(0.8, "VST transformation...")
        vsd <- vst(dds, blind=FALSE)

        rv$dds=dds; rv$res_df=res_df; rv$vsd=vsd; rv$ready=TRUE

        n_up   <- sum(res_df$significance=="Up in Mutant",   na.rm=TRUE)
        n_down <- sum(res_df$significance=="Down in Mutant",  na.rm=TRUE)

        rv$log_text <- paste0(
          "✓ DESeq2 complete (apeglm shrinkage)!\n",
          "  Comparison:    ",lvls[2]," vs ",lvls[1],"\n",
          "  Genes tested:  ",nrow(res_df),"\n",
          "  Up in Mutant:  ",n_up, " (padj<0.05, |LFC|>1)\n",
          "  Down in Mutant:",n_down," (padj<0.05, |LFC|>1)\n\n",
          "→ Explore Volcano / PCA / Heatmap tabs!"
        )
        showNotification("Analysis complete!", type="message")
        incProgress(1.0)

      }, error=function(e) {
        rv$log_text <- paste("ERROR:",e$message)
        showNotification(paste("Failed:",e$message),
                         type="error", duration=15)
      })
    })
  })

  # ---- UI helpers ----
  output$data_status <- renderUI({
    if (!is.null(rv$counts) && !is.null(rv$metadata)) {
      div(style="color:green;font-weight:bold;margin-top:10px;",
          icon("circle-check"),
          paste(" Ready:",nrow(rv$counts),"genes,",
                ncol(rv$counts),"samples |",
                sum(rv$metadata$condition=="Control"),"Control,",
                sum(rv$metadata$condition=="Mutant"),"Mutant"))
    } else {
      div(style="color:orange;margin-top:10px;",
          icon("clock")," Waiting for data...")
    }
  })

  output$analysis_log <- renderText({ rv$log_text })

  output$home_stats <- renderUI({
    if (!rv$ready) return(NULL)
    n_up  <- sum(rv$res_df$significance=="Up in Mutant",   na.rm=TRUE)
    n_down<- sum(rv$res_df$significance=="Down in Mutant",  na.rm=TRUE)
    fluidRow(
      infoBox("Up in Mutant",   n_up,           icon=icon("arrow-up"),  color="red"),
      infoBox("Down in Mutant", n_down,          icon=icon("arrow-down"),color="blue"),
      infoBox("Genes Tested",   nrow(rv$res_df), icon=icon("dna"),       color="green")
    )
  })

  output$n_up  <- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Up in Mutant",   na.rm=TRUE) else "—"
    valueBox(n,"Up in Mutant",  icon=icon("arrow-up"),  color="red")
  })
  output$n_down<- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Down in Mutant", na.rm=TRUE) else "—"
    valueBox(n,"Down in Mutant",icon=icon("arrow-down"),color="blue")
  })
  output$n_ns  <- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Not Significant", na.rm=TRUE) else "—"
    valueBox(n,"Not Significant",icon=icon("minus"),    color="gray")
  })

  # ---- Volcano ----
  make_volcano <- reactive({
    req(rv$ready)
    df <- rv$res_df |>
      filter(!is.na(padj), !is.na(log2FoldChange)) |>
      mutate(
        sig = dplyr::case_when(
          padj < input$padj_cutoff & log2FoldChange >  input$lfc_cutoff ~ "Up in Mutant",
          padj < input$padj_cutoff & log2FoldChange < -input$lfc_cutoff ~ "Down in Mutant",
          TRUE ~ "Not Significant"),
        neg_log10_padj = -log10(padj + 1e-300)
      )
    top_lab <- df |> filter(sig != "Not Significant") |>
      arrange(padj) |> head(input$top_n_genes)

    ggplot(df, aes(x=log2FoldChange, y=neg_log10_padj,
                   color=sig, text=gene)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=c("Up in Mutant"    ="#E74C3C",
                                  "Down in Mutant"  ="#3498DB",
                                  "Not Significant" ="grey70")) +
      geom_vline(xintercept=c(-input$lfc_cutoff, input$lfc_cutoff),
                 linetype="dashed", color="grey40") +
      geom_hline(yintercept=-log10(input$padj_cutoff),
                 linetype="dashed", color="grey40") +
      geom_text_repel(data=top_lab, aes(label=gene),
                      size=2.8, max.overlaps=25, color="black") +
      labs(title   ="Volcano Plot: IL10R-Mutant vs Control (48h post-LPS)",
           subtitle=paste0("V2 — True Raw Counts | apeglm shrinkage | ",
                           "padj<", input$padj_cutoff,
                           "  |log2FC|>", input$lfc_cutoff),
           x="Log2 Fold Change (Mutant / Control)",
           y="-Log10 Adjusted P-value",
           color="Regulation") +
      theme_bw(base_size=13) +
      theme(plot.title=element_text(face="bold"),
            panel.grid.minor=element_blank())
  })

  output$volcano_plot <- renderPlotly({
    req(rv$ready); ggplotly(make_volcano(), tooltip=c("text","x","y"))
  })
  output$dl_volcano_png <- downloadHandler(
    filename="volcano_v2.png",
    content=function(f) ggsave(f, make_volcano(), width=10, height=7, dpi=300))
  output$dl_volcano_pdf <- downloadHandler(
    filename="volcano_v2.pdf",
    content=function(f) ggsave(f, make_volcano(), width=10, height=7))

  # ---- PCA ----
  make_pca <- reactive({
    req(rv$ready)
    pca_data <- plotPCA(rv$vsd, intgroup="condition", returnData=TRUE)
    # Flip PC1 axis for visual consistency with paper Figure 3B
    pca_data$PC1 <- -pca_data$PC1
    pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

    ggplot(pca_data, aes(x=PC1, y=PC2, color=condition,
                         shape=condition, text=name)) +
      geom_point(size=4, alpha=0.85) +
      scale_color_manual(values=c("Control"="#4575b4",
                                  "Mutant" ="#d73027")) +
      scale_shape_manual(values=c("Control"=16, "Mutant"=17)) +
      stat_ellipse(aes(group=condition), linetype="dashed", level=0.8) +
      geom_text_repel(aes(label=name), size=2.5, max.overlaps=20) +
      labs(title   ="PCA: IL10R-Mutant vs Control (48h post-LPS)",
           subtitle="V2 — True Raw Counts | VST-transformed",
           x=paste0("PC1: ", pct_var[1], "% variance"),
           y=paste0("PC2: ", pct_var[2], "% variance"),
           caption ="PC1 axis direction flipped for visual consistency with Figure 3B") +
      theme_bw(base_size=13) +
      theme(plot.title=element_text(face="bold"),
            panel.grid.minor=element_blank())
  })

  output$pca_plot <- renderPlotly({
    req(rv$ready); ggplotly(make_pca(), tooltip=c("text","x","y"))
  })
  output$dl_pca_png <- downloadHandler(
    filename="pca_v2.png",
    content=function(f) ggsave(f, make_pca(), width=9, height=6, dpi=300))
  output$dl_pca_pdf <- downloadHandler(
    filename="pca_v2.pdf",
    content=function(f) ggsave(f, make_pca(), width=9, height=6))

  # ---- Heatmap ----
  make_heatmap <- reactive({
    req(rv$ready)
    top_genes <- rv$res_df |>
      filter(!is.na(padj), significance != "Not Significant") |>
      arrange(padj) |> head(input$top_n_heat) |> pull(gene)
    if (length(top_genes) < 5)
      top_genes <- rv$res_df |> filter(!is.na(pvalue)) |>
        arrange(pvalue) |> head(input$top_n_heat) |> pull(gene)

    valid_genes    <- top_genes[top_genes %in% rownames(rv$vsd)]
    heat_mat       <- assay(rv$vsd)[valid_genes, ]
    heat_scaled    <- t(scale(t(heat_mat)))
    heat_scaled[heat_scaled >  3] <-  3
    heat_scaled[heat_scaled < -3] <- -3

    col_annot <- data.frame(
      Condition=as.character(rv$metadata[colnames(heat_mat), "condition"]),
      row.names=colnames(heat_mat))
    annot_colors <- list(Condition=c("Control"="#4575b4",
                                     "Mutant" ="#d73027"))

    pheatmap(heat_scaled,
             annotation_col   = col_annot,
             annotation_colors= annot_colors,
             color            = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
             cluster_rows=TRUE, cluster_cols=TRUE,
             show_rownames=TRUE,
             fontsize_row=max(5, 9 - input$top_n_heat/20),
             main=paste0("Top ", input$top_n_heat,
                         " DEGs | Mutant vs Control | 48h | V2"),
             border_color=NA, silent=TRUE)
  })

  output$heatmap_plot <- renderPlot({
    req(rv$ready)
    grid::grid.newpage(); grid::grid.draw(make_heatmap()$gtable)
  })
  output$dl_heat_png <- downloadHandler(
    filename="heatmap_v2.png",
    content=function(f){
      png(f, width=2400, height=3200, res=300)
      grid::grid.draw(make_heatmap()$gtable); dev.off()})
  output$dl_heat_pdf <- downloadHandler(
    filename="heatmap_v2.pdf",
    content=function(f){
      pdf(f, width=10, height=13)
      grid::grid.draw(make_heatmap()$gtable); dev.off()})

  # ---- Results Table ----
  filtered_res <- reactive({
    req(rv$ready)
    df <- rv$res_df |>
      filter(!is.na(padj)) |>
      mutate(across(where(is.numeric), ~round(., 4)))
    switch(input$table_filter,
           "sig"  = df |> filter(significance != "Not Significant"),
           "up"   = df |> filter(significance == "Up in Mutant"),
           "down" = df |> filter(significance == "Down in Mutant"),
           "all"  = df)
  })

  output$results_table <- renderDT({
    req(rv$ready)
    datatable(
      filtered_res() |> select(gene, baseMean, log2FoldChange,
                                pvalue, padj, significance),
      options=list(pageLength=25, scrollX=TRUE,
                   order=list(list(4, "asc"))),
      rownames=FALSE
    ) |> formatStyle("significance",
                     backgroundColor=styleEqual(
                       c("Up in Mutant","Down in Mutant","Not Significant"),
                       c("#FADBD8","#D6EAF8","#F8F9F9")))
  })

  output$dl_results_csv <- downloadHandler(
    filename="DESeq2_results_v2.csv",
    content=function(f) write.csv(filtered_res(), f, row.names=FALSE))
}

shinyApp(ui=ui, server=server)
