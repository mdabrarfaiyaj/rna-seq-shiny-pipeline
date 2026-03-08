options(repos = BiocManager::repositories())

# ============================================================
# app_final.R - RNA-Seq DESeq2 Shiny Dashboard
# Dataset: GSE157234 | Shemer et al., Immunity 2020
#
# KEY FIX: Demo button now reads pre-saved CSV files from
#          analysis_final.R instead of re-processing RAW.tar.
#
# PREREQUISITE: Run analysis_final.R first so these exist:
#   data/count_matrix_48h_clean.csv
#   data/metadata_48h_clean.csv
# ============================================================



shiny_pkgs <- c("shiny","shinydashboard","DT","plotly",
                "ggrepel","pheatmap","dplyr","RColorBrewer",
                "ggplot2","DESeq2")
for (pkg in shiny_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}

library(shiny);        library(shinydashboard)
library(DESeq2);       library(ggplot2);  library(ggrepel)
library(pheatmap);     library(dplyr);    library(DT)
library(RColorBrewer); library(plotly)


# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = span(icon("dna"), " RNA-Seq Explorer"),
    titleWidth = 260
  ),
  
  dashboardSidebar(
    width = 260,
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
    "))),
    
    tabItems(
      
      # HOME
      tabItem("home",
              fluidRow(
                box(width=12, status="primary",
                    h2(icon("dna")," RNA-Seq Differential Expression Dashboard",
                       style="color:#2C3E50;font-weight:bold;"),
                    p("Professional RNA-Seq analysis using ",strong("DESeq2"),
                      " — the gold standard."),
                    p("Dataset: ",strong("Shemer et al., Immunity 2020"),
                      " | Comparison: ",
                      strong("IL10R-Mutant vs Control microglia at 48h post-LPS.")),
                    hr(),
                    fluidRow(
                      column(4,div(class="stat-card",
                                   icon("flask",style="font-size:30px;color:#3498DB;"),
                                   br(),h4("GSE157234"),p("Real published dataset"))),
                      column(4,div(class="stat-card",
                                   icon("chart-line",style="font-size:30px;color:#E74C3C;"),
                                   br(),h4("DESeq2"),p("Industry-standard analysis"))),
                      column(4,div(class="stat-card",
                                   icon("upload",style="font-size:30px;color:#2ECC71;"),
                                   br(),h4("Your Data"),p("Upload your own counts")))
                    ),
                    hr(),
                    h4("How to use:"),
                    tags$ol(
                      tags$li("Run ",strong("analysis_final.R")," in RStudio first"),
                      tags$li("Click ",strong("Upload Data")," → Load Demo"),
                      tags$li("Click ",strong("Run DESeq2 Analysis")),
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
                    p("CSV: rows=genes, columns=samples, values=counts"),
                    fileInput("count_file","Choose Count Matrix (.csv)",
                              accept=c(".csv",".txt")),
                    tags$small(style="color:grey;",
                               "Gene names in first column, sample IDs as header")
                ),
                box(width=6, title="Upload Metadata", status="primary",
                    p("CSV: rows=samples. Must have a ",strong("condition"),
                      " column (2 groups)."),
                    fileInput("meta_file","Choose Metadata (.csv)",accept=".csv")
                )
              ),
              fluidRow(
                box(width=12,
                    title="OR Load Demo Dataset (GSE157234 — 48h post-LPS)",
                    status="success",
                    p("Reads pre-processed files saved by ",
                      strong("analysis_final.R"),":"),
                    tags$ul(
                      tags$li(code("data/count_matrix_48h_clean.csv")),
                      tags$li(code("data/metadata_48h_clean.csv"))
                    ),
                    p(icon("triangle-exclamation",style="color:orange;"),
                      " Run ",strong("analysis_final.R")," first if these don't exist."),
                    actionButton("load_demo","Load Demo Dataset (GSE157234)",
                                 icon=icon("play"),class="btn-success btn-lg"),
                    br(),br(),
                    uiOutput("data_status")
                )
              ),
              fluidRow(
                box(width=12, title="Run Analysis", status="warning",
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
                    p("Red=upregulated in Mutant (Tnf, Ccl5, Il6, Il12b…). ",
                      "Blue=downregulated (P2ry12, Sall1, Tmem119…). ",
                      "Adjust thresholds in sidebar."),
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
                    p("Mutant (red) should separate clearly from Control (green)",
                      " along PC1 at 48h — confirms strong hyperactivation effect."),
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
                      "Expect Tnf, Ccl5, Il6 cluster high in Mutant."),
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
                    h4("RNA-Seq Differential Expression Dashboard"),
                    p("Built with R + Shiny + DESeq2"),
                    hr(),
                    h4("Dataset"),
                    p(strong("GEO:"),"GSE157234"),
                    p(strong("Paper:"),"Shemer et al., Immunity 53, 1033–1049, 2020"),
                    p(strong("Comparison:"),"IL10R-Mutant vs Control microglia"),
                    p(strong("Timepoint:"),"48h post-LPS (peak hyperactivation)"),
                    p(strong("Key finding:"),
                      "Without IL-10 signalling, microglia hyperactivate",
                      "and produce toxic TNF, damaging neurons."),
                    hr(),
                    h4("Pipeline"),
                    tags$ul(
                      tags$li("Input: UTAP-normalized counts, GSE157234"),
                      tags$li("Subset: 48h post-LPS only, DKO excluded"),
                      tags$li("DESeq2 design: ~ condition"),
                      tags$li("Filter: genes ≥10 counts in ≥2 samples"),
                      tags$li("Significance: padj<0.05 AND |log2FC|>1"),
                      tags$li("VST transformation for PCA and heatmap")
                    )
                ),
                box(width=4, title="Developer", status="success",
                    h4("Md. Abrar Faiyaj"),
                    p("Bioinformatics Analyst"),
                    p(icon("github"),
                      a("github.com/mdabrarfaiyaj",
                        href="https://github.com/mdabrarfaiyaj",
                        target="_blank")),
                    hr(),
                    p(strong("Services:")),
                    tags$ul(
                      tags$li("Bulk RNA-Seq (DESeq2, edgeR)"),
                      tags$li("scRNA-Seq (Seurat)"),
                      tags$li("Custom Shiny dashboards"),
                      tags$li("Pathway enrichment (GSEA, GO, KEGG)")
                    )
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
      "Step 1: Run analysis_final.R in RStudio\n",
      "Step 2: Click 'Load Demo Dataset' here\n",
      "Step 3: Click 'Run DESeq2 Analysis'\n",
      "Step 4: Explore the plot tabs"
    )
  )
  
  # ---- Load Demo: reads CSV files saved by analysis_final.R ----
  observeEvent(input$load_demo, {
    withProgress(message="Loading pre-processed demo data...", {
      tryCatch({
        count_file <- "data/count_matrix_48h_clean.csv"
        meta_file  <- "data/metadata_48h_clean.csv"
        
        if (!file.exists(count_file))
          stop(paste("File not found:", count_file,
                     "\nPlease run analysis_final.R first."))
        if (!file.exists(meta_file))
          stop(paste("File not found:", meta_file,
                     "\nPlease run analysis_final.R first."))
        
        incProgress(0.3, "Reading count matrix...")
        mat <- read.csv(count_file, row.names=1, check.names=FALSE)
        mat <- round(as.matrix(mat))
        storage.mode(mat) <- "integer"
        mat[mat < 0] <- 0
        
        incProgress(0.6, "Reading metadata...")
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
          "✓ Demo data loaded!\n",
          "  Genes:   ", nrow(mat),"\n",
          "  Samples: ", ncol(mat),"\n",
          "  Control: ", n_ctrl," samples (48h post-LPS)\n",
          "  Mutant:  ", n_mut, " samples (48h post-LPS)\n\n",
          "→ Click 'Run DESeq2 Analysis' to proceed!"
        )
        showNotification("Demo data loaded!", type="message")
        incProgress(1.0)
        
      }, error=function(e) {
        rv$log_text <- paste0("ERROR: ",e$message,
                              "\n\nFix: Run analysis_final.R first.")
        showNotification(paste("Error:",e$message),
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
      mat <- round(as.matrix(df))
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
  
  # ---- Run DESeq2 ----
  observeEvent(input$run_analysis, {
    req(rv$counts, rv$metadata)
    withProgress(message="Running DESeq2...", value=0, {
      tryCatch({
        shared     <- intersect(colnames(rv$counts), rownames(rv$metadata))
        counts_sub <- rv$counts[, shared]
        meta_sub   <- rv$metadata[shared, , drop=FALSE]
        
        # Remove NA rows
        na_rows    <- rowSums(is.na(counts_sub)) > 0
        counts_sub <- counts_sub[!na_rows, ]
        
        incProgress(0.1, "Creating DESeq2 object...")
        dds  <- DESeqDataSetFromMatrix(countData=counts_sub,
                                       colData=meta_sub,
                                       design=~condition)
        keep <- rowSums(counts(dds) >= 10) >= 2
        dds  <- dds[keep, ]
        
        incProgress(0.3, "Running DESeq2 (1-3 min)...")
        dds <- DESeq(dds)
        
        incProgress(0.6, "Extracting results...")
        lvls  <- levels(meta_sub$condition)
        res   <- results(dds, contrast=c("condition",lvls[2],lvls[1]),
                         alpha=0.05)
        res_df       <- as.data.frame(res)
        res_df$gene  <- rownames(res_df)
        
        res_df$significance <- "Not Significant"
        res_df$significance[!is.na(res_df$padj) &
                              res_df$padj < input$padj_cutoff &
                              res_df$log2FoldChange > input$lfc_cutoff]  <- "Upregulated"
        res_df$significance[!is.na(res_df$padj) &
                              res_df$padj < input$padj_cutoff &
                              res_df$log2FoldChange < -input$lfc_cutoff] <- "Downregulated"
        res_df$significance <- factor(res_df$significance,
                                      levels=c("Upregulated","Downregulated","Not Significant"))
        
        incProgress(0.8, "VST transformation...")
        vsd <- vst(dds, blind=FALSE)
        
        rv$dds=dds; rv$res_df=res_df; rv$vsd=vsd; rv$ready=TRUE
        
        n_up   <- sum(res_df$significance=="Upregulated",  na.rm=TRUE)
        n_down <- sum(res_df$significance=="Downregulated",na.rm=TRUE)
        
        rv$log_text <- paste0(
          "✓ DESeq2 complete!\n",
          "  Comparison:    ",lvls[2]," vs ",lvls[1],"\n",
          "  Genes tested:  ",nrow(res_df),"\n",
          "  Upregulated:   ",n_up, " (padj<0.05, LFC>1)\n",
          "  Downregulated: ",n_down," (padj<0.05, LFC< -1)\n\n",
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
    n_up  <- sum(rv$res_df$significance=="Upregulated",  na.rm=TRUE)
    n_down<- sum(rv$res_df$significance=="Downregulated",na.rm=TRUE)
    fluidRow(
      infoBox("Upregulated",  n_up,           icon=icon("arrow-up"),  color="red"),
      infoBox("Downregulated",n_down,          icon=icon("arrow-down"),color="blue"),
      infoBox("Genes Tested", nrow(rv$res_df),icon=icon("dna"),       color="green")
    )
  })
  
  output$n_up  <- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Upregulated",  na.rm=TRUE) else "—"
    valueBox(n,"Upregulated",  icon=icon("arrow-up"),  color="red")
  })
  output$n_down<- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Downregulated",na.rm=TRUE) else "—"
    valueBox(n,"Downregulated",icon=icon("arrow-down"),color="blue")
  })
  output$n_ns  <- renderValueBox({
    n <- if(rv$ready) sum(rv$res_df$significance=="Not Significant",na.rm=TRUE) else "—"
    valueBox(n,"Not Significant",icon=icon("minus"),   color="gray")
  })
  
  # ---- Volcano ----
  make_volcano <- reactive({
    req(rv$ready)
    df <- rv$res_df |>
      filter(!is.na(padj),!is.na(log2FoldChange)) |>
      mutate(sig=case_when(
        padj<input$padj_cutoff & log2FoldChange> input$lfc_cutoff ~"Upregulated",
        padj<input$padj_cutoff & log2FoldChange< -input$lfc_cutoff~"Downregulated",
        TRUE~"Not Significant"),
        neg_log10_padj=-log10(padj+1e-300))
    top_lab <- df |> filter(sig!="Not Significant") |>
      arrange(padj) |> head(input$top_n_genes)
    
    ggplot(df,aes(x=log2FoldChange,y=neg_log10_padj,color=sig,text=gene))+
      geom_point(alpha=0.6,size=1.5)+
      scale_color_manual(values=c("Upregulated"="#E74C3C",
                                  "Downregulated"="#3498DB",
                                  "Not Significant"="grey70"))+
      geom_vline(xintercept=c(-input$lfc_cutoff,input$lfc_cutoff),
                 linetype="dashed",color="grey40")+
      geom_hline(yintercept=-log10(input$padj_cutoff),
                 linetype="dashed",color="grey40")+
      geom_text_repel(data=top_lab,aes(label=gene),
                      size=2.8,max.overlaps=25,color="black")+
      labs(title="Volcano Plot: IL10R-Mutant vs Control (48h post-LPS)",
           subtitle=paste0("padj<",input$padj_cutoff,
                           "  |log2FC|>",input$lfc_cutoff),
           x="Log2 Fold Change (Mutant/Control)",
           y="-Log10 Adjusted P-value",color="Regulation")+
      theme_bw(base_size=13)+
      theme(plot.title=element_text(face="bold"),
            panel.grid.minor=element_blank())
  })
  
  output$volcano_plot <- renderPlotly({
    req(rv$ready); ggplotly(make_volcano(),tooltip=c("text","x","y"))
  })
  output$dl_volcano_png <- downloadHandler(
    filename="volcano_plot.png",
    content=function(f) ggsave(f,make_volcano(),width=10,height=7,dpi=300))
  output$dl_volcano_pdf <- downloadHandler(
    filename="volcano_plot.pdf",
    content=function(f) ggsave(f,make_volcano(),width=10,height=7))
  
  # ---- PCA ----
  make_pca <- reactive({
    req(rv$ready)
    pca_data <- plotPCA(rv$vsd,intgroup="condition",returnData=TRUE)
    pct_var  <- round(100*attr(pca_data,"percentVar"),1)
    ggplot(pca_data,aes(x=PC1,y=PC2,color=condition,
                        shape=condition,text=name))+
      geom_point(size=4,alpha=0.85)+
      scale_color_manual(values=c("Control"="#2ECC71","Mutant"="#E74C3C"))+
      scale_shape_manual(values=c("Control"=16,"Mutant"=17))+
      stat_ellipse(aes(group=condition),linetype="dashed",level=0.8)+
      geom_text_repel(aes(label=name),size=2.5,max.overlaps=20)+
      labs(title="PCA: IL10R-Mutant vs Control (48h post-LPS)",
           x=paste0("PC1: ",pct_var[1],"% variance"),
           y=paste0("PC2: ",pct_var[2],"% variance"))+
      theme_bw(base_size=13)+
      theme(plot.title=element_text(face="bold"),
            panel.grid.minor=element_blank())
  })
  
  output$pca_plot <- renderPlotly({
    req(rv$ready); ggplotly(make_pca(),tooltip=c("text","x","y"))
  })
  output$dl_pca_png <- downloadHandler(
    filename="pca_plot.png",
    content=function(f) ggsave(f,make_pca(),width=9,height=6,dpi=300))
  output$dl_pca_pdf <- downloadHandler(
    filename="pca_plot.pdf",
    content=function(f) ggsave(f,make_pca(),width=9,height=6))
  
  # ---- Heatmap ----
  make_heatmap <- reactive({
    req(rv$ready)
    top_genes <- rv$res_df |>
      filter(!is.na(padj),significance!="Not Significant") |>
      arrange(padj) |> head(input$top_n_heat) |> pull(gene)
    if (length(top_genes)<5)
      top_genes <- rv$res_df |> filter(!is.na(pvalue)) |>
      arrange(pvalue) |> head(input$top_n_heat) |> pull(gene)
    
    heat_mat        <- assay(rv$vsd)[top_genes,]
    heat_mat_scaled <- t(scale(t(heat_mat)))
    heat_mat_scaled[heat_mat_scaled> 3] <-  3
    heat_mat_scaled[heat_mat_scaled< -3]<- -3
    
    lvls      <- levels(rv$metadata$condition)
    col_annot <- data.frame(
      Condition=as.character(rv$metadata[colnames(heat_mat),"condition"]),
      row.names=colnames(heat_mat))
    pal_named    <- setNames(c("#2ECC71","#E74C3C")[seq_along(lvls)],lvls)
    annot_colors <- list(Condition=pal_named)
    
    pheatmap(heat_mat_scaled,
             annotation_col=col_annot, annotation_colors=annot_colors,
             color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
             cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE,
             fontsize_row=max(5,9-input$top_n_heat/20),
             main=paste0("Top ",input$top_n_heat,
                         " DEGs | Mutant vs Control | 48h"),
             border_color=NA, silent=TRUE)
  })
  
  output$heatmap_plot <- renderPlot({
    req(rv$ready)
    grid::grid.newpage(); grid::grid.draw(make_heatmap()$gtable)
  })
  output$dl_heat_png <- downloadHandler(
    filename="heatmap.png",
    content=function(f){
      png(f,width=2400,height=3200,res=300)
      grid::grid.draw(make_heatmap()$gtable); dev.off()})
  output$dl_heat_pdf <- downloadHandler(
    filename="heatmap.pdf",
    content=function(f){
      pdf(f,width=10,height=13)
      grid::grid.draw(make_heatmap()$gtable); dev.off()})
  
  # ---- Results Table ----
  filtered_res <- reactive({
    req(rv$ready)
    df <- rv$res_df |> filter(!is.na(padj)) |>
      mutate(across(where(is.numeric),~round(.,4)))
    switch(input$table_filter,
           "sig" =df|>filter(significance!="Not Significant"),
           "up"  =df|>filter(significance=="Upregulated"),
           "down"=df|>filter(significance=="Downregulated"),
           "all" =df)
  })
  
  output$results_table <- renderDT({
    req(rv$ready)
    datatable(
      filtered_res()|>select(gene,baseMean,log2FoldChange,
                             pvalue,padj,significance),
      options=list(pageLength=25,scrollX=TRUE,
                   order=list(list(4,"asc"))),
      rownames=FALSE
    )|>formatStyle("significance",
                   backgroundColor=styleEqual(
                     c("Upregulated","Downregulated","Not Significant"),
                     c("#FADBD8","#D6EAF8","#F8F9F9")))
  })
  
  output$dl_results_csv <- downloadHandler(
    filename="DESeq2_results.csv",
    content=function(f) write.csv(filtered_res(),f,row.names=FALSE))
}

shinyApp(ui=ui, server=server)