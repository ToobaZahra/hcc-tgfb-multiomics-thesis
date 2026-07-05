# app.R
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(survival)
library(survminer)
library(SummarizedExperiment)

# Load data
integrated <- read.csv("../../03_results/dge/multiomics_integrated.csv")
cox <- read.csv("../../03_results/survival/cox_univariate.csv")
variants <- read.csv("../../03_results/variants/nexus_variants_summary.csv")
risk_df <- read.csv("../../03_results/survival/glmnet_risk_scores.csv")
vsd_app <- readRDS("../../data/raw/tcga_lihc_vsd.rds")
res_app <- readRDS("../../03_results/dge/dge_tcga_tumor_vs_normal.rds")
beta_app <- readRDS("../../data/raw/tcga_lihc_meth_filtered.rds")
meth_se_app <- readRDS("../../data/raw/tcga_lihc_meth450_se.rds")
probe_annot_app <- as.data.frame(rowData(meth_se_app))
probe_annot_app$probe_id <- rownames(probe_annot_app)

nexus_genes <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7",
                 "MAPK8","MAPK9","MAPK10","DUSP1","DUSP4","DUSP10",
                 "MYC","CDKN1A","CDKN2B","SNAI1","TNF","IL6","IL10","IL37",
                 "HDAC11","NPC1","CCDC110","TGFBRAP1","KLF4")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "TGF-β HCC Multi-Omics"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Multi-Omics Table", tabName = "table", icon = icon("table")),
      menuItem("Survival", tabName = "survival", icon = icon("heartbeat")),
      menuItem("Variants", tabName = "variants", icon = icon("dna")),
      menuItem("Expression", tabName = "expression", icon = icon("chart-bar")),
      menuItem("Methylation", tabName = "methylation", icon = icon("flask"))
    ),
    selectInput("gene", "Select Gene:", choices = nexus_genes)
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                valueBox(26, "TGF-β Nexus Genes", icon = icon("gene"), color = "blue"),
                valueBox(16, "Replicated Genes", icon = icon("check"), color = "green"),
                valueBox("HDAC11", "Top Hit (4 omics layers)", icon = icon("star"), color = "red")
              ),
              fluidRow(
                box(title = "Study Summary", width = 12,
                    p("Multi-omics analysis of TGF-β pathway dysregulation in HCC."),
                    p("Primary cohort: TCGA-LIHC (n=462). Replication: GSE14520, GSE36376, GSE76427 (n=1,088)."),
                    p("Analyses: DESeq2, limma, methylation 450K, Cox regression, glmnet signature, variant annotation.")
                )
              )
      ),
      tabItem(tabName = "table",
              fluidRow(
                box(title = "Integrated Multi-Omics Summary", width = 12,
                    DTOutput("omics_table")
                )
              )
      ),
      tabItem(tabName = "survival",
              fluidRow(
                box(title = "Cox Regression Results", width = 6,
                    DTOutput("cox_table")
                ),
                box(title = "Risk Score Distribution", width = 6,
                    plotOutput("risk_plot")
                )
              )
      ),
      tabItem(tabName = "variants",
              fluidRow(
                box(title = "Variant Summary", width = 12,
                    DTOutput("variant_table")
                )
              )
      ),
      tabItem(tabName = "expression",
              fluidRow(
                box(title = "Gene Expression — Tumor vs Normal", width = 12,
                    plotOutput("expr_plot")
                )
              )
      ),
      tabItem(tabName = "methylation",
              fluidRow(
                box(title = "Gene Methylation — Tumor vs Normal", width = 12,
                    plotOutput("meth_plot")
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output) {
  output$omics_table <- renderDT({
    datatable(integrated, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$cox_table <- renderDT({
    datatable(cox, options = list(pageLength = 15))
  })
  output$risk_plot <- renderPlot({
    ggplot(risk_df, aes(x = risk_score, fill = risk_group)) +
      geom_histogram(bins = 30, alpha = 0.7) +
      scale_fill_manual(values = c("Low" = "#2E7C4A", "High" = "#D9534F")) +
      labs(title = "Risk Score Distribution", x = "Risk Score", y = "Count") +
      theme_bw()
  })
  output$variant_table <- renderDT({
    datatable(variants, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$expr_plot <- renderPlot({
    gene <- input$gene
    ens <- res_app$ensembl_id[res_app$symbol == gene & !is.na(res_app$symbol)][1]
    row_match <- grep(paste0("^", ens), rownames(assay(vsd_app)))
    if (length(row_match) == 0) return(NULL)
    df <- data.frame(
      expression = assay(vsd_app)[row_match[1], ],
      sample_type = vsd_app$sample_type
    )
    ggplot(df, aes(sample_type, expression, fill = sample_type)) +
      geom_boxplot() +
      scale_fill_manual(values = c("Primary Tumor" = "#D9534F",
                                   "Solid Tissue Normal" = "#2E7C4A")) +
      labs(title = paste(gene, "— VST Expression"), x = "", y = "VST") +
      theme_bw() + theme(legend.position = "none")
  })
  output$meth_plot <- renderPlot({
    gene <- input$gene
    probes <- probe_annot_app$probe_id[probe_annot_app$gene == gene]
    probes <- probes[probes %in% rownames(beta_app)]
    if (length(probes) == 0) return(NULL)
    keep <- meth_se_app$shortLetterCode %in% c("TP", "NT")
    keep_idx <- which(keep)[1:ncol(beta_app)]
    meth_vals <- colMeans(beta_app[probes, , drop = FALSE], na.rm = TRUE)
    df <- data.frame(
      beta = meth_vals,
      sample_type = meth_se_app$shortLetterCode[keep_idx]
    )
    ggplot(df, aes(sample_type, beta, fill = sample_type)) +
      geom_boxplot() +
      scale_fill_manual(values = c("NT" = "#2E7C4A", "TP" = "#D9534F")) +
      labs(title = paste(gene, "— Methylation (Beta)"), x = "", y = "Beta value") +
      theme_bw() + theme(legend.position = "none")
  })
}

shinyApp(ui = ui, server = server)