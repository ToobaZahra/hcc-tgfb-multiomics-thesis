# app.R
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(survival)
library(survminer)

# Load data
integrated <- read.csv("../../03_results/multiomics_integrated.csv")
cox <- read.csv("../../03_results/cox_nexus_genes.csv")
variants <- read.csv("../../03_results/nexus_variants_summary.csv")
risk_df <- read.csv("../../03_results/glmnet_risk_scores.csv")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "TGF-β HCC Multi-Omics"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Multi-Omics Table", tabName = "table", icon = icon("table")),
      menuItem("Survival", tabName = "survival", icon = icon("heartbeat")),
      menuItem("Variants", tabName = "variants", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(
      # Overview
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
      # Multi-omics table
      tabItem(tabName = "table",
              fluidRow(
                box(title = "Integrated Multi-Omics Summary", width = 12,
                    DTOutput("omics_table")
                )
              )
      ),
      # Survival
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
      # Variants
      tabItem(tabName = "variants",
              fluidRow(
                box(title = "Variant Summary", width = 12,
                    DTOutput("variant_table")
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
      labs(title = "Risk Score Distribution",
           x = "Risk Score", y = "Count") +
      theme_bw()
  })
  
  output$variant_table <- renderDT({
    datatable(variants, options = list(pageLength = 15, scrollX = TRUE))
  })
}

# Run
shinyApp(ui = ui, server = server)