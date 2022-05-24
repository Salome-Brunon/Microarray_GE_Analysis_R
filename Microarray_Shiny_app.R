######WELCOME TO STUDENT B118056 APP TO VISUALISE AFFYMETRIX DATA ANALYSIS#################

#load required packages
library(shiny)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(pheatmap)
library(shinyjs)
library(ggplot2)
library(ggrepel)
library(DT)
graphics.off()

#load required data
load("expression.Rdata")
load("limma_genes.Rdata")
load("enrichment_camera.Rdata")
load("enrichment_romer.Rdata")
load("enrichment_roast.Rdata")
experiment <- data.frame(limma_genes[,-1], row.names = limma_genes[,1])

ui <- fluidPage(
    tabsetPanel(
        tabPanel("Volcano Plot", fluid = TRUE,
                 
                 sidebarLayout(
                     #first panel - volcano plot
                     sidebarPanel(
                         #slider for LogFc Threshold
                         sliderInput("LogFC_slider",
                                     "LogFc Threshold:",
                                     #minimum value
                                     min = 0,
                                     #maximum value
                                     max = 6,
                                     #value input when starting application
                                     value = 3),
                         sliderInput("pval_slider",
                                     "P value threshold",
                                     min = 0.001,
                                     max = 0.05,
                                     value = 0.05),
                         sliderInput("AveExpr_data",
                                     "Minimum Average expression",
                                     min = 0,
                                     max = 15,
                                     value = 0),
                         #if genes very close together, their names will overlap
                         #avoid this with overlap option
                         sliderInput("overlap",
                                     "Overlapping gene names",
                                     min = 10,
                                     max = 50,
                                     value = 10),
                         checkboxInput("names", "Show Gene Names", TRUE),
                         checkboxInput("lines", "Show Threshold Lines", TRUE)
                     ),
                     mainPanel(
                         plotOutput("VolcanoPlot")
                     )
                 )
        ),
        tabPanel("Heatmap", fluid = TRUE,
                 # second panel - Heatmap of gene expression data
                 sidebarLayout(
                     sidebarPanel(
                         sliderInput("font_row",
                                     "Font size row:",
                                     min = 6,
                                     max = 14,
                                     value = 10),
                         sliderInput("font_col",
                                     "Font size col:",
                                     min = 6,
                                     max = 14,
                                     value = 10),
                         checkboxInput("srownames", "Show Row Names", TRUE),
                         checkboxInput("logtansform", "Log transform values", TRUE),
                         radioButtons("norm", "Scale by", choices=c("none","row","column"))
                         
                     ),
                     mainPanel(
                         plotOutput("HeatMap")
                     )
                 )
        ),
        tabPanel("Functional Enrichment", fluid = TRUE,
                 #third panel - functional enrichment analysis
                 sidebarLayout(
                     sidebarPanel(
                         helpText("Carry out functional enrichment analysis of choice"),
                         selectInput("user_input", 
                                     label = "Choose a test to conduct",
                                     choices = c("Camera (default)", 
                                                 "Romer",
                                                 "Roast"),
                                     selected = "Camera (default)"),
                     ),
                     mainPanel(
                         DT::dataTableOutput("FuncEnrich")
                     )
                 )
        )
    )
)    

server <- function(input, output,session) {
    output$VolcanoPlot <- renderPlot({
        #code adjusted from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
        limma<-limma_genes[limma_genes$AveExpr > input$AveExpr_data,]
        #The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
        # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
        #classify all genes as non-differentially expressed
        limma$diffexpressed <- "NO"
        # Up Differentially Expressed genes (adj.p.val <0.05 or slider input value   AND Fold change > 3 or slider)
        limma$diffexpressed[limma$logFC > input$LogFC_slider & limma$adj.P.Val < input$pval_slider] <- "UP"
        # Down Differentially Expressed genes (adj.p.val <0.05 or slider input value   AND Fold change < than 3 or slider)
        limma$diffexpressed[limma$logFC < -input$LogFC_slider & limma$adj.P.Val < input$pval_slider] <- "DOWN"
        
        p <- ggplot(data=limma, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()
        # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
        p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
        #annotating genes fold change going up or down
        colours_annotation <- c("blue", "red", "black")
        names(colours_annotation) <- c("DOWN", "UP", "NO")
        p3 <- p2 + scale_colour_manual(values = colours_annotation)
        # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
        limma$label <- NA
        limma$label[limma$diffexpressed != "NO"] <- limma$Symbol[limma$diffexpressed != "NO"]
        ggplot(data=limma, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label='')) + 
            geom_point() + 
            theme_minimal() +
            geom_text_repel() +
            scale_color_manual(values=c("blue", "black", "red")) +
            geom_vline(xintercept=c(-input$LogFC_slider, input$LogFC_slider), col="red") +
            geom_hline(yintercept=-log10(input$pval_slider), col="red")
    },execOnResize = F)
    output$HeatMap <- renderPlot({
        #log-transformation of data if user checked the option
        if(input$logtansform){
            expression <- log2(expression + 1)
        }
        pheatmap(expression,
                 fontsize_row = input$font_row,
                 fontsize_col=input$font_col,
                 show_rownames=input$srownames,
                 scale=input$norm)
    }, execOnResize = F)
    #display dataframe of functional enrichment data depending on user choice of analysis
    output$FuncEnrich <- DT::renderDataTable({
        if(input$user_input == "Camera (default)"){
            enrichment_camera
        } else if (input$user_input == "Romer") {
            enrichment_romer
        } else if (input$user_input == "Roast") {
            enrichment_roast
        }
    })
}
shinyApp(ui = ui, server = server)
shinyApp(ui = ui, server = server)
shinyApp(ui = ui, server = server)

            