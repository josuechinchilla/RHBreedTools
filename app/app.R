require(shiny)
require(viridis)
require(BIGr)
require(scales)
require(tidyverse)
require(openxlsx)

#### Helper Functions ####
format_percent <- function(x) {
  percent_format(accuracy = 0.1)(x)
}

# Load static reference panel
reference <- as.data.frame(
  read.table(
    "https://raw.githubusercontent.com/josuechinchilla/RHBreedTools/main/data/2010_reference.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  ) %>%
    dplyr::select(-c(1, 2)) %>%
    t()
)

reference_ids <- read.table(
  "https://raw.githubusercontent.com/josuechinchilla/RHBreedTools/main/data/ref_ids.txt",
  header = TRUE,
  sep = "\t"
)
ref_ids <- lapply(as.list(reference_ids), as.character)

#### UI ####
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .header-img {
        text-align: center;
        margin-bottom: 20px;
      }
      .header-img img {
        max-width: 80%;
        height: auto;
      }
    "))
  ),
  
  div(class = "header-img",
      img(src = "logos.png", alt = "Logos")
  ),
  
  titlePanel("Russian Honeybee (RHB) content estimation"),
  
  fluidRow(
    column(12,
           wellPanel(
             HTML('
  <ul>
    <li>This tool was developed by <strong>Breeding Insight</strong> in collaboration with the USDA Honey Bee Breeding, Genetics, and Physiology Research Lab.</li>
    <li>It estimates the proportion of <strong>Russian Honey Bee (RHB)</strong> ancestry in genotype samples using methods from 
      <a href="https://www.animalsciencepublications.org/publications/tas/articles/1/1/36" target="_blank">Funkhouser et al. (2017)</a>.
    </li>
    <li><strong>Input format:</strong> Upload a genotype matrix (.txt) with SNPs in rows and samples in columns. The first three columns must be <code>ID</code>, <code>ref</code>, and <code>alt</code>. Following columns contain genotype calls (0,1,2).</li>
    <li>Example format:
      <pre>
ID    ref  alt  Sample1  Sample2  Sample3 ...
SNP1  T    A    0        0        0
SNP2  G    C    0        0        0
SNP3  A    C    1        1        1
      </pre>
    </li>
    <li>Results can be downloaded as an Excel file.</li>
  </ul>
  <p><strong>For publication, please cite:</strong></p>
  <ul>
    <li><a href="https://academic.oup.com/tas/article/1/1/36/4636602" target="_blank">Funkhouser et al.</a> for original breedTools methods</li>
    <li><a href="https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.70067" target="_blank">Sandercock et al.</a> for methods expansion to polyploidy</li>
  </ul>
')
           )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("validation_file", "Upload Genotypes to test (.txt)", accept = ".txt"),
      actionButton("run", "Run Estimation"),
      br(), br(),
      downloadButton("download_results", "Download Excel Results")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Results Table", DT::dataTableOutput("preview")),
        tabPanel("Ancestry Plot", plotOutput("bar_plot"))
      ),
      verbatimTextOutput("status")
    )
  )
)

#### Server ####
server <- function(input, output, session) {
  result_data <- reactiveVal(NULL)
  result_filename <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    req(input$validation_file)
    
    output$status <- renderText("Running estimation...")
    
    # Read validation data
    validation <- as.data.frame(
      read.table(
        input$validation_file$datapath,
        header = TRUE,
        row.names = 1,
        sep = "\t"
      ) %>%
        dplyr::select(-c(1, 2)) %>%
        t()
    )
    
    # Allele frequency
    freq <- BIGr:::allele_freq_poly(reference, ref_ids, ploidy = 2)
    
    # Breed prediction
    prediction <- as.data.frame(BIGr:::solve_composition_poly(validation, freq, ploidy = 2)) %>%
      select(-R2) %>%
      rename(
        `RHB` = RHB,
        `non-RHB` = non.RHB
      )
    
    columns_to_select <- c("RHB", "non-RHB")
    
    pred_results <- prediction %>%
      as.data.frame() %>%
      rownames_to_column(var = "ID") %>%
      mutate(
        across(c("RHB", "non-RHB"), ~ round(.x * 100, 0))  # multiply by 100 and round 0 decimals
      )
    
    # Rename columns to include % in header but keep numeric values raw
    colnames(pred_results)[colnames(pred_results) == "RHB"] <- "RHB (%)"
    colnames(pred_results)[colnames(pred_results) == "non-RHB"] <- "non-RHB (%)"
    
    # Add predicted line column based on max value in RHB and non-RHB columns
    pred_results <- pred_results %>%
      mutate(
        `Predicted line` = columns_to_select[max.col(select(., all_of(c("RHB (%)", "non-RHB (%)"))), ties.method = "first")]
      )
    
    result_data(pred_results)
    
    # Save Excel file to temp
    date_str <- format(Sys.Date(), "%Y-%m-%d")
    filename <- paste0("RHB_estimation_", date_str, ".xlsx")
    temp_path <- file.path(tempdir(), filename)
    result_filename(temp_path)
    write.xlsx(pred_results, file = temp_path, rowNames = FALSE)
    
    # Table
    output$preview <- DT::renderDataTable({
      DT::datatable(pred_results, options = list(pageLength = 10))
    })
    
    # Plot
    pred_results_long <- prediction %>%
      as.data.frame() %>%
      rownames_to_column(var = "ID") %>%
      pivot_longer(
        cols = where(is.numeric), 
        names_to = "category", 
        values_to = "percent"
      )
    
    output$bar_plot <- renderPlot({
      ggplot(pred_results_long, aes(x = ID, y = percent, fill = category)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d(option = "D") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
        labs(
          x = "Individual ID",
          y = "Ancestry Proportion",
          fill = "Line"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    })
    
    output$status <- renderText("Estimation complete. File ready for download.")
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      basename(result_filename())
    },
    content = function(file) {
      file.copy(result_filename(), file)
    }
  )
}

shinyApp(ui, server)