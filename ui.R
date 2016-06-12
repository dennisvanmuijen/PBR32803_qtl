library(shiny)
shinyUI(
  navbarPage("", inverse = TRUE, theme = shinytheme("cerulean"),
    tabPanel(div(h4("Upload data")),
          sidebarLayout(
           sidebarPanel(
             uiOutput("poptype"),
                 fileInput('file1', 'Choose file to upload',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    )
                    )),
            mainPanel(
              verbatimTextOutput(outputId = "inputSummary")
                    )
          )
      ),
    tabPanel(div(h4("Data exploration")),
             tabsetPanel(
              tabPanel("Genetic Map",
               wellPanel(
               ggiraphOutput('raw_plot', width = "100%", height = "600px")
               )
               )
             )
    ),
       tabPanel(div(h4("About")),
             mainPanel(
               includeMarkdown("about.Rmd")
             )
    )
  )
)