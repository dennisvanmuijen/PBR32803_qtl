library(shiny)
shinyUI(
  navbarPage("", inverse = TRUE, theme = shinytheme("cerulean"),
    tabPanel(div(h4("Interval Mapping")),
          sidebarLayout(
           sidebarPanel(
             # uiOutput("estmap"),
             # br(),

             fileInput('file1', 'Choose file to upload',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    )
                    ),
             uiOutput("poptype"),
             br(),
             uiOutput("pheno"),      
             br(),
             h4("Select chromosome"),
             uiOutput("chromSelect")
             ),
            mainPanel(
              tabPanel("Interval mapping",
                  tabsetPanel("IM",
                        tabPanel("Input summary",
                          verbatimTextOutput(outputId = "inputSummary")
                          ),
                        tabPanel("LOD Profile",
                                 helpText("Interval mapping - white dashed = 5% signficance level, orange dashed = 1% significance level"),
                           ggiraphOutput("distPlot",width="100%", height="800px")
                        ),
                        tabPanel("Fit qtl",
                                 verbatimTextOutput(outputId = "IMsummary")
                        )
                       )
              )
            )
          )
      ),
    # tabPanel(div(h4("Data exploration")),
    #          tabsetPanel(
    #           tabPanel("Genetic Map",
    #            wellPanel(
    #            ggiraphOutput('raw_plot', width = "100%", height = "600px")
    #            )
    #            )
    #          )
    # ),
    # tabPanel(div(h4("Data exploration")),
    #          tabsetPanel(
    #            tabPanel("Genetic Map",
    # sidebarPanel(print("Hello World")),
    # wellPanel(
    #   ggiraphOutput("distPlot",width="100%", height="400px")
    # )
    #            )
    #          )
    # ),
    tabPanel(div(h4("About")),
             mainPanel(
               includeMarkdown("about.Rmd")
             )
    )
  )
)