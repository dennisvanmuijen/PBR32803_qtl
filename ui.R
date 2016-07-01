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
                                 helpText("Interval mapping - white dashed = 5% genome-wide signficance level, orange dashed = 1% genome-wide significance level"),
                                 verbatimTextOutput(outputId = "lodthr1"),
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
    tabPanel(div(h4("Composite Interval Mapping")),
             sidebarLayout(
               sidebarPanel(
                 uiOutput("poptype2"),
                 br(),
                 uiOutput("pheno2"),      
                 br(),
                 h4("Select chromosome"),
                 uiOutput("chromSelect2")
               ),
               mainPanel(
                 tabPanel("Composite interval mapping",
                          tabsetPanel("CIM",
                                      tabPanel("LOD Profile CIM",
                                               helpText("Composite Interval mapping - white dashed = 5% signficance level, orange dashed = 1% significance level"),
                                               verbatimTextOutput(outputId = "lodthr2"),
                                               ggiraphOutput("distPlot2",width="100%", height="800px")
                                      ),
                                      tabPanel("Fitted qtl",
                                               verbatimTextOutput(outputId = "CIMsummary")
                                      )
                          )
                 )
               )
             )
    ),
    tabPanel(div(h4("Export Results")),
             sidebarLayout(
               sidebarPanel(
                 helpText("If there issues opening in Excel, select the semicolon as delimiter"),
                 selectInput("sep_type", "Select delimiter",
                             multiple = FALSE, choices = c(",",";")),
                 htmlOutput("save")
               ),
               mainPanel(
                 helpText("Select a dataset to export and click download")
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