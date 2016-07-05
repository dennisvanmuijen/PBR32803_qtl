library(shiny)
library(ggplot2)
library(ggiraph)

ui <- shinyUI(fluidPage(
  column(3,
         wellPanel(
           h4("dataselection")
         ),
         
         wellPanel(
           h4("Select chromosome"),
           uiOutput("chromSelect")
         )
  ),
  column(9,
         tabsetPanel(type = "tabs",
                     tabPanel("Input Summary"
                     ),
                     tabPanel("QTL distance Plot",
                              ggiraphOutput("distPlot",width="100%", height="800px")
                     ),
                     tabPanel("foobar"
                     )
         )
  )
))

server <- shinyServer(function(input, output) {
  
  geno <- reactive({
    cross <- read.cross("data/cross.csv", dir = "",
                        format = "csv", genotypes = c("A","H","B"),
                        alleles = c("A","B"))
    return(cross)
  })
  
  plotdata <- reactive({
    plotdata <- scanone(geno(), pheno.col = 1)
    return(plotdata)
  })
  
  output$chromSelect <- renderUI({
    #     if(input$mapactivator == 0){
    selectizeInput(inputId="chromSelect", 
                   label = "Select Chromosome", 
                   choices = names(geno()$geno),
                   multiple = TRUE
    )
    #     } else {
    #       selectizeInput(inputId="chromSelect",
    #                      label = "Select Chromosome",
    #                      choices = names(mstresult()$geno)
    #       )
    #     }
  })  
  
  output$distPlot <- renderggiraph({
    if(is.null(input$chromSelect))
      return(NULL)
    #     validate(
    #       need(input$file1 != "", "Upload a cross file to begin")
    #     )
    p <- ggplot(plotdata()[which(plotdata()$chr%in%input$chromSelect),], aes(x=pos, y=lod, tooltip = lod))+
      geom_line()+geom_rug(sides = "b")+
      facet_wrap(~chr, nrow=3)+
      geom_point_interactive(color="orange",size=0.5)+
      theme_dark()
    
    return(ggiraph(code = {print(p)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;"))
    
  })
})

shinyApp(ui = ui, server = server)