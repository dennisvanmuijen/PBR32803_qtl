library(shiny)
library(ggplot2)
library(ggiraph)

ui <- shinyUI(fluidPage(
  sidebarPanel(print("Hello World")),
  wellPanel(
    ggiraphOutput("distPlot",width="100%", height="400px")
  )
))

server <- shinyServer(function(input, output) {
  
  output$distPlot <- renderggiraph({
    #Dit moet ff gekoppeld worden aan de input van de crossfile lader
    cross <- read.cross("data1.geg.csv", dir = "",
                        format = "csv", genotypes = c("A","H","B"),
                        alleles = c("A","B"))
    plotdata <- scanone(cross, pheno.col = 1)
    
    p <- ggplot(plotdata, aes(x=pos, y=lod, tooltip = lod))+
      geom_line()+geom_rug(sides = "b")+
      geom_point_interactive(color="orange",size=2)+
      theme_dark()
    return(ggiraph(code = {print(p)},zoom_max = 2))
  })
})

shinyApp(ui = ui, server = server)