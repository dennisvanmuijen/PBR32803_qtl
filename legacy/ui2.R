library(shiny)
shinyUI(pageWithSidebar(  
  headerPanel("Spatially Balanced Sampling Tool"),
  sidebarPanel(
    conditionalPanel(condition="input.conditionedPanels == 'Tab1'",       
                     helpText("TAB 1")
    ),    
    conditionalPanel(condition="input.conditionedPanels == 'Tab2'", 
                     helpText("TAB 2 SELECTED")
    )
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Tab1",
               p(paste("Tab 1 text")                 
               )
      ),
      tabPanel("Tab2",  helpText("Map of input polygons"),
               p(paste("polygons by strata")
               )
      ), 
      id = "conditionedPanels"                
    )
  )  
))