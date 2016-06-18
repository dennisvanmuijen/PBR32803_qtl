library(shiny)
shinyServer(function(input, output, session) {
    ##################################
  ###### Upload data
  ##################################
  output$poptype <- renderUI({
      selectInput("poptype", label = "Cross type",
                  choices = c("F2","RIL","DH"),
                  selected = 1)
  })
  
    geno <- reactive({
    inFile <- input$file1
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
      # need(input$poptype != "", "Upload a cross file to begin")
    )
    if(input$poptype == "F2"){
      cross <- read.cross(
        format = "csv",
        file = inFile$datapath,
        genotypes = c("A","H","B"),
        alleles = c("A","B"),
        estimate.map = TRUE
      )
      return(cross)
    } 
    if(input$poptype == "DH"){
      cross <- read.cross(
        format = "csv",
        file = inFile$datapath,
        genotypes = c("A","B"),
        alleles = c("A","B"),
        estimate.map = TRUE
      )
      return(cross)
    } 
    if(input$poptype == "RIL"){
      cross <- read.cross(
        format = "csv",
        file = inFile$datapath,
        genotypes = c("A","H","B"),
        alleles = c("A","B"),
        BC.gen = 0,
        F.gen = 6,
        estimate.map = TRUE
      )
      cross <- cross %>% convert2riself()
      return(cross)
    } 
  })
  
  output$inputSummary <- renderPrint({
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
    )
    summary(geno())
  })
  
###############################
###### Data exploration
###############################
	
## Plot of raw data
output$raw_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  if(input$mapactivator == 0){
  mymap <- pull.map(geno(), as.table = T)
  mymap <- data.frame(marker = row.names(mymap), mymap)
  mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
    geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
    theme_dark(base_size = 12)  + ylab("cM") +
    xlab("Linkage Group") + scale_y_reverse() +
    geom_point_interactive(size=4.5, col="orange", shape = 95) 
  ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
  } else {
    mymap <- pull.map(mstresult(), as.table = T)
    mymap <- data.frame(marker = row.names(mymap), mymap)
    mapplot <- ggplot(aes(x = factor(chr), y = pos, tooltip=marker, data_id=marker), data = mymap)  +
      geom_path(data = mymap, aes(x = factor(chr), y = pos, group = chr),size=4.5, lineend="round", col = colors()[284]) +
      theme_dark(base_size = 12)  + ylab("cM") +
      xlab("Linkage Group") + scale_y_reverse() +
      geom_point_interactive(size=4.5, col="orange", shape = 95)
    ggiraph(code = {print(mapplot)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;")
  }
  })

#######################
####QTL visualisation##
#######################

output$distPlot <- renderggiraph({
  plotdata <- scanone(geno(), pheno.col = 1)
  
  p <- ggplot(plotdata, aes(x=pos, y=lod, tooltip = lod))+
    geom_line()+geom_rug(sides = "b")+
    geom_point_interactive(color="orange",size=0.1)+
    theme_dark() +
    facet_wrap(~chr)
  return(ggiraph(code = {print(p)},zoom_max = 2))
})



########################
#####Export results ####
########################

# Ui to save files
output$save <- renderUI({
  list(
    selectInput("datatype", label = "Select data to export",
                choices = c("Genotype file","Genetic map","Segregation distortion","Recombination fractions"),
                selected = 1),
    downloadButton('downloadData', 'Save')
  )
})

output$DownloadData <- downloadHandler(
  filename = function() {
    paste(input$Download, format(Sys.time(), "%a %b %d %X"), '.csv', sep='', col.names = NA)
  },
  content = function(con) {
    write.csv(data, con)
  }
)

# Download file handler
observeEvent(input$datatype, {
  if (input$datatype == 'Genotype file'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('Genotype_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          genodata <- mstresult() %>% pull.geno %>% as.data.frame
          row.names(genodata) <- mstresult()$pheno$RILs
          genodata[genodata == 1] <- "A"
          genodata[genodata == 2] <- "H"
          genodata[genodata == 3] <- "B"
          genodata <- cbind(Rils = rownames(genodata), genodata)
          row.names(genodata) <- NULL
          write.table(genodata, file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Genetic map'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('Genetic map_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          mymap <- pull.map(mstresult(), as.table = T)
          mymap <- data.frame(marker = row.names(mymap), mymap)
          mymap$bp <-lapply(mymap$marker %>% as.character %>% strsplit(split = "-", fixed = T),"[",1) %>% as.numeric
          mymap <- mymap %>% arrange(chr,pos)
          write.table(mymap, file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Segregation distortion'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('SegDist_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          segdist_data <- mstresult() %>% geno.table
          colnames(segdist_data)[3:4] <- c("A","B")
          segdist_data <- cbind(marker = rownames(segdist_data), segdist_data )
          map <- pull.map(mstresult(), as.table = T)
          map <- cbind(marker = row.names(map),map)
          result <- left_join(map[,c(1,3)], segdist_data, by = "marker") %>% arrange(chr,pos)
          write.table(result[,1:6], file, sep = ",", row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'Recombination fractions'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('RecFraction_', format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          rf_data <- mstresult()$rf %>% as.data.frame()
          rf_data <- cbind(marker = row.names(rf_data), rf_data )
          write.table(rf_data, file, sep = ",", row.names = FALSE)
        }
      )}
  }
})
})


