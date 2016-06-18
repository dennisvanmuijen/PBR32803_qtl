library(shiny)
shinyServer(function(input, output, session) {
  ##################################
  ###### Upload data
  ##################################
  output$poptype <- renderUI({
      selectInput("ngen", label = "F Generation",
                  choices = 2:10,
                  selected = 2)
  })
  
  output$chromSelect <- renderUI({
    #     if(input$mapactivator == 0){
    selectizeInput(inputId="chromSelect", 
                   label = "Select Chromosome", 
                   choices = names(geno()$geno),
                   multiple = TRUE,
                   selected =  names(geno()$geno)
    )
    #     } else {
    #       selectizeInput(inputId="chromSelect",
    #                      label = "Select Chromosome",
    #                      choices = names(mstresult()$geno)
    #       )
    #     }
  })  
  
  
  output$pheno <- renderUI({
    selectInput("phenosel", label = "Select Phenotype",
                choices = geno()$pheno %>% names)
  })
  
  # output$estmap <- renderUI({
  #   checkboxInput("mapest", label = "Estimate map",
  #               value = TRUE)
  # })
  # 
  
    geno <- reactive({
    inFile <- input$file1
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
      # need(input$poptype != "", "Upload a cross file to begin")
    )
      cross <- read.cross(
        format = "csv",
        file = inFile$datapath,
        genotypes = c("A","H","B"),
        alleles = c("A","B"),
        estimate.map = TRUE,
        F.gen = input$ngen %>% as.numeric,
        BC.gen = 0
      )
      cross <- jittermap(cross)
      cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
      return(cross)
  })
  
  output$inputSummary <- renderPrint({
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
    )
    summary(geno())
  })
  
  
###############################
###### Interval mapping
IMapping <- reactive({
  out.s1perm <- scanone(geno(), pheno.col = input$phenosel, n.perm = 50, n.cluster = 4, method = "hk")
  out.s1  <- scanone(geno(), pheno.col = input$phenosel, method = "hk")###############################
  thrs <- summary(out.s1perm , alpha=c(0.05, 0.01))
  # QTL detect at alpha = 0.05
  res <- summary(out.s1, perms=out.s1perm, alpha=0.05)
  if(dim(res)[1] != 0){
  mqtl <- makeqtl(geno(),res[,1], res[,2], what="prob")
  QTLnames <- unlist(mqtl$altname)
  form <- paste(QTLnames, collapse="+")
  form <- paste("y~",form, sep="")
  fit <- fitqtl(geno(), pheno.col = input$phenosel, qtl = mqtl, method = "hk", model = "normal",
                formula = form, get.ests = T)
  estimate.out <- summary(fit)
  mylist <- list()
  mylist[[1]] <- out.s1
  mylist[[2]] <- thrs
  mylist[[3]] <- estimate.out
  } 
  if(dim(res)[1] == 0){
  mylist <- list()
  mylist[[1]] <- 0
  mylist[[2]] <- 0
  mylist[[3]] <- 0
}
  return(mylist)
})

output$IMsummary <- renderPrint({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(IMapping()[[1]] != 0, "No QTL detected")
  )
  IMapping()[[3]]
})

###############################
###### Data exploration
###############################
	
## Plot of raw data
output$raw_plot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(IMapping()[[1]] != 0, "No QTL detected")
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
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(IMapping()[[1]] != 0, "No QTL detected")
  )
  plotdata <- IMapping()[[1]]
  thr <- IMapping()[[2]]
    p <- ggplot(plotdata[which(plotdata$chr%in%input$chromSelect),], aes(x=pos, y=lod, tooltip = lod))+
      geom_line()+geom_rug(sides = "b")+
      facet_wrap(~chr, nrow=3)+
      geom_point_interactive(color="orange",size=0.5)+
      theme_dark() +
      geom_hline(yintercept = thr[1], lwd = 0.5, lty = 2, col = "white") +
      geom_hline(yintercept = thr[2], lwd = 0.5, lty = 2, col = "orange")
    
    return(ggiraph(code = {print(p)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;"))
  })
  
})