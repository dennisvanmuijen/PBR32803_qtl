library(shiny)
shinyServer(function(input, output, session) {
  ##################################
  ###### Upload data
  ##################################
  output$poptype <- renderUI({
      selectInput("ngen", label = "F Generation",
                  choices = 2:10,
                  selected = input$ngen2)
  })
  
  output$poptype2 <- renderUI({
    selectInput("ngen2", label = "F Generation",
                choices = 2:10,
                selected = input$ngen)
  })

  output$chromSelect2 <- renderUI({
    selectizeInput(inputId="chromSelect2", 
                   label = "Select Chromosome", 
                   choices = names(geno()$geno),
                   multiple = TRUE,
                   selected =  names(geno()$geno)
    )
  })  

  output$chromSelect <- renderUI({
    selectizeInput(inputId="chromSelect", 
                   label = "Select Chromosome", 
                   choices = names(geno()$geno),
                   multiple = TRUE,
                   selected =  names(geno()$geno)
    )
  })  
  
  output$pheno <- renderUI({
    selectInput("phenosel", label = "Select Phenotype",
                choices = geno()$pheno %>% names, selected = input$phenosel2
                )
  })
  
  output$pheno2 <- renderUI({
    selectInput("phenosel2", label = "Select Phenotype",
                choices = geno()$pheno %>% names, selected = input$phenosel
    )
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
  if(input$ngen > 5){
  # cross <- geno() %>% convert2riself()
    cross <- geno()
    cross <- cross %>% jittermap()
  cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
  } else {
    cross <- geno()
    cross <- cross %>% jittermap()
    cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
  }
  out.s1perm <- scanone(cross, pheno.col = input$phenosel, n.perm = 2500, n.cluster = 4, method = "hk")
  out.s1  <- scanone(cross, pheno.col = input$phenosel, method = "hk")###############################
  thrs <- summary(out.s1perm , alpha=c(0.05, 0.01))
  # QTL detect at alpha = 0.05
  res <- summary(out.s1, perms=out.s1perm, alpha=0.05)
  if(dim(res)[1] != 0){
  mqtl <- makeqtl(cross,res[,1], res[,2], what="prob")
  QTLnames <- unlist(mqtl$altname)
  form <- paste(QTLnames, collapse="+")
  form <- paste("y~",form, sep="")
  fit <- fitqtl(cross, pheno.col = input$phenosel, qtl = mqtl, method = "hk", model = "normal",
                formula = form, get.ests = T)
  estimate.out <- summary(fit)
  mylist <- list()
  mylist[[1]] <- out.s1
  mylist[[2]] <- thrs
  mylist[[3]] <- estimate.out
  } 
  if(dim(res)[1] == 0){
  mylist <- list()
  mylist[[1]] <- out.s1
  mylist[[2]] <- thrs
  mylist[[3]] <- "No QTL detected"
}
  return(mylist)
})



output$IMsummary <- renderPrint({
    validate(
      need(input$file1 != "", "Upload a cross file to begin")
    )
   IMapping()[[3]]
  })

waitlist <- reactive({
  vals <- c("Counting backwards from infinity",
                "So, do you come here often?",
                "Just stalling to simulate activity...",
                "Waiting for approval from Bill Gates...",
                "Transporting you into the future one second at a time...",
                "About your data... I lost it... in a volcano.",
                "Loading the enchanted bunny...",
                "Please be patient. The program should finish loading in six to eight weeks.")
  indexnr <- runif(8,1,9) %>% trunc()
  vals[indexnr[1]]
})

###############################
###### Composite Interval mapping
############################
CIMapping <- reactive({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(input$phenosel2 != "", waitlist()
         )
  )
  if(input$ngen > 5){
    # cross <- geno() %>% convert2riself()
    cross <- geno()
    cross <- cross %>% jittermap()
    cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
  } else {
    cross <- geno()
    cross <- cross %>% jittermap()
    cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
  }
  out.s1perm <- scanone(cross, pheno.col = input$phenosel2, n.perm = 2500, n.cluster = 4, method = "hk")
  out.s1  <- scanone(cross, pheno.col = input$phenosel2, method = "hk")###############################
  thrs <- summary(out.s1perm , alpha=c(0.05, 0.01))
  # QTL detect at alpha = 0.05
  res <- summary(out.s1, perms=out.s1perm, alpha=0.05)
  ## Cim using n.covar = result IM + 2
  if(dim(res)[1] == 0){
    n.mcovar <- 1 
  }
  if(dim(res)[1] != 0){
    n.mcovar <- dim(res)[1] 
  }
  out.s2 <- cim(cross, pheno.col = input$phenosel2, method = "hk", n.marcovar = n.mcovar)
  # out.s2 <-  stepwiseqtl(cross, pheno.col = input$phenosel2, method = "hk", penalties = summary(out.s1perm, alpha = 0.05)[1], 
  #                        max.qtl = 7, additive.only = T, keeplodprofile = TRUE)
  
  res2 <- summary(out.s2, perms=out.s1perm, alpha=0.05)
  # res2 <- summary(out.s2)
  # if(attributes(out.s2)$pLOD != 0){
  if(dim(res2)[1] > 0){
    mqtl <- makeqtl(cross,res2[,1], res2[,2], what="prob")
    QTLnames <- unlist(mqtl$altname)
    form <- paste(QTLnames, collapse="+")
    form <- paste("y~",form, sep="")
    fit <- fitqtl(cross, pheno.col = input$phenosel2, qtl = mqtl, method = "hk", model = "normal",
                  formula = form, get.ests = T)
    # mqtl <- makeqtl(cross,res2$chr, res2$pos, what="prob")
    # fit <- fitqtl(cross, pheno.col = input$phenosel2, qtl = mqtl, method = "hk", model = "normal",
    #               formula = attributes(res2)$formula, get.ests = T)
    estimate.out <- summary(fit)
    # lodprofile <- do.call(rbind,attributes(out.s2)$lodprofile %>% unname) 
    mylist <- list()
    mylist[[1]] <- out.s2
    # mylist[[1]] <- lodprofile
    mylist[[2]] <- thrs
    mylist[[3]] <- estimate.out
    mylist[[4]] <- out.s1
  } 
  # if(attributes(out.s2)$pLOD == 0){
  if(dim(res2)[1] == 0){
    mylist <- list()
    # mylist[[1]] <- out.s2
    mylist[[1]] <- out.s2
    mylist[[2]] <- thrs
    mylist[[3]] <- "No QTL detected"
    mylist[[4]] <- out.s1
  }
  return(mylist)
})

output$CIMsummary <- renderPrint({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(CIMapping()[[3]] != "No QTL detected","No QTL detected")
  )
  CIMapping()[[3]]
})


#######################
####QTL visualisation##
#######################

output$distPlot <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  plotdata <- IMapping()[[1]]
  plotdata$data_id <- tolower(row.names(plotdata))
  plotdata$tooltip <- row.names(plotdata)
  thr <- IMapping()[[2]]
  plotdata$marker <- row.names(IMapping()[[1]])
  plotdata_set <- filter(plotdata, chr %in% input$chromSelect)
     p <- ggplot(data = plotdata_set, aes(x=pos, y=lod, tooltip = tooltip, data_id = data_id))+
      geom_line(aes(group = 1))+geom_rug(data = plotdata_set[which(!grepl("loc",plotdata_set$marker)),], sides = "b")+
      facet_wrap(~chr, nrow=3, scales = "free_x")+
      geom_point_interactive(color="orange", size=0.1)+
      theme_dark() +
      geom_hline(yintercept = thr[1], lwd = 0.5, lty = 2, col = "white") +
      geom_hline(yintercept = thr[2], lwd = 0.5, lty = 2, col = "orange") +
       xlab("Position (cM)") + ylab("LOD")
    return(ggiraph(code = {print(p)}, zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;"))
  })
  
output$distPlot2 <- renderggiraph({
  validate(
    need(input$file1 != "", "Upload a cross file to begin"),
    need(CIMapping()[[1]] != "No qtl detected", "No QTL detected")
  )
  data1  <- CIMapping()[[1]]
  data2  <- CIMapping()[[4]]
  data1$method <- "CIM"
  data2$method <- "IM"
  ## Get data in right format for plot
  # data1$marker <- row.names(data1)
  # data2$marker <- row.names(data2)
  # extra_df <- filter(data2, !marker %in% data1$marker)
  # extra_df$lod <- 0
  # extra_df$method <- "CIM"
  # data3 <-  rbind(extra_df, data2) %>% arrange(chr,pos)
  # plotdata <- rbind(data1,data3)
  plotdata <- rbind(data1,data2)
  plotdata$data_id <- tolower(row.names(plotdata))
  plotdata$tooltip <- row.names(plotdata)
  plotdata$marker <- row.names(plotdata)
  thr <- CIMapping()[[2]]
  plotdata_set <- filter(plotdata, chr %in% input$chromSelect2)
  p <- ggplot(data = plotdata_set, aes(x=pos, y=lod, colour = method, tooltip = tooltip, data_id = data_id, group = method))+
    geom_line()+geom_rug(data = plotdata_set[which(!grepl("loc",plotdata_set$marker)),], sides = "b")+
    facet_wrap(~chr, nrow=3, scales = "free_x")+
    geom_point_interactive(color="black", size=0.1, alpha = 1)+
    theme_dark() +
    geom_hline(yintercept = thr[1], lwd = 0.5, lty = 2, col = "white") +
    geom_hline(yintercept = thr[2], lwd = 0.5, lty = 2, col = "orange") + 
    scale_color_manual(values = c("darkred","darkblue")) +
    xlab("Position (cM)") + ylab("LOD")
  return(ggiraph(code = {print(p)}, width =12, height = 6,zoom_max = 2, tooltip_offx = 20, tooltip_offy = -10, hover_css = "fill:black;stroke-width:1px;stroke:wheat;cursor:pointer;alpha:1;"))
})

output$lodthr1 <- renderPrint({
  cat("LOD threshold at alpha 0.05 = ",IMapping()[[2]][1])
})

output$lodthr2 <- renderPrint({
  cat("LOD threshold at alpha 0.05 = ",CIMapping()[[2]][1])
})

output$save <- renderUI({
  validate(
    need(input$file1 != "", "Upload a cross file to begin")
  )
  list(
    selectInput("datatype", label = "Select data to export", 
                choices = c("IM LOD","CIM LOD"), 
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
  if (input$datatype == 'IM LOD'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('IM_LOD_',input$phenosel,"_", format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
        lodscores <- IMapping()[[1]] %>% as.data.frame()
        lodscores <- data.frame(marker = row.names(lodscores), lodscores)
        sep <- input$sep_type
        write.table(lodscores, file, sep = sep, row.names = FALSE)
        }
      )}
  }
  if (input$datatype == 'CIM LOD'){
    output$downloadData <- {
      downloadHandler(
        filename = function() {paste0('CIM_LOD_',input$phenosel2,"_", format(Sys.time(), "%a %b %d %X"), '.csv') },
        content = function(file) {
          lodscores <- CIMapping()[[1]] %>% as.data.frame()
          lodscores <- data.frame(marker = row.names(lodscores), lodscores)
          sep <- input$sep_type
          write.table(lodscores, file, sep = sep, row.names = FALSE)
        }
      )}
  }
})



})