library(shiny)
library(data.table)
library(qpcR)
library(ggplot2)
library(googledrive)
library(reshape2)
library(dplyr)

########################### SHINY APP FOR QPCR ANALYSIS ############################################
####################################################################################################

ui <- fluidPage(
  titlePanel(strong("qPCR Analysis")),
  sidebarLayout(
    sidebarPanel(
      checkboxInput("taqman", "Taqman"),
      ###### Sidebar panel appears only of Taqman is checked ########
      conditionalPanel(
        condition = "input.taqman == true",
        fileInput("taqmancsv",
                  label = "Upload CSV here",
                  multiple = FALSE),
        numericInput("ct", "Enter PCR number of cycles", value = 40),
        textInput("endoC", "Enter endogenous control", value = "RNAseP"),
        fileInput("dataendoC", "Enter qPCR results for endoC"),
        fileInput("taqwellid", label = "Upload ID_well.csv here"),
        fileInput("taqidres", label = "Upload ID_results.csv here")
      ),
      ###### Sidebar panel appears only of SYBR is checked ########
      checkboxInput("sybr", "SYBR"),
      conditionalPanel(
        condition = "input.sybr == true",
        fileInput("sybrcsv",
                  label="Upload CSV here",
                  multiple = FALSE),
        radioButtons("isderiv", "Isderiv?",
                     choices = c(True = TRUE,
                                 False = FALSE),
                     selected = TRUE),
        fileInput("wellid",
                  label="Upload ID_well.csv here"),
      )
    ),
    mainPanel(
      ###### Main panel appears only of Taqman is checked ########
      conditionalPanel(
        condition = "input.taqman == true",
        tabsetPanel(
          id = "taq",
          type = "tabs",
          tabPanel("Data", tableOutput("incycles")),
          tabPanel("ID_Well", tableOutput("taqmanwell")),
          tabPanel("ID_Result", tableOutput("taqmanidres")),
          tabPanel("Info", tableOutput("info"))
        )
      ),
      ###### Main panel appears only of SYBR is checked ########
      conditionalPanel(
        condition = "input.sybr == true",
        tabsetPanel(
        id = "syb",
        type = "tabs",
        tabPanel("Gene", tableOutput("tables")),
        tabPanel("Plate Map", tableOutput("well")),
        tabPanel("TM Plots", plotOutput("TM"),
                 downloadButton("downloadplot", "Download PNG")),
        tabPanel("TM Table", tableOutput("tmtable"),
                 downloadButton("downloadtable", "Download"),
                 actionButton("send2drive", "Send to Google Drive", icon = icon("google-drive", lib="font-awesome"))
        )
      )),
    )
  )
)

server <- function(input, output) {
  
  ################################## READ AMPLIFICATION INPUT #################################
  ## 1) Read file and display table
  # a) Read file
  cyclesInput <- reactive({
    dat <- input$taqmancsv
    if (!is.null(dat)){
      data <- read.csv(dat$datapath, sep = ';', header = TRUE, dec = ',');
      data$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$taqmancsv$name))
      colnames(data)[1] <- "Gene"
      return(data)
    }else{
      return (NULL);
    }  
  })
  
  #b) Display table
  output$incycles <- renderTable({
    newData = cyclesInput()
    if(is.null(newData)){
      return();
    }else{
      newData;
    }
    })
  
  ################################## READ TM FILE #############################################
  ## 2) Display table 
  ## a) Define reading function
  rawInputData <- reactive({
    rawData <- input$sybrcsv
    if (!is.null(rawData)){
      'tbl_list <- lapply(input$files$datapath, read.csv, header=TRUE, sep=";")
      for (tbl in tbl_list){
        tbl$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$csvs$name))
      }
      df <- do.call(rbind, tbl_list)
      return(df)'
      data <- read.csv(rawData$datapath, sep = ';', header = TRUE, dec = ',');
      data$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$sybrcsv$name))
      colnames(data)[1] <- "Gene"
      return(data)
    }else{
      return (NULL);
    }  
  })
  ## b) Call function and display in app
  output$tables <- renderTable({
    newData = rawInputData()
    if(is.null(newData)){
      return();
    }else{
      newData;
    }
  })
  
  ############################## READ TAQMAN ID-WELL FILE ############################################
  ## 3) Well id
  taqwellID <- reactive({
    raw <- input$taqwellid
    if (!is.null(raw)){
      data = read.csv(raw$datapath, sep = ',', header = TRUE, dec=',');
      colnames(data)[3] <- "ID"
      return(data)
    }else{
      return (NULL);
    }  
  })
  ## b) Call function and display in app
  output$taqmanwell <- renderTable({
    newWell = taqwellID()
    if(is.null(newWell)){
      return();
    }else{
      newWell;
    }
  })
  
  
  ################################## READ SYBR ID-WELL FILE ###########################################
  ## 3) Well id
  wellID <- reactive({
    raw <- input$wellid
    if (!is.null(raw)){
      data = read.csv(raw$datapath, sep = ',', header = TRUE, dec=',');
      ## Correct Well names, adds new column called Well2 to data 
      data$Well2<-sub('(^[A-Z]+)0', "\\1", data$Well, perl=TRUE) #replace f.i. A01 by A1
      return(data)
    }else{
      return (NULL);
    }  
  })
  ## b) Call function and display in app
  output$well <- renderTable({
    newWell = wellID()
    if(is.null(newWell)){
      return();
    }else{
      newWell;
    }
  })
  
  ############################# READ TAQMAN ID_RESULTS FILE ###################################
  ## Read file
  taqIDRes <- reactive({
    idres <- input$taqidres
    if (!is.null(idres)){
      data <- read.csv(idres$datapath, sep = ',', header = TRUE, dec=',');
      data <- data[,-3] 
      return(data)
    }
  })
  
  ## Display table 
  output$taqmanidres <- renderTable({
    res <- taqIDRes()
    if(is.null(res)){
      return();
    }else{
      res;
    }
  })
  
  ############################# COMBINE ID_WELL AND ID_RESULTS FILE #################################
  constVarendoC <- reactive({
    wellid <- taqwellID()
    idres <- taqIDRes()
    info <-merge(idres,wellid,by="ID")
    info$Well2<-sub('(?<![0-9])0*(?=[0-9])', '', info$Well, perl=TRUE)
    info<-as.data.frame(info)
    
    # set constant variables from endogenous control data
    C<-info[which(info$Target==input$endoC),] # select target
    C<-C[order(C$Well2),]
    return(C)
  })
  
  output$info <- renderTable({
    inf <- constVarendoC()
    if(is.null(inf)){
      return();
    }else{
      inf;
    }
  })
  
  ############################## READ ENDOC QPCR DATA ##########################################
  ## Read file
  taqDim <- reactive({
    d.endoC <- input$dataendoC
    if (!is.null(d.endoC)){
      raw <- read.csv(d.endoC$datapath, sep = ',', header = TRUE, dec=',');
      raw<-raw[,-1]
      raw<-raw[,C$Well2]
      l<-round_any(rowMaxs(as.matrix(raw),value=T)[input$ct], 100, f = ceiling)
      l<-as.integer(l)
      r<-dim(raw)
      r<-r[2]
      r<-as.integer(r)
      lr <- list(l,r)
      print(lr)
      return(lr)
      }
  })
  
  
  ############################## PICK GENE COLUMNS FROM FLUORESCENCE FILE ###############################3
  
  ## 4) Select target gene and match columns in fluorescence file
  # a) Define function
  matchTarget <- reactive({
    Well <- wellID()
    melt_gene <- rawInputData()
    target <- unique(melt_gene$Gene)
    #Select matching target rows from well_id
    well_gene <- Well[which(Well$Target==target),] # select target
    #Select fluorescence matching columns
    fluos_gene<-melt_gene[,well_gene$Well2]
    return(fluos_gene)
    })
  
  ## Prepare data
  fluoTemp <- reactive({
    fluos_gene <- matchTarget()
    melt_gene <- rawInputData()
    fluo_temp <- cbind(melt_gene$Temperature, fluos_gene)
    colnames(fluo_temp)[1] <- "Temperature"
    return(fluo_temp)
  })

################################# MELTING CURVE PLOTS ##############################################
  ## 5) Plot output
  ## a) Define function
  plotInput <- reactive({
    data.gene <- fluoTemp()
    data.gene.m <- melt(data.gene, "Temperature")
    ggplot(data.gene.m, aes(Temperature, value)) + 
      geom_line(color="darkgreen", size=1) +
      facet_wrap(~variable, scales = "free")
  })
  ## b) Call in app
  output$TM <- renderPlot({
    print(plotInput())
      })
  
  ## 6) Download plot
  output$downloadplot <- downloadHandler(
    filename = paste0(gsub(".csv", "", basename(input$sybrcsv$name)),"_plots.png"),
    content = function(file){
      ggsave(file, plotInput(), device = "png")
    },
    contentType = "png"
    )
  
  ############################## TM RESULTS ##################################
  rep.col<-function(x,n) {
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  
  ## 7) Load data and generate output table
  ## a) Generate table
  TMTable <- reactive({
    melt_gene <- rawInputData()
    fluos_gene <- matchTarget()
    a<-dim(fluos_gene)
    b<-a[2]
    c<-as.integer(b+1)
    d<-as.integer(2*b)
    data_gene <- cbind(fluos_gene,rep.col(melt_gene$Temperature,b))
    res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
    ## Wells showing unique Tm with Tm/area info
    dfList <- list()
    for (w in 1:length(res_gene)){
      df<-as.data.frame(res_gene[[w]])
      df<-df[!is.na(df$Tm),6:7]
      dfList[[w]] <- df
    }
    
    names(dfList)<-names(res_gene)
    
    multiTm <- which(sapply(dfList, nrow) != 1)
    length(multiTm) #don't run next step if length multiTm = 0
    if (length(multiTm >0)) {
      dfList<-dfList[-multiTm]
    }
    
    results_gene<-do.call(rbind, dfList)
    results_gene$Target <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$sybrcsv$name))
    return(results_gene)
  })
  
  # b) Display table in app
  output$tmtable <-renderTable({
    t <- TMTable()
    return(t)
  })
  
  ## 8) Download table
  output$downloadtable <- downloadHandler(
    filename = paste0(gsub(".csv", "", basename(input$sybrcsv$name)), "_results.csv"),
    content = function(fname){
      write.csv(TMTable(), fname)}
    )
  'fullpath<-paste("/drive/my-drive", "results_SYBR", sep="/")
  drive_mkdir(fullpath, overwrite = FALSE)#the new folder will be called "results"
  
  output$send2drive <-observeEvent(
    t <- TMTable(),
    drive_upload(t, path = fullpath)
  )'
}
  
shinyApp(ui = ui, server = server)