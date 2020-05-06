library(shiny)
library(data.table)
library(qpcR)
library(ggplot2)

ui <- fluidPage(
  titlePanel(strong("qPCR Analysis")),
  sidebarLayout(
    sidebarPanel(
      fileInput("csvs",
                label="Upload CSVs here",
                multiple = TRUE),
      radioButtons("isderiv", "Isderiv?",
                   choices = c(True = TRUE,
                               False = FALSE),
                   selected = TRUE),
      fileInput("wellid",
                label="Upload well_id.csv here")
    ),
    mainPanel(
      tabsetPanel(
        id = "abcd",
        type = "tabs",
        tabPanel("Gene", tableOutput("tables")),
        tabPanel("Well ID", tableOutput("well")),
        tabPanel("Match", tableOutput("match")),
        tabPanel("Plots", plotOutput("TM"),
                 downloadButton("downloadplot", "Download Plots")),
        tabPanel("Summary", tableOutput("summary")),
        tabPanel("TM Table", tableOutput("tmtable"),
                 downloadButton("downloadtable", "Download Table"))
      )
    )
  )
)

server <- function(input, output) {
  ################################## TM FILE #############################################
  ## 2) Display table 
  ## a) Define reading function
  rawInputData <- reactive({
    rawData <- input$csvs
    if (!is.null(rawData)){
      'tbl_list <- lapply(input$files$datapath, read.csv, header=TRUE, sep=";")
      for (tbl in tbl_list){
        tbl$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$csvs$name))
      }
      df <- do.call(rbind, tbl_list)
      return(df)'
      data <- read.csv(rawData$datapath, sep = ';', header = TRUE, dec = ',');
      data$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$csvs$name))
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
  
  ################################## WELL-ID FILE ###########################################
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
  
  ############################## PICK GENE COLUMNS FROM FLUORESCENCE FILE ###############################3
  
  ## 4) Select target gene and match columns in fluorescence file
  # a) Define function
  matchTarget <- reactive({
    Well <- wellID()
    melt_gene <- rawInputData()
    target <- unique(melt_gene$X)
    #Select matching target rows from well_id
    well_gene <- Well[which(Well$Target==target),] # select target
    #Select fluorescence matching columns
    fluos_gene<-melt_gene[,well_gene$Well2]
    return(fluos_gene)
    })
  
  #b) Display in app
  output$match <- renderTable({
    matchFluo = matchTarget()
    if(is.null(matchFluo)){
      return();
    }else{
      matchFluo;
    }
  })
  
  ## Prepare data
  fluoTemp <- reactive({
    fluos_gene <- matchTarget()
    melt_gene <- rawInputData()
    fluo_temp <- cbind(melt_gene$Temperature, fluos_gene)
    colnames(fluo_temp)[1] <- "Temperature"
    return(fluo_temp)
  })

  output$summary <-renderTable({
    d <- fluoTemp()
    return(d)
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
    filename = paste0(gsub(".csv", "", basename(input$csvs$name)),"_plots.png"),
    content = function(file){
      ggsave(file, plotInput(), device = "png")
    },
    contentType = "png"
    )
  
  ############################## TM RESULTS ##################################
  rep.col<-function(x,n) {
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  
  TMTable <- reactive({
    melt_gene <- rawInputData()
    fluos_gene <- matchTarget()
    a<-dim(fluos_gene)
    b<-a[2]
    c<-as.integer(b+1)
    d<-as.integer(2*b)
    data_gene <- cbind(fluos_gene,rep.col(melt_gene$Temperature,b))
    #print(data_gene)
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
    results_gene$Target <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$csvs$name))
    return(results_gene)
  })
  
  output$tmtable <-renderTable({
    t <- TMTable()
    return(t)
  })
  
  ## 6) Download table
  output$downloadtable <- downloadHandler(
    filename = paste0(gsub(".csv", "", basename(input$csvs$name)), "_results.csv"),
    content = function(fname){
      write.csv(TMTable(), fname)}
    )
}
  
shinyApp(ui = ui, server = server)