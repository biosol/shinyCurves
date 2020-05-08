
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
  ## a) Combine ID well and ID results
  combIDwellIDres <- reactive({
    wellid <- taqwellID()
    idres <- taqIDRes()
    info <-merge(idres,wellid,by="ID")
    info$Well2<-sub('(?<![0-9])0*(?=[0-9])', '', info$Well, perl=TRUE)
    info<-as.data.frame(info)
    return(info)
  })
    
    # b) Set constant variables from endogenous control data
  constVarendoC <- reactive({
    info <- combIDwellIDres()
    C<-info[which(info$Target==input$endoC),] # select target
    C<-C[order(C$Well2),]
    return(C)
  })
  
  output$c <- renderTable({
    C <- constVarendoC()
    C
  })
  
  # c) Display info table
  output$info <- renderTable({
    inf <- combIDwellIDres()
    if(is.null(inf)){
      return();
    }else{
      inf;
    }
  })
  
  ############################## READ ENDOC QPCR DATA ##########################################
  ## Read file
  taqDim <- reactive({
    C <- constVarendoC()
    if (!is.null(input$dataendoC)){
      raw <- read.csv(input$dataendoC$datapath, sep = ";", header = TRUE, dec=",")
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
  
  ############################# AMPLIFICATION CURVE PLOTS FOR EACH GENE (INDET) ################################# 
  indetPlots <- reactive({
    ## Define nrows for plot
    result <- taqIDRes()
    indet<-as.integer(dim(result[result$Interpretation=="Repeat",])[1])
    if (indet<=15) {
      indet=as.integer(5)
    } else {
      indet=as.integer(7)
    }
    
    info <- combIDwellIDres()
    lr <- taqDim()
    ## Select target gene wells
    df<-info[which(info$Target==gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$taqmancsv$name))),] # select target
    df<-df[order(df$Well2),]
    
    ## Read fluorescence data for target gene
    raw <- cyclesInput()
    raw <- raw[,-1]
    raw <- raw[,df$Well2]
    raw_s <- stack(raw)
    raw_s$cycles<-rep(1:input$ct,lr[2]) # rep(1:cycle number,sample number x2 replicates)
    colnames(raw_s)<-c("RFU","Well2","Cycles")
    merge<-merge(raw_s,df,by="Well2")
    
    ## one-by-one 
    OBO<-split(merge,merge$Interpretation,drop=T)
    merge_pn<-as.data.frame(rbind(OBO[[1]],OBO[[2]]))
    merge_r<-as.data.frame(OBO[[3]])
    merge_r_persample<-split(merge_r,merge_r$ID,drop=T)
    
    ## Plots
    pltList <- list()
    for (z in 1:length(merge_r_persample)){
      pltName <- paste(gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$taqmancsv$name)),"_",names(merge_r_persample)[[z]])
      pltList[[pltName]] <- ggplot(data = merge_pn, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
        geom_line()+
        geom_line(data = merge_r_persample[[z]], aes(x=Cycles, y=RFU, by=Well2, color=Interpretation))+
        ylim(-100, lr[[1]])+ 
        theme(legend.position = "none")+ 
        ggtitle(paste(gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$taqmancsv$name)),"_",names(merge_r_persample)[[z]]))
    }
    plot_grid(plotlist = pltList, ncol = 2, nrow=indet)
  })
  
  ## Render plots in app
  output$indetplots <- renderPlot({
    indetPlots()
  })
  
  ## Download plots
  output$downloadIndet <- downloadHandler(
    filename = paste0(gsub(".*[_]([^.]+)[.].*", "\\1", basename(input$taqmancsv$name)),"_IndetPlots.png"),
    content = function(file){
      ggsave(file, indetPlots(), device = "png")
    },
    contentType = "png"
  )
  
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