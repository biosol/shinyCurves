server <- function(input, output) {
  ################################## READ AMPLIFICATION INPUT #################################
  ## 1) Read file and display table
  # a) Read file
  cyclesInput <- reactive({
    dat <- input$taqmancsv
    if (!is.null(dat)){
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=";", dec=",")
      tbl_list <- lapply(tbl_list, function(x){x[1]<-NULL;x})
      return(tbl_list)
    }
  })
  
  ## Gene list
  geneList <- reactive({
    dat <- input$taqmancsv
    lst <- lapply(dat$name, FUN = function(x) gsub(pattern = ".*[_]([^.]+)[.].*", replacement = "\\1",
                                                   basename(dat$name)))
    lst <- unlist(lst)
    out <- unique(lst)
    return(out)
  })
  
  'cyclesInput <- reactive({
    dat <- input$taqmancsv
    if (!is.null(dat)){
      data <- read.csv(dat$datapath, sep = ';', header = TRUE, dec = ",");
      data$X <- gsub(".*[_]([^.]+)[.].*", "\\1", basename(dat$name))
      colnames(data)[1] <- "Gene"
      return(data)
    }else{
      return (NULL);
    }  
  })'
  
  '#b) Display table
  output$incycles <- renderTable({
    newData = cyclesInput()
    if(is.null(newData)){
      return();
    }else{
      newData;
    }
  })'
  
  ################################## READ TM FILE #############################################
  ## 2) Display table
  
  # a) Read TM files (multiple upload will be saved in list)
  TMinput <- reactive({
    dat <- input$sybrcsv
    if (!is.null(dat)){
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=";", dec=",")
      tbl_list <- lapply(tbl_list, function(x){x[1]<-NULL;x})
      return(tbl_list)
    }
  })
  
  ## Gene list
  SybrGeneList <- reactive({
    dat <- input$sybrcsv
    lst <- lapply(dat$name, FUN = function(x) gsub(pattern = ".*[_]([^.]+)[.].*", replacement = "\\1",
                                                   basename(dat$name)))
    lst <- unlist(lst)
    out <- unique(lst)
    return(out)
  })
  
  ## a) Define reading function
  'rawInputData <- reactive({
    rawData <- input$sybrcsv
    if (!is.null(rawData)){
      data <- read.csv(rawData$datapath, sep = ';', header = TRUE, dec = ",");
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
  })'
  
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
      return(lr)
    }
  })
  
  ############################ GENERAL PLOTS FOR AMPLIFICATION ###############################################
  generalPlots <- function(){
    info <- combIDwellIDres()
    lr <- taqDim()
    genes <- geneList()
    
    pltList <- list()
    for (i in 1:length(genes)){
      df<-info[which(info$Target==genes[i]),] # select target
      df<-df[order(df$Well2),]
    
      ## load raw data (fluorescence-RFU per temperature) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ## Read fluorescence data for target gene
      inp <- cyclesInput()
      raw <- inp[i][[1]]
      raw<-raw[,-1]
      raw<-raw[,df$Well2]
      raw_s <- stack(raw)
      raw_s$cycles<-rep(1:input$ct,lr[[2]]) # rep(1:cycle number,sample number x2 replicates)
      colnames(raw_s)<-c("RFU","Well2","Cycles")
      merge<-merge(raw_s,df,by="Well2")
      
      ###plot
      pltName <- paste(genes[i])
      pltList[[pltName]] = ggplot(data = merge, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
        geom_line()+
        ylim(-100, lr[[1]])+ # manually adjust after seeing plot scale - keep same for all plots
        ggtitle(paste(genes[i]))
    }
    return(pltList)
  }
  
  # Render
  output$genplots <- renderPlot({
    plots <- generalPlots()
    plot_grid(plotlist = plots, ncol = 1)
  })
  
  output$downgen <- downloadHandler(
    filename = "GeneralPlots.pdf",
    content = function(file){
      plots <- generalPlots()
      p <- plot_grid(plotlist = plots, ncol = 1)
      save_plot(file, p, ncol = 1, base_height = 8, base_width = 6)
    }
  )
  
  ############################# AMPLIFICATION CURVE PLOTS FOR EACH GENE (INDET) ################################# 
  indetPlots <- reactive({
    info <- combIDwellIDres()
    lr <- taqDim()
    genes <- geneList()
    
    defPltLs <- list()
    for (i in 1:length(genes)){
      ## Select target gene wells
      df<-info[which(info$Target==genes[i]),] # select target
      df<-df[order(df$Well2),]
      ## Read fluorescence data for target gene
      inp <- cyclesInput()
      raw <- inp[i][[1]]
      raw <- raw[,-1]
      raw <- raw[,df$Well2]
      raw_s <- stack(raw)
      raw_s$cycles<-rep(1:input$ct,lr[[2]]) # rep(1:cycle number,sample number x2 replicates)
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
        pltName <- paste(genes[i],"_",names(merge_r_persample)[[z]])
        pltList[[pltName]] <- ggplot(data = merge_pn, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
          geom_line()+
          geom_line(data = merge_r_persample[[z]], aes(x=Cycles, y=RFU, by=Well2, color=Interpretation))+
          ylim(-100, lr[[1]])+ 
          theme(legend.position = "none")+ 
          ggtitle(paste(genes[i],"_",names(merge_r_persample)[[z]]))
      }
      defPltLs[[i]] <- pltList
    } 
    return(defPltLs)
  })
    
  ## Render plots in app
  output$n1 <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[1]], ncol = 2)
  })
  
  ## Download plots
  output$downln1 <- downloadHandler(
    filename = ("N1_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[1]], ncol = 3)
      save_plot(file, p, ncol = 3, base_height = 8, base_width = 4)
    }
  )
  
  # N2 tab
  output$n2 <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[2]], ncol = 2)
  })
  
  ## Download plots
  output$downln2 <- downloadHandler(
    filename = ("N2_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[2]], ncol = 3)
      save_plot(file, p, ncol = 3, base_height = 8, base_width = 4)
    }
  )
  
  #RNAsep tab
  output$rnasep <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[3]], ncol = 2)
  })
  
  ## Download plots
  output$downrnasep <- downloadHandler(
    filename = ("RNAseP_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[3]], ncol = 3)
      save_plot(file, p, ncol = 3, base_height = 8, base_width = 4)
    }
  )
  
  
  ############################## PICK GENE COLUMNS FROM FLUORESCENCE FILE ###############################3
  ## 4) Select target gene and match columns in fluorescence file
  # a) Define function
  matchTarget <- reactive({
    Well <- wellID()
    genes <- SybrGeneList()
    melt_gene <- TMinput()
    #target <- unique(melt_gene$Gene)
    matchLs <- list()
    for (i in 1:length(genes)){
      #Select matching target rows from well_id
      well_gene <- Well[which(Well$Target==genes[i]),] # select target
      #Select fluorescence matching columns
      fluos_gene<-melt_gene[[i]][,well_gene$Well2]
      matchLs[[i]] <- fluos_gene
    }
    return(matchLs)
  })
  
  output$fluosgene <- renderTable({
    matchTarget()
  })
  
  ## Prepare data
  fluoTemp <- reactive({
    fluos_gene <- matchTarget()
    melt_gene <- TMinput()
    genes <- SybrGeneList()
    LsforPlot <- list()
    for (i in 1:length(genes)){
      fluo_temp <- cbind(melt_gene[[i]]$Temperature, fluos_gene[[i]])
      colnames(fluo_temp)[1] <- "Temperature"
      LsforPlot[[i]] <- fluo_temp
    }
    return(LsforPlot)
  })
  
  output$fluotemp <-renderTable({
    fluoTemp()
  })
  
  ################################# MELTING CURVE PLOTS ##############################################
  # Function to repeat temperature columns
  ## 5) Plot  TM curves
  # a) Define function
  rep.col<-function(x,n) {
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  
  TMPlots <- reactive({
    melt_gene <- TMinput()
    fluos_gene <- matchTarget()
    genes <- SybrGeneList()
    
    pltList <- list()
    for (i in 1:length(genes)){
      a<-dim(fluos_gene[[i]])
      b<-a[2]
      c<-as.integer(b+1)
      d<-as.integer(2*b)
      data_gene <- cbind(fluos_gene[[i]],rep.col(melt_gene[[i]]$Temperature,b))
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
      pltList[[i]] <- res_gene
    }
    return(pltList)
  })
  
  'TMPlots <- function(){
    #Combine all input genes
    melt_gene <- TMinput()
    melt_gene_df <- list.cbind(melt_gene)
    melt_gene_df<- melt_gene_df[, !duplicated(colnames(melt_gene_df))]
    
    # Combine gene fluo data
    fluos_gene <- matchTarget()
    fluos_gene_df <- list.cbind(fluos_gene)
    fluos_gene_df <- fluos_gene_df[, !duplicated(colnames(fluos_gene_df))]
    
    # Plot all Tm together (for all genes)
    a<-dim(fluos_gene_df)
    b<-a[2]
    c<-as.integer(b+1)
    d<-as.integer(2*b)
    data_gene <- cbind(fluos_gene_df,rep.col(melt_gene_df$Temperature,b))
    res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
  }'
  
  output$sybrN <- renderPlot({
    TMPlots()[[1]]
  })
  
  output$sybrRdrp <- renderPlot({
    TMPlots()[[2]]
  })
  
  output$sybrRpp30 <- renderPlot({
    TMPlots()[[3]]
  })
  
  output$sybrS <- renderPlot({
    TMPlots()[[4]]
  })
  
  output$downsybrN <- downloadHandler(
    filename = ("N_TMPlots.pdf"),
    content = function(file){
      p <- TMPlots()[[1]]
      save_plot(file, p, base_height = 8, base_width = 4)
    })
  
  # b) Display output in app
  'output$TM <- renderPlot({
    TMPlots()
  })'
  
  ## c) Download plot
  output$downloadTMplots <- downloadHandler(
    filename = "TMplots.pdf",
    content = function(file){
      pdf(file)
      TMPlots()
      dev.off()
    }
  )
  
  ############################## TM RESULTS ##################################
  ## 7) Load data and generate output table
  ## a) Generate table for all genes
  TMTable <- reactive({
    #Combine all input genes
    melt_gene <- TMinput()
    fluos_gene <- matchTarget()
    genes <- SybrGeneList()
    
    results_all <- list()
    for (i in 1:length(genes)){
      a<-dim(fluos_gene[[i]])
      b<-a[2]
      c<-as.integer(b+1)
      d<-as.integer(2*b)
      data_gene <- cbind(fluos_gene[[i]],rep.col(melt_gene[[i]]$Temperature,b))
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
      
      names(res_gene)<-colnames(fluos_gene[[i]])
      ## Wells showing unique Tm with Tm/area info
      dfList <- list()
      for (j in 1:length(res_gene)){
        df<-as.data.frame(res_gene[[j]])
        df<-df[!is.na(df$Tm),6:7]
        #df$Target <- genes[[i]]
        dfList[[j]] <- df
      }
      
      names(dfList)<-names(res_gene)
      
      multiTm <- which(sapply(dfList, nrow) != 1)
      length(multiTm) #dont run next step if length multiTm = 0
      if (length(multiTm > 0)) {
        dfList<-dfList[-multiTm]
      }
      results_gene<-do.call(rbind, dfList)
      results_gene$Target <- rep(genes[[i]],dim(results_gene)[1])
      results_all[[i]] <- results_gene
    }
    results <- list.stack(results_all)
    return(results)
  })
  
  # b) Display table in app
  output$tmtable <-renderTable({
    TMTable()
    })
  
  ## c) Download table
  output$downloadtable <- downloadHandler(
    filename = "TMtable.csv",
    content = function(fname){
      write.csv(TMTable(), fname, quote = F, row.names = F)}
  )
  
  'fullpath<-paste("/drive/my-drive", "results_SYBR", sep="/")
  drive_mkdir(fullpath, overwrite = FALSE)#the new folder will be called "results"
  
  output$send2drive <-observeEvent(
    t <- TMTable(),
    drive_upload(t, path = fullpath)
  )'
}