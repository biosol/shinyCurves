server <- function(input, output) {
  ################################## READ AMPLIFICATION INPUT #################################
  ## 1) Read file and display table
  # a) Read file
  cyclesInput <- reactive({
    dat <- input$taqmancsv
    if (!is.null(dat)){
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=",", dec=".")
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
  
  ############################## READ TAQMAN ID-WELL FILE ############################################
  ## 3) Well id
  taqwellID <- reactive({
    raw <- input$taqwellid
    if (!is.null(raw)){
      data = read.csv(raw$datapath, sep = ';', header = TRUE, dec=',');
      data[is.na(data)]<-""
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
      data <- read.csv(idres$datapath, sep = ';', header = TRUE, dec=',');
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
      raw <- read.csv(input$dataendoC$datapath, header = TRUE, sep = ",", dec=".")
      raw<-raw[,-1]
      raw<-raw[,C$Well2]
      l<-round_any(rowMaxs(as.matrix(raw),value=T), 100, f = ceiling)
      l <-max(l)
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
    #genes <- c("N1", "N2", "RNAseP")
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
      pltList[[ pltName ]] = ggplot(data = merge, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
        geom_line()+
        ylim(-100, lr[[1]])+ # manually adjust after seeing plot scale - keep same for all plots
        ggtitle(pltName)
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
      merge_pn<-merge[merge$Interpretation!="Repeat",]
      merge_r<-merge[merge$Interpretation=="Repeat",]
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
      p <- plot_grid(plotlist = plots[[1]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
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
      p <- plot_grid(plotlist = plots[[2]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
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
      p <- plot_grid(plotlist = plots[[3]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
    }
  )
  
  
  ############################## PICK GENE COLUMNS FROM FLUORESCENCE FILE ###############################3
  ## 4) Select target gene and match columns in fluorescence file
  # a) Define function
  ## matchTarget por archivo independiente
  matchTarget <- function(tm, gene){
    Well <- wellID()
    well_gene <- Well[which(Well$Target==gene),] # select target
    #Select fluorescence matching columns
    fluos_gene<-tm[,well_gene$Well2]
    return(fluos_gene)
  }
  
  ## Misma función que matchTarget pero con un loop para guardar todos los archivos en la misma tabla
  ## Bastante guarro, pero funciona
  matchAllTarget <- reactive({
    Well <- wellID()
    genes <- SybrGeneList()
    melt_gene <- TMinput()
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
    matchAllTarget()
  })
  
  ## Prepare data
  fluoTemp <- reactive({
    fluos_gene <- matchAllTarget()
    print(fluos_gene)
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
  
  TMPlots <- function(tm, gene){
    fluos_gene <- matchTarget(tm, gene)
    a<-dim(fluos_gene)
    b<-a[2]
    print(b)
    c<-as.integer(b+1)
    d<-as.integer(2*b)
    data_gene <- cbind(fluos_gene,rep.col(tm$Temperature,b))
    if (input$isderiv == TRUE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
    } else if (input$isderiv == FALSE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5))
    }
  }
  
  ##### b) Render plots
  output$sybrN <- renderPlot({
    tm <- TMinput()
    genes <- SybrGeneList()
    TMPlots(tm[[1]], genes[[1]])
  })
  
  output$sybrRdrp <- renderPlot({
    genes <- SybrGeneList()
    tm <- TMinput()
    TMPlots(tm[[2]], genes[[2]])
  })
  
  output$sybrRpp30 <- renderPlot({
    genes <- SybrGeneList()
    tm <- TMinput()
    TMPlots(tm[[3]], genes[[3]])
  })
  
  output$sybrS <- renderPlot({
    genes <- SybrGeneList()
    tm <- TMinput()
    TMPlots(tm[[4]], genes[[4]])
  })
  
  ## c) Downloads
  output$downsybrN <- downloadHandler(
    filename = ("N_TMPlots.pdf"),
    content = function(file){
      pdf(file)
      tm <- TMinput()
      genes <- SybrGeneList()
      p <- TMPlots(tm[[1]], genes[[1]])
      dev.off()
    })
  
  output$downsybrRdrp <- downloadHandler(
    filename = ("Rdrp_TMPlots.pdf"),
    content = function(file){
      pdf(file)
      tm <- TMinput()
      genes <- SybrGeneList()
      p <- TMPlots(tm[[2]], genes[[2]])
      dev.off()
    })
  
  output$downsybrRpp30 <- downloadHandler(
    filename = ("Rdrp_TMPlots.pdf"),
    content = function(file){
      pdf(file)
      tm <- TMinput()
      genes <- SybrGeneList()
      p <- TMPlots(tm[[3]], genes[[3]])
      dev.off()
    })
  
  output$downsybrS <- downloadHandler(
    filename = ("S_TMPlots.pdf"),
    content = function(file){
      pdf(file)
      tm <- TMinput()
      genes <- SybrGeneList()
      p <- TMPlots(tm[[4]], genes[[4]])
      dev.off()
    })
  
  ############################## TM RESULTS ##################################
  ## 7) Load data and generate output table
  ## a) Generate table for all genes
  TMTable <- reactive({
    #Combine all input genes
    melt_gene <- TMinput()
    fluos_gene <- matchAllTarget()
    print(fluos_gene)
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
  
  ## b) Render table
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
  
  ################################# TAQMAN EXCEL ##################################
  #################################################################################
  
  ######### PRINT EMPTY PLATE ###########
  platesTable <- function(cnames, controls){
    col1 <- rep(1:8, each=2)
    col2 <- rep(9:16, each=2)
    col3 <- rep(17:24, each = 2)
    col4 <- rep(25: 32, each=2)
    col5 <- rep(33:40, each=2)
    col6 <- rep(41:48, each=2)
    col7 <- rep(49:56, each=2)
    col8 <- c(rep(57:60, each=2),controls)
    
    df <- cbind(col1,col2,col3,col4,col5,col6,col7, col8)
    rownames(df) <- LETTERS[1:16]
    colnames(df) <- cnames
    return(df)
  }
  
  platePrint <- reactive({
    p1 <- platesTable(1:8, c("NTC","NTC","EXP-5","EXP-5","EXP-4","EXP-4", "EXP-2","EXP-2"))
    p2 <- platesTable(9:16, c("NTC","NTC","EXP-5","EXP-5","EXP-4","EXP-4", "EXP-2","EXP-2"))
    p3 <- platesTable(17:24, c("NTC","NTC"," "," "," "," ", "HeLa","HeLa"))
    
    def <- cbind(p1,p2,p3)
    defdt <- datatable(def, options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:8,
        backgroundColor = "yellow"
        ) %>%
      formatStyle(
        columns = 9:16,
        backgroundColor = "orange"
        ) %>%
      formatStyle(
        columns = 17:24,
        backgroundColor = "lightblue")
    return(f)
  })
  
  output$plate <- DT::renderDataTable(
    platePrint()
  )
  
  ######################### BIORAD RESULTS #################################
  ##### Run Information #####
  bioradruninfo <- reactive({
    res <- input$biorad
    runinfo <- read_xlsx(res$datapath, sheet = "Run Information")
    df <- data.frame(runinfo)
    return(runinfo)
  })

    output$BRruninfo <- renderTable({
    bioradruninfo()
  })
  
  #### Sample Column ####
  bioradsample <- reactive({
    res <- input$biorad
    samples <- read_xlsx(res$datapath, sheet = "Data")
    df <- data.frame(samples)
    return(samples)
  })
  
  ## Well column ###
  inputDF <- function(){
    well <- list()
    counter <- 0
    for (let in LETTERS[1:16]){
      for (nb in 1:24){
        counter = counter +1
        i <- paste(let,nb,sep="")
        well[[counter]] <- i
      }
    }
    
    ## FLuor column ##
    welldf <- data.table(well)
    welldf$fluor <- 'FAM'
    
    ## Target column ###
    n1 <- rep("N1", each=8)
    n2 <- rep("N2", each=8)
    rnasep <- rep("RNAseP", each=8)
    all <- list(n1,n2,rnasep)
    d <- do.call(rbind, Map(data.frame, Target=all))
    welldf$Target <- rep(d$Target, len=384)
    
    ## Content column ###
    welldf$Content <- 'Unknown'
    
    ## Sample column ###
    bio <- bioradsample()
    welldf$Sample <- bio$Sample
    welldf$Cq <- bio$Cq
    
    # Return data table
    return(welldf)
  }
  
  output$inputdf <- renderTable(
    inputDF()
  )
  
  cqPlate <- reactive({
    inp <- inputDF()
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$Cq[s:e]
      s <- s + 24
      e <- e + 24
    }
    dt <- dt %>% 
      mutate_if(is.numeric, round, digits = 3)
    return(dt)
  })
  
    ## Convert to datatable to render nicely in shiny
  printCqPlate <- reactive({
    dt <- cqPlate()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:8,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 9:16,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 17:24,
        backgroundColor = "lightblue")
    return(f)
  })
  
  output$cqplate<- DT::renderDataTable(
    printCqPlate()
  )
  
  #################### Sample Plate ############################
  samplePlate <- reactive({
    inp <- inputDF()
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$Sample[s:e]
      s <- s + 24
      e <- e + 24
    }
    return(dt)
  })
  ## Convert to datatable to render nicely in shiny
  printSamplePlate <- reactive({
    dt <- samplePlate()
    defdt <- datatable(dt, rownames = F,  options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:8,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 9:16,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 17:24,
        backgroundColor = "lightblue")
    return(f)
  })
  
  output$sampleplate<- DT::renderDataTable(
    printSamplePlate()
  )
  
  #################### Check sample plates #######################
  checkSamples <- reactive({
    s <- samplePlate()
    first <- data.frame(d = unlist(s[1:8], use.names = FALSE))
    second <- data.frame(d = unlist(s[9:16], use.names = FALSE))
    third <- data.frame(d = unlist(s[17:24], use.names = FALSE))
    all <- cbind(first, second, third)
    colnames(all) <- c("first", "second", "third")
    def <- head(all, 120)
    
    ## Add "coinciden" to script 
    for (i in 1:nrow(def)){
      if (def[i,"first"] == def[i, "second"] && def[i,"first"] == def[i, "third"]){
        def$Check[i] <- "Coinciden"
      } else {
        def$Check[i] <- "No coinciden"
      }
    }
    return(def)
  })
    
    ## Convert to datatable to render nicely in shiny
  checkSamplesDT <- reactive({
    def <- checkSamples()
    defdt <- datatable(def, rownames = F,  options = list(pageLength = 60))
    f <- defdt %>% 
      formatStyle(
        columns = 1,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 2,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 3,
        backgroundColor = "lightblue"
      ) %>%
      formatStyle(
        columns = 4,
        backgroundColor = styleEqual(c("Coinciden", "No coinciden"), c('seagreen', 'red'))
      )
    return(f)
  })
  
  ##### Render Table ####
  output$samplecheck <- DT::renderDataTable(
    checkSamplesDT()
  )
  
  ################### Analysis Standard Curve #####################
  
  stCurve <- reactive({
    a <- data.frame(matrix(0, nrow = 5, ncol = 19))
    colnames(a)<-c("Sample","Dilution","Copies","log(copies)","N1_dup1","N2_dup1","RNAseP_dup1","N1_dup2","N2_dup2","RNAseP_dup2","N1_avg","N2_avg","RNAseP_avg","N1_logcop","N2_logcop","RNAseP_logcop", "N1_copies", "N2_copies", "RNAseP_copies")
    a$Sample <- c("NTC","C(-)","C(+)10-2","C(+)10-4","C(+)10-5")
    a$Dilution <- c("-","-",100,1000,100000)
    a$Copies <- c("-","-",200000*2/as.numeric(a$Dilution[3]), 200000*2/as.numeric(a$Dilution[4]),200000*2/as.numeric(a$Dilution[5]))
    a$`log(copies)`<- c("-","-", 3.602, 1.602, 0.602)
    
    cq <- cqPlate()
    a$N1_dup1 <- c(cq[9,8],'-',cq[15,8],cq[13,8], cq[11,8])
    a$N2_dup2 <- c(cq[9,16],'-',cq[15,16],cq[13,16], cq[11,16])
    a$RNAseP_dup1 <- c(cq[9,24],cq[15,24],'-','-','-')
    a$N1_dup2 <- c(cq[10,8],'-',cq[16,8],cq[14,8], cq[12,8])
    a$N2_dup2 <- c(cq[10,16],'-',cq[16,16],cq[14,16], cq[12,16])
    a$RNAseP_dup2 <- c(cq[10,24],cq[16,24],'-','-','-')
    a$N1_avg <- c(mean(c(cq[9,8],cq[10,8])), "-", mean(c(cq[15,8], cq[16,8])), mean(c(cq[13,8], cq[14,8])), mean(c(cq[11,8], cq[12,8])))
    a$N2_avg <- c(mean(c(cq[9,16],cq[10,16])), "-", mean(c(cq[15,16], cq[16,16])), mean(c(cq[13,16], cq[14,16])), mean(c(cq[11,16], cq[12,16])))
    a$RNAseP_avg <- c(mean(c(cq[9,24],cq[10,24])), "-", mean(c(cq[15,24], cq[16,24])), mean(c(cq[13,24], cq[14,24])), mean(c(cq[11,24], cq[12,24])))
    return(a)
  })
  
  #### Convert to datatable ####
  stCurveDT <- reactive({
    a <- stCurve()
    ## Transform into data table
    defdt <- datatable(a, rownames = F)
    dt <- defdt %>%
      formatStyle(
        columns = c(5,8,11),
        backgroundColor = "yellow"
      )%>%
      formatStyle(
        columns = c(6,9,12),
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = c(7,10,13),
        backgroundColor = "lightblue"
      )
    
    return(dt)
  })
  
  ##### Render Table ####
  output$stdcurve <- DT::renderDataTable(
    stCurveDT()
  )
  
  ################### Standard Curve Plots ###########################
  stdCoeffs <- function(){
    a <- stCurve()
    n1 <- a$N1_avg[3:5]
    n2 <- a$N2_avg[3:5]
    cp <- a$`log(copies)`[3:5]
    p <- as.data.frame(cbind(as.numeric(n1),as.numeric(n2) , as.numeric(cp)))
    colnames(p) <- c("n1", "n2", "cp")
    
    model_n1 <- lm(cp ~ n1, p)
    n1_coeff1 <- as.numeric(model_n1$coefficients[1])
    n1_coeff2 <- as.numeric(model_n1$coefficients[2])
    n1_coeff <- as.data.frame(cbind(n1_coeff1, n1_coeff2))
    colnames(n1_coeff) <- c("Intercept", "Slope")
    
    model_n2 <- lm(cp ~ n2, p)
    n2_coeff1 <- as.numeric(model_n2$coefficients[1])
    n2_coeff2 <- as.numeric(model_n2$coefficients[2])
    n2_coeff <- as.data.frame(cbind(n2_coeff1, n2_coeff2))
    colnames(n2_coeff) <- c("Intercept", "Slope")
    
    coeff <- as.data.frame(rbind(n1_coeff, n2_coeff))
    rownames(coeff) <- c("N1", "N2")
    
    return(coeff)
  }
  
  stdPlots <- reactive({
    a <- stCurve()
    
    n11 <- a$N1_avg[3:5]
    n22 <- a$N2_avg[3:5]
    cpp <- a$`log(copies)`[3:5]
    p <- as.data.frame(cbind(as.numeric(n11),as.numeric(n22),as.numeric(cpp)))
    colnames(p) <- c("n1", "n2", "cp")
    
    p1 <- ggplot(p, aes(x=n1, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = cp ~ n1,
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(p)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(p)[1]))) +
      ylab("log(Copies)")
    
    p2 <- ggplot(p, aes(x=n2, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = cp ~ n2, 
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(p)[2])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(p)[2]))) +
      ylab("log(Copies)")
    
    s <- grid.arrange(p1, p2, nrow = 1)
    return(s)
  
  })
  
  output$std <- renderPlot(
    stdPlots()
  )
  
  IDWELLtab <- reactive({
    inp <- inputDF()
    wellid <- cbind(inp$well, as.character(inp$Target), inp$Sample)
    colnames(wellid) <- c("Well", "Target", "Sample")
    return(wellid)
  })
  
  output$IDWELL <- renderTable(
    IDWELLtab()
  )
  
  
  ############ Plate Setup MultiChanel #################### 
  setupMultiC <- reactive({
    def <- checkSamples()
    col1 <- rep(paste("S",1:60, sep = ""),each =2)
    col2 <- rep(paste(col1, "-", 1:2, sep=""))
    col3 <- unique(def$first)
    multic <- as.data.frame(cbind(col1, col2, col3))
    colnames(multic)<- c("Sample", "Replicate", "Real_ID")
    return(multic)
  })
  
  setupMultiCDT <- reactive({
    a <- setupMultiC()
    defdt <- datatable(a, rownames = F, 
                       options = list(pageLength = 50))
    f <- defdt %>%
      formatStyle(
        columns = c(1,3),
        backgroundColor = "seagreen"
      )%>%
      formatStyle(
        columns = 2,
        backgroundColor = "lightgreen"
      )
    return(f)
  })
  
  output$setupmultic <- renderDataTable(
    setupMultiCDT()
  )
 
  ############# Hidrolysis Probe Plate ###############
  AnalysisSamples <- reactive({
    
    ## Real Id, Sample, Replicate columns
    a <- setupMultiC()
    
    ## Well columns
    wells_n1 <- vector()
    for (n in 1:8){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n1 <- append(wells_n1, p)
    }
    wells_n1_def <- as.data.frame(head(wells_n1, 120))
    colnames(wells_n1_def) <- "Wells_N1"
    
    wells_n2 <- vector()
    for (n in 9:16){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n2 <- append(wells_n2, p)
    }
    wells_n2_def <- as.data.frame(head(wells_n2, 120))
    colnames(wells_n2_def) <- "Wells_N2"
    
    wells_rnasep <- vector()
    for (n in 9:16){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_rnasep <- append(wells_rnasep, p)
    }
    wells_rnasep_def <- as.data.frame(head(wells_n2, 120))
    colnames(wells_rnasep_def) <- "Wells_RNAseP"
    
    ###### Ct columns 
    cq <- cqPlate()
    cq_n1 <- stack(cq[,1:8])
    cq_n1 <- as.data.frame(head(cq_n1[,1], 120))
    colnames(cq_n1) <- "Cq_N1"
    
    cq_n2 <- stack(cq[,9:16])
    cq_n2 <- as.data.frame(head(cq_n2[,1], 120))
    colnames(cq_n2) <- "Cq_N2"
    
    cq_rnasep <- stack(cq[,17:24])
    cq_rnasep <- as.data.frame(head(cq_rnasep[,1], 120))
    colnames(cq_rnasep) <- "Cq_RNAseP"
    
    #### log(copies) Column
    coeffs <- stdCoeffs()
    
    lgcop_n1 <- vector()
    for (i in 1:nrow(cq_n1)){
      if (is.na(cq_n1[i,])){
        lgcop_n1 <- append(lgcop_n1, NA)
      }else{
        t <- i*coeffs[1,2]+coeffs[1,1]
        print(t)
        lgcop_n1 <- append(lgcop_n1, t)
      }
    }
    lgcop_n1 <- as.data.frame(lgcop_n1)
    colnames(lgcop_n1)<- "LogCopies(N1)"
    
    lgcop_n2 <- vector()
    for (i in 1:nrow(cq_n2)){
      print(i)
      if (is.na(cq_n2[i,])){
        lgcop_n2 <- append(lgcop_n2, NA)
      }else{
        t <- i*coeffs[1,2]+coeffs[1,1]
        lgcop_n2 <- append(lgcop_n2, t)
      }
    }
    lgcop_n2 <- as.data.frame(lgcop_n2)
    colnames(lgcop_n2)<- "LogCopies(N2)"
    
    #### Copies column
    cop_n1 <- vector()
    for (d in 1:nrow(lgcop_n1)){
      if (is.na(lgcop_n1[d,])){
        cop_n1 <- append(cop_n1, NA)
      }else{
        v <- 10^lgcop_n1[d,]
        cop_n1 <- append(cop_n1, v)
      }
    }
    cop_n1 <- as.data.frame(cop_n1)
    colnames(cop_n1)<- "Copies(N1)"
    
    #### Copies column
    cop_n2 <- vector()
    for (d in 1:nrow(lgcop_n2)){
      if (is.na(lgcop_n2[d,])){
        cop_n2 <- append(cop_n2, NA)
      }else{
        v <- 10^lgcop_n2[d,]
        cop_n2 <- append(cop_n2, v)
      }
    }
    cop_n2 <- as.data.frame(cop_n2)
    colnames(cop_n2)<- "Copies(N2)"
    
    final <- as.data.frame(cbind(a, wells_n1_def, wells_n2_def, wells_rnasep_def, cq_n1, cq_n2, cq_rnasep, lgcop_n1, lgcop_n2, cop_n1, cop_n2))
    return(final)
  })
  
  AnalysisSamplesDT <- reactive({
    df <- AnalysisSamples()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    
    f <- defdt %>%
      formatStyle(
        columns = 1:3,
        backgroundColor = "lemonchiffon"
      ) %>%
      formatStyle(
        columns = 4:6,
        backgroundColor = "khaki"
      )%>%
      formatStyle(
        columns = 7:9,
        backgroundColor = "lemonchiffon"
      )%>%
      formatStyle(
        columns = 10:11,
        backgroundColor = "khaki"
      )%>%
      formatStyle(
        columns = 12:13,
        backgroundColor = "lemonchiffon"
      )
    
    return(f)
    
  })
  
  output$analysis <- renderDataTable({
    AnalysisSamplesDT()
  })
  
  ###########################################################################
  ######################### APPLIED QUANT STUDIO ############################
  readApplied <- reactive({
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Raw Data")
    df <- data.frame(samples)
    df <- df[-c(1:44),1:4]
    colnames(df) <- c("Well", "Well_Position", "Cycle","Fluorescence")
    return(df)
  })
  
  output$prueba <- renderTable(
    readApplied()
  )
  
  aggApplied <- reactive({
    df <- readApplied()
    final <- df %>%
      pivot_wider(names_from = Cycle, values_from = Fluorescence)
    
    final <- as.data.frame(t(as.matrix(final)))
    f <- final %>%
      row_to_names(row_number = 2)
    return(f)
  })
  
  output$trans <- renderTable(
    aggApplied()
  )
  
  getGenes <- reactive({
    all <- aggApplied()
    
    n1 <- vector()
    n2 <- vector()
    rnasep <- vector()
    for (let in LETTERS[1:16]){
      n1 <- append(n1,rep(paste(let,1:8,sep="")))
      n2 <- append(n2,rep(paste(let,9:16,sep="")))
      rnasep <- append(rnasep,rep(paste(let,17:24,sep="")))
    }
    
    allgenes <- list(n1,n2,rnasep)
    
    n1_nam <- match(unlist(allgenes[1]), names(all))
    n1_def <- all[,n1_nam]
    n1_def <- cbind(Cycle = 1:40, n1_def)
    
    n2_nam <- match(unlist(allgenes[2]), names(all))
    n2_def <- all[,n2_nam]
    n2_def <- cbind(Cycle = 1:40, n2_def)
    
    rnasep_nam <- match(unlist(allgenes[3]), names(all))
    rnasep_def <- all[,rnasep_nam]
    rnasep_def <- cbind(Cycle = 1:40, rnasep_def)
    
    final <- list(n1_def, n2_def, rnasep_def)
    return(final)
  })
  
  output$transN1 <- renderTable(
    getGenes()[1]
  )
  
  ## Download plots
  output$downtransN1 <- downloadHandler(
    filename = "N1.csv",
    content = function(fname){
      write.csv(getGenes()[1], fname, quote = F, row.names = F)}
  )
  
  output$transN2 <- renderTable(
    getGenes()[2]
  )
  
  output$downtransN2 <- downloadHandler(
    filename = "N2.csv",
    content = function(fname){
      write.csv(getGenes()[2], fname, quote = F, row.names = F)}
  )
  
  output$transRNAsep <- renderTable(
    getGenes()[3]
  )
  
  output$downtransRNAseP <- downloadHandler(
    filename = "RNAseP.csv",
    content = function(fname){
      write.csv(getGenes()[3], fname, quote = F, row.names = F)}
  )
  
  ### Run Info ###
  readAppliedRunInfo <- reactive({
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Results")
    df <- data.frame(samples)
    df <- df[1:44,]
    return(df)
  })
  
  output$appliedruninfo <- renderTable({
    readAppliedRunInfo()
  })
  
  ## Read Applied Results sheet ##
  readAppliedResults <- reactive({
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Results")
    df <- data.frame(samples)
    df <- df[-c(1:44),c(1,2,4,7,9)]
    colnames(df) <- c("Well", "Well_Position","Sample", "Fluor","Cq")
    return(df)
  })
  
  output$appliedres <- renderTable({
    readAppliedResults()
  })
  
  #### Cq values Plate for Applied ###
  cqPlateApp <- reactive({
    inp <- readAppliedResults()
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$Cq[s:e]
      s <- s + 24
      e <- e + 24
    }
    dt <- dt %>% 
      mutate_if(is.numeric, round, digits = 3)
    return(dt)
  })
  
  ## Convert to datatable to render nicely in shiny
  printCqPlateApp <- reactive({
    dt <- cqPlateApp()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:8,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 9:16,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 17:24,
        backgroundColor = "lightblue")
    return(f)
  })
  
  output$cqplateapp<- DT::renderDataTable(
    printCqPlateApp()
  )
  
  ############# Sample Plate for Applied ##########################
  samplePlateApp <- reactive({
    inp <- readAppliedResults()
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$Sample[s:e]
      s <- s + 24
      e <- e + 24
    }
    return(dt)
  })
  ## Convert to datatable to render nicely in shiny
  printSamplePlateApp <- reactive({
    dt <- samplePlateApp()
    defdt <- datatable(dt, rownames = F,  options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:8,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 9:16,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 17:24,
        backgroundColor = "lightblue")
    return(f)
  })
  
  output$sampleplateapp<- renderDataTable(
    printSamplePlateApp()
  )
  
  #################### Check Sample plates for Applied #######################
  checkSamplesApp <- reactive({
    s <- samplePlateApp()
    first <- data.frame(d = unlist(s[1:8], use.names = FALSE))
    second <- data.frame(d = unlist(s[9:16], use.names = FALSE))
    third <- data.frame(d = unlist(s[17:24], use.names = FALSE))
    all <- cbind(first, second, third)
    colnames(all) <- c("first", "second", "third")
    #def <- head(all, 120)
    def <- all[complete.cases(all),]
    
    ## Add "coinciden" to script 
    for (i in 1:nrow(def)){
      if (def[i,"first"] == def[i, "second"] && def[i,"first"] == def[i, "third"]){
        def$Check[i] <- "Coinciden"
      } else {
        def$Check[i] <- "No coinciden"
      }
    }
    return(def)
  })
  
  ## Convert to datatable to render nicely in shiny
  checkSamplesAppDT <- reactive({
    def <- checkSamplesApp()
    defdt <- datatable(def, rownames = F,  options = list(pageLength = 60))
    f <- defdt %>% 
      formatStyle(
        columns = 1,
        backgroundColor = "yellow"
      ) %>%
      formatStyle(
        columns = 2,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 3,
        backgroundColor = "lightblue"
      ) %>%
      formatStyle(
        columns = 4,
        backgroundColor = styleEqual(c("Coinciden", "No coinciden"), c('seagreen', 'red'))
      )
    return(f)
  })
  
  ##### Render Table ####
  output$samplecheckapp <- renderDataTable(
    checkSamplesAppDT()
  )
  
  ################### Analysis Standard Curve for Applied #####################
  stCurveApp <- reactive({
    a <- data.frame(matrix(0, nrow = 5, ncol = 19))
    colnames(a)<-c("Sample","Dilution","Copies","log(copies)","N1_dup1","N2_dup1","RNAseP_dup1","N1_dup2","N2_dup2","RNAseP_dup2","N1_avg","N2_avg","RNAseP_avg","N1_logcop","N2_logcop","RNAseP_logcop", "N1_copies", "N2_copies", "RNAseP_copies")
    a$Sample <- c("NTC","C(-)","C(+)10-2","C(+)10-4","C(+)10-5")
    a$Dilution <- c("-","-",100,1000,100000)
    a$Copies <- c("-","-",200000*2/as.numeric(a$Dilution[3]), 200000*2/as.numeric(a$Dilution[4]),200000*2/as.numeric(a$Dilution[5]))
    a$`log(copies)`<- c("-","-", 3.602, 1.602, 0.602)
    
    cq <- cqPlateApp()
    a$N1_dup1 <- c(cq[9,8],'-',cq[15,8],cq[13,8], cq[11,8])
    a$N2_dup2 <- c(cq[9,16],'-',cq[15,16],cq[13,16], cq[11,16])
    a$RNAseP_dup1 <- c(cq[9,24],cq[15,24],'-','-','-')
    a$N1_dup2 <- c(cq[10,8],'-',cq[16,8],cq[14,8], cq[12,8])
    a$N2_dup2 <- c(cq[10,16],'-',cq[16,16],cq[14,16], cq[12,16])
    a$RNAseP_dup2 <- c(cq[10,24],cq[16,24],'-','-','-')
    a$N1_avg <- c(mean(c(cq[9,8],cq[10,8])), "-", mean(c(cq[15,8], cq[16,8])), mean(c(cq[13,8], cq[14,8])), mean(c(cq[11,8], cq[12,8])))
    a$N2_avg <- c(mean(c(cq[9,16],cq[10,16])), "-", mean(c(cq[15,16], cq[16,16])), mean(c(cq[13,16], cq[14,16])), mean(c(cq[11,16], cq[12,16])))
    a$RNAseP_avg <- c(mean(c(cq[9,24],cq[10,24])), "-", mean(c(cq[15,24], cq[16,24])), mean(c(cq[13,24], cq[14,24])), mean(c(cq[11,24], cq[12,24])))
    return(a)
  })
  
  #### Convert to datatable ####
  stCurveAppDT <- reactive({
    a <- stCurveApp()
    ## Transform into data table
    defdt <- datatable(a, rownames = F)
    dt <- defdt %>%
      formatStyle(
        columns = c(5,8,11),
        backgroundColor = "yellow"
      )%>%
      formatStyle(
        columns = c(6,9,12),
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = c(7,10,13),
        backgroundColor = "lightblue"
      )
    
    return(dt)
  })
  
  ##### Render Table ####
  output$stdcurveapp <- renderDataTable(
    stCurveAppDT()
  )
  
  ################### Standard Curve Plots ###########################
  stdCoeffsApp <- function(){
    a <- stCurveApp()
    n1 <- a$N1_avg[3:5]
    n2 <- a$N2_avg[3:5]
    cp <- a$`log(copies)`[3:5]
    p <- as.data.frame(cbind(as.numeric(n1),as.numeric(n2) , as.numeric(cp)))
    colnames(p) <- c("n1", "n2", "cp")
    
    model_n1 <- lm(cp ~ n1, p)
    n1_coeff1 <- as.numeric(model_n1$coefficients[1])
    n1_coeff2 <- as.numeric(model_n1$coefficients[2])
    n1_coeff <- as.data.frame(cbind(n1_coeff1, n1_coeff2))
    colnames(n1_coeff) <- c("Intercept", "Slope")
    
    model_n2 <- lm(cp ~ n2, p)
    n2_coeff1 <- as.numeric(model_n2$coefficients[1])
    n2_coeff2 <- as.numeric(model_n2$coefficients[2])
    n2_coeff <- as.data.frame(cbind(n2_coeff1, n2_coeff2))
    colnames(n2_coeff) <- c("Intercept", "Slope")
    
    coeff <- as.data.frame(rbind(n1_coeff, n2_coeff))
    rownames(coeff) <- c("N1", "N2")
    
    return(coeff)
  }
  
  stdPlotsApp <- reactive({
    a <- stCurveApp()
    
    n11 <- a$N1_avg[3:5]
    n22 <- a$N2_avg[3:5]
    cpp <- a$`log(copies)`[3:5]
    p <- as.data.frame(cbind(as.numeric(n11),as.numeric(n22),as.numeric(cpp)))
    colnames(p) <- c("n1", "n2", "cp")
    
    p1 <- ggplot(p, aes(x=n1, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = cp ~ n1,
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(p)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(p)[1]))) +
      ylab("log(Copies)")
    
    p2 <- ggplot(p, aes(x=n2, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = cp ~ n2, 
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(p)[2])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(p)[2]))) +
      ylab("log(Copies)")
    
    s <- grid.arrange(p1, p2, nrow = 1)
    return(s)
    
  })
  
  output$stdapp <- renderPlot(
    stdPlotsApp()
  )
  
  ############## ID Well for Applied ###############
  IDWELLtabApp <- reactive({
    inp <- readAppliedResults()
    wellid <- cbind(inp$well, as.character(inp$Target), inp$Sample)
    colnames(wellid) <- c("Well", "Target", "Sample")
    return(wellid)
  })
  
  output$IDWELLApp <- renderTable(
    IDWELLtabApp()
  )
  
  
}
  