################################# SHINY APP FOR QPCR ANALYSIS ############################################
########################### Sonia Olaechea-Lázaro (UPV/EHU, May 2020) ############################################  

server <- function(input, output) {
  
  ################## BIORAD ##################
  ############################################
  
  
  ################## Empty Plate Tab: Biorad #####################
  platesTable <- function(cnames, controls){
    col1 <- rep(1:8, each=2)
    col2 <- rep(9:16, each=2)
    col3 <- rep(17:24, each = 2)
    col4 <- rep(25: 32, each=2)
    col5 <- rep(33:40, each=2)
    col6 <- rep(41:48, each=2)
    col7 <- rep(49:56, each=2)
    col8 <- c(rep(57:59, each=2),controls)
    
    df <- cbind(col1,col2,col3,col4,col5,col6,col7, col8)
    rownames(df) <- LETTERS[1:16]
    colnames(df) <- cnames
    return(df)
  }
  
  platePrint <- function(){
    p1 <- platesTable(1:8, c("H2O","H2O","EXP-5","EXP-5","EXP-4","EXP-4","EXP-3","EXP-3","EXP-2","EXP-2"))
    p2 <- platesTable(9:16, c("H2O","H2O","EXP-5","EXP-5","EXP-4","EXP-4", "EXP-3","EXP-3","EXP-2","EXP-2"))
    p3 <- platesTable(17:24, c("H2O","H2O"," "," "," "," "," "," ", "HeLa","HeLa"))
    
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
  }
  
  output$plate <- renderDataTable(
    platePrint()
  )
  
  ################## Read "Run Information" and "Data" Sheet: Biorad ##########################
  # This function allows the user to upload either 1 Excel with 2 sheets (Data and Run Information)
  # or to upload 2 independent CSVs with the Run Information and Quantification Summary (with Sample,
  # Fluor, Target, Cq, etc columns)
  readBiorad <- function(){
    res <- input$biorad
    tbl_list <- list()
    ## Read data from CSVs
    if (grepl("csv",res$datapath[1]) == TRUE){
      tbl_list <- lapply(res$datapath, read.csv)
      tmp <- as.data.frame(tbl_list[1]) %>%
        select(Well, Fluor, Target, Content,Sample, Cq)
      colnames(tmp) <- c("Well", "Fluor", "Target", "Content", "ID", "Cq")
      tbl_list[[1]] <- tmp
      
      ## Read Data from 1 Excel
      } else if (grepl("xlsx", res$datapath[1]) == TRUE){
      dat <- read_xlsx(res$datapath, sheet = "Data")
      df <- as.data.frame(dat)
      df <- cbind(df$Well, df$Fluor, df$Target, df$Content ,df$Sample, df$Cq)
      df <- as.data.frame(df)
      colnames(df) <- c("Well", "Fluor", "Target", "Content","Sample", "Cq")
      
      ## Read Run Information
      run <- read_xlsx(res$datapath, sheet = "Run Information")
      df2 <- as.data.frame(run)
      
      ## Append to list
      tbl_list[[1]] <- df
      tbl_list[[2]] <- df2
    }
    return(tbl_list)
  }
  
  output$BRruninfo <- renderTable(
    readBiorad()[[2]]
  )
  
  output$inputdf <- renderTable(
    readBiorad()[[1]]
  )
  
  ################## New_Cq Tab: Biorad ##################
  newDataBiorad <- function(){
    inp <- readBiorad()[[1]]
    dat <- inp
    dat$Sample <- as.character(dat$Sample)
    for (i in 1:length(dat$Sample)){
      ## Interm column
      if (grepl("o", dat$Sample[i]) == TRUE){
        dat$interm <- dat$Sample[i]
      } else{
        dat$interm <- ""
      }
      ## NewSample column
      if (dat$Sample[i] == dat$interm[i]){
        dat$news <- ""
      } else {
        dat$news <- dat$Sample[i]
      }
      ## Formula column
      form <- rep(c(0,16,32,48,64,80,96,112), times = length(dat$Sample)/8)
      dat$Formula <- as.character(form)
      ## NewCq column
      if (dat$Cq[i] == "Undetermined"){
        dat$newcq[i] <- as.character(0)
      } else if (is.na(dat$Cq[i]) == TRUE){
        dat$newcq[i] <- as.character(0)
      } else if (dat$Cq[i] == ""){
        dat$newcq[i] <- as.character(0)
      } else {
        dat$newcq[i] <- dat$Cq[i]
      }
    }
    return(dat)
  }
  
  output$newdatabiorad <- renderTable({
    newDataBiorad()
  })
  
  ################## Cq Plate Tab: Biorad ######################
  cqPlate <- function(){
    #inp <- inputDF()
    inp <- readBiorad()[[1]]
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
  }
  
  ## Convert to datatable to render nicely in shiny
  printCqPlate <- function(){
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
  }
  
  #### Render Table in App ### 
  output$cqplate<- renderDataTable(
    printCqPlate()
  )
  
  ################## Sample Plate Tab: Biorad ############################
  samplePlate <- function(){
    #inp <- inputDF()
    inp <- readBiorad()[[1]]
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$ID[s:e]
      s <- s + 24
      e <- e + 24
    }
    return(dt)
  }
  
  ## Convert to datatable to render nicely in shiny ##
  printSamplePlate <- function(){
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
  }
  
  ###### Render Table in App
  output$sampleplate<- renderDataTable(
    printSamplePlate()
  )
  
  ################## Plate Setup MultiChanel: Biorad ########################
  setupMultiC <- function(){
    def <- checkSamples()
    col1 <- rep(paste("S",1:59, sep = ""),each =2)
    col2 <- rep(paste(col1, "-", 1:2, sep=""))
    col3 <- def[1:118,"first"]
    multic <- as.data.frame(cbind(col1, col2, col3))
    colnames(multic)<- c("Sample", "Replicate", "Real_ID")
    return(multic)
  }
  
  setupMultiCDT <- function(){
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
  }
  
  output$setupmultic <- renderDataTable(
    setupMultiCDT()
  )
  
  ################## Check Sample Order Tab: Biorad #######################
  checkSamples <- function(){
    s <- samplePlate()
    first <- data.frame(d = unlist(s[1:8], use.names = FALSE))
    second <- data.frame(d = unlist(s[9:16], use.names = FALSE))
    third <- data.frame(d = unlist(s[17:24], use.names = FALSE))
    third[ third == "HeLa"] <- NA
    all <- cbind(first, second, third)
    colnames(all) <- c("first", "second", "third")
    def <- all[complete.cases(all),]
    
    ## Add "coinciden" to script 
    for (i in 1:nrow(def)){
      if (def[i,"first"] == def[i, "second"] && def[i,"first"] == def[i, "third"]){
        def$Check[i] <- "Coinciden"
      } else {
        def$Check[i] <- "No coinciden"
      }
    }
    def <- head(def,59*2)
    return(def)
  }
  
  ## Convert to datatable to render nicely in shiny
  checkSamplesDT <- function(){
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
        backgroundColor = styleEqual(c("Coinciden", "No coinciden"), c('seagreen', 'tomato'))
      )
    return(f)
  }
  
  ####### Render Table in App
  output$samplecheck <- renderDataTable(
    checkSamplesDT()
  )
  
  ################## Standard Curve Tab (Table): Biorad #######################
  stCurve <- function(){
    a <- data.frame(matrix(0, nrow = 6, ncol = 19))
    colnames(a)<-c("Sample","Dilution","Copies","log(copies)","gen1_dup1","gen2_dup1","gen3_dup1","gen1_dup2","gen2_dup2","gen3_dup2","gen1_avg","gen2_avg","gen3_avg","gen1_logcop","gen2_logcop","gen3_logcop", "gen1_copies", "gen2_copies", "gen3_copies")
    a$Sample <- c("NTC","C(-)","C(+)10-2","C(+)10-3","C(+)10-4","C(+)10-5")
    a$Dilution <- c("-","-",100,1000,10000,100000)
    a$Copies <- c("-","-",200000*2/as.numeric(a$Dilution[3]), 200000*2/as.numeric(a$Dilution[4]),200000*2/as.numeric(a$Dilution[5]),200000*2/as.numeric(a$Dilution[6]))
    a$`log(copies)`<- c("-","-", 3.602, 2.602,1.602, 0.602)
    
    cq <- cqPlate()
    a$gen1_dup1 <- as.numeric(c(cq[7,8],NA,cq[15,8],cq[13,8], cq[11,8], cq[9,8]))
    a$gen2_dup1 <- as.numeric(c(cq[7,16],NA,cq[15,16],cq[13,16], cq[11,16], cq[9,16]))
    a$gen3_dup1 <- as.numeric(c(cq[9,24],cq[15,24],NA,NA,NA,NA))
    a$gen1_dup2 <- as.numeric(c(cq[8,8],NA,cq[16,8],cq[14,8], cq[12,8], cq[10,8]))
    a$gen2_dup2 <- as.numeric(c(cq[8,16],NA,cq[16,16],cq[14,16], cq[12,16], cq[10,16]))
    a$gen3_dup2 <- as.numeric(c(cq[10,24],cq[16,24],NA,NA,NA, NA))
    a$gen1_avg <- c(NA,NA,mean(c(a$gen1_dup1[3], a$gen1_dup2[3]), na.rm = T), mean(c(a$gen1_dup1[4], a$gen1_dup2[4]), na.rm = T), mean(c(a$gen1_dup1[5], a$gen1_dup2[5]), na.rm = T), mean(c(a$gen1_dup1[6], a$gen1_dup2[6]), na.rm = T))
    a$gen2_avg <- c(NA,NA,mean(c(a$gen2_dup1[3], a$gen2_dup2[3]), na.rm = T), mean(c(a$gen2_dup1[4], a$gen2_dup2[4]), na.rm = T), mean(c(a$gen2_dup1[5], a$gen2_dup2[5]), na.rm = T),mean(c(a$gen2_dup1[6], a$gen2_dup2[6]), na.rm = T))
    a$gen3_avg <- c(mean(c(a$gen3_dup1[1], a$gen3_dup2[1]), na.rm = T),mean(c(a$gen3_dup1[2], a$gen3_dup2[2]), na.rm = T),NA,NA,NA,NA)
    
    ################### Coefficients ####################
    gen1 <- a$gen1_avg[3:6]
    gen2 <- a$gen2_avg[3:6]
    cp <- a$`log(copies)`[3:6]
    p <- as.data.frame(cbind(as.numeric(gen1),as.numeric(gen2) , as.numeric(cp)))
    colnames(p) <- c("gen1", "gen2", "cp")
    
    model_gen1 <- lm(cp ~ gen1, p)
    gen1_coeff1 <- as.numeric(model_gen1$coefficients[1])
    gen1_coeff2 <- as.numeric(model_gen1$coefficients[2])
    gen1_coeff <- as.data.frame(cbind(gen1_coeff1, gen1_coeff2))
    colnames(gen1_coeff) <- c("Intercept", "Slope")
    
    model_gen2 <- lm(cp ~ gen2, p)
    gen2_coeff1 <- as.numeric(model_gen2$coefficients[1])
    gen2_coeff2 <- as.numeric(model_gen2$coefficients[2])
    gen2_coeff <- as.data.frame(cbind(gen2_coeff1, gen2_coeff2))
    colnames(gen2_coeff) <- c("Intercept", "Slope")
    
    coeff <- as.data.frame(rbind(gen1_coeff, gen2_coeff))
    rownames(coeff) <- c("gen1", "gen2")
    ###############################################################
    
    ## log(Copies) ##
    for (i in 1:length(a$gen1_avg)){
      if(is.na(a$gen1_avg[i])== TRUE){
        a$gen1_logcop[i] <- NA
      }else{
        a$gen1_logcop[i] <- a$gen1_avg[i]*coeff[1,2]+coeff[1,1]
      }
      if(is.na(a$gen2_avg[i])== TRUE){
        a$gen2_logcop[i] <- NA
      }else{
        a$gen2_logcop[i] <- a$gen2_avg[i]*coeff[2,2]+coeff[2,1]
      }
      if(is.na(a$gen3_avg[i])== TRUE){
        a$gen3_logcop[i] <- NA
      }else{
        a$gen3_logcop[i] <- a$gen3_avg[i]*coeff[2,2]+coeff[2,1]
      }
    }
    
    ##### Copies #####
    for (i in 1:length(a$gen1_logcop)){
      if (is.na(a$gen1_logcop[i])){
        a$gen1_copies[i] <- NA
      }else{
        a$gen1_copies[i] <-10^a$gen1_logcop[i]
      }
      if (is.na(a$gen2_logcop[i])){
        a$gen2_copies[i] <- NA
      }else{
        a$gen2_copies[i] <-10^a$gen2_logcop[i]
      }
      if (is.na(a$gen3_logcop[i])){
        a$gen3_copies[i] <- NA
      }else{
        a$gen3_copies[i] <-10^a$gen3_logcop[i]
      }
    }
    return(a)
  }
  
  #### Convert to datatable
  stCurveDT <- function(){
    a <- stCurve()
    ## Transform into data table
    defdt <- datatable(a, rownames = F)
    dt <- defdt %>%
      formatStyle(
        columns = c(5,8,11,14,17),
        backgroundColor = "yellow"
      )%>%
      formatStyle(
        columns = c(6,9,12,15,18),
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = c(7,10,13,16,19),
        backgroundColor = "lightblue"
      )
    
    return(dt)
  }
  
  ##### Render Table
  output$stdcurve <- renderDataTable(
    stCurveDT()
  )
  
  ####### Coefficients function ########
  stdCoeffs <- function(){
    a <- stCurve()
    gen1 <- a$gen1_avg[3:6]
    gen2 <- a$gen2_avg[3:6]
    cp <- a$`log(copies)`[3:6]
    p <- as.data.frame(cbind(as.numeric(gen1),as.numeric(gen2) , as.numeric(cp)))
    colnames(p) <- c("gen1", "gen2", "cp")
    
    model_gen1 <- lm(cp ~ gen1, p)
    gen1_coeff1 <- as.numeric(model_gen1$coefficients[1])
    gen1_coeff2 <- as.numeric(model_gen1$coefficients[2])
    gen1_coeff <- as.data.frame(cbind(gen1_coeff1, gen1_coeff2))
    colnames(gen1_coeff) <- c("Intercept", "Slope")
    
    model_gen2 <- lm(cp ~ gen2, p)
    gen2_coeff1 <- as.numeric(model_gen2$coefficients[1])
    gen2_coeff2 <- as.numeric(model_gen2$coefficients[2])
    gen2_coeff <- as.data.frame(cbind(gen2_coeff1, gen2_coeff2))
    colnames(gen2_coeff) <- c("Intercept", "Slope")
    
    coeff <- as.data.frame(rbind(gen1_coeff, gen2_coeff))
    rownames(coeff) <- c("gen1", "gen2")
    return(coeff)
  }
  
  ################## Standard Curve Tab (Plots): Biorad ##########################
  stdPlots <- function(){
    a <- stCurve()
    
    gen11 <- a$gen1_avg[3:6]
    gen22 <- a$gen2_avg[3:6]
    cpp <- a$`log(copies)`[3:6]
    
    pl1 <- as.data.frame(cbind(as.numeric(gen11),as.numeric(cpp)))
    colnames(pl1) <- c("gen1", "cp")
    pl1 <- pl1[complete.cases(pl1),]
    
    pl2 <- as.data.frame(cbind(as.numeric(gen22),as.numeric(cpp)))
    colnames(pl2) <- c("gen2", "cp")
    pl2 <- pl2[complete.cases(pl2),]
    
    form1 <- pl1$cp ~ pl1$gen1
    
    p1 <- ggplot(pl1, aes(x=gen1, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form1,
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl1)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl1)[1]))) +
      ylab("log(Copies)")
    
    form2 <- pl2$cp ~ pl2$gen2
    
    p2 <- ggplot(pl2, aes(x=gen2, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form2, 
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl2)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl2)[1]))) +
      ylab("log(Copies)")
    
    s <- grid.arrange(p1, p2, nrow = 1)
    return(s)
    
  }
  
  ##### Render Plots in App
  output$std <- renderPlot(
    stdPlots()
  )
  
  ################## Analysis Samples Tab: Biorad ###################
  AnalysisSamples <- function(){
    
    ## Real Id, Sample, Replicate columns
    a <- setupMultiC()
    
    ## Well columns
    wells_n1 <- vector()
    for (n in 1:8){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n1 <- append(wells_n1, p)
    }
    wells_n1_def <- as.data.frame(head(wells_n1, length(a[,1])))
    colnames(wells_n1_def) <- "Wells_N1"
    
    wells_n2 <- vector()
    for (n in 9:16){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n2 <- append(wells_n2, p)
    }
    wells_n2_def <- as.data.frame(head(wells_n2, length(a[,1])))
    colnames(wells_n2_def) <- "Wells_N2"
    
    wells_rnasep <- vector()
    for (n in 17:24){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_rnasep <- append(wells_rnasep, p)
    }
    wells_rnasep_def <- as.data.frame(head(wells_rnasep, length(a[,1])))
    colnames(wells_rnasep_def) <- "Wells_gen3"
    
    ###### Ct columns 
    cq <- cqPlate()
    cq_n1 <- stack(cq[,1:8])
    cq_n1 <- as.data.frame(head(cq_n1[,1], length(a[,1])))
    colnames(cq_n1) <- "Cq_N1"
    
    cq_n2 <- stack(cq[,9:16])
    cq_n2 <- as.data.frame(head(cq_n2[,1], length(a[,1])))
    colnames(cq_n2) <- "Cq_N2"
    
    cq_rnasep <- stack(cq[,17:24])
    cq_rnasep <- as.data.frame(head(cq_rnasep[,1], length(a[,1])))
    colnames(cq_rnasep) <- "Cq_RNAseP"
    
    #### log(copies) Column
    coeffs <- stdCoeffs()
    
    lgcop_n1 <- vector()
    for (i in 1:nrow(cq_n1)){
      if (is.na(cq_n1[i,])){
        lgcop_n1 <- append(lgcop_n1, NA)
      }else{
        t <- as.numeric(as.character(cq_n1[i,1]))*coeffs[1,2]+coeffs[1,1]
        lgcop_n1 <- append(lgcop_n1, t)
      }
    }
    lgcop_n1 <- as.data.frame(lgcop_n1)
    colnames(lgcop_n1)<- "LogCopies(N1)"
    
    lgcop_n2 <- vector()
    for (i in 1:nrow(cq_n2)){
      if (is.na(cq_n2[i,])){
        lgcop_n2 <- append(lgcop_n2, NA)
      }else{
        t <- as.numeric(as.character(cq_n2[i,1]))*(coeffs[2,2])+coeffs[2,1]
        lgcop_n2 <- append(lgcop_n2, t)
      }
    }
    lgcop_n2 <- as.data.frame(lgcop_n2)
    colnames(lgcop_n2)<- "LogCopies(N2)"
    
    #### Copies column N1
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
    
    #### Copies column N2
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
    
    ### Hidden columns
    c1 <- vector()
    c2 <- vector()
    c3 <- vector()
    c4 <- vector()
    c5 <- vector()
    c6 <- vector()
    c7 <- vector()
    c8 <- vector()
    c9 <- vector()
    for (i in 1:length(cq_rnasep[,1])){
      ## c1 col ##
      if(is.na(as.numeric(as.character(cq_rnasep[i,1]))) == TRUE){
        c1 <- append(c1,"OJO")
      }else if (cq_rnasep[i,1] > 35){
        c1 <- append(c1,"OJO")
      }else{
        c1 <- append(c1,"OK")
      }
      ## c2 col ##
      if(is.na(cop_n1[i,1]) == TRUE){
        c2 <- append(c2, "neg")
      } else if (cop_n1[i,1] > 4){
        c2 <- append(c2, "pos")
      }else{
        c2 <- append(c2, "neg")
      }
      ## c3 col ##
      if(is.na(cop_n2[i,1]) == TRUE){
        c3 <- append(c3, "neg")
      }else if (cop_n2[i,1] > 4){
        c3 <- append(c3, "pos")
      }else{
        c3 <- append(c3, "neg")
      }
      ## c4 col ##
      if(is.na(cq_n1[i,1]) == TRUE){
        c4 <- append(c4, "neg")
      }else if(cq_n1[i,1] < 40){
        c4 <- append(c4, "pos")
      }else{
        c4 <- append(c4, "neg")
      }
      ## c5 col ##
      if(is.na(cq_n2[i,1]) == TRUE){
        c5 <- append(c5, "neg")
      }else if(cq_n2[i,1] < 40){
        c5 <- append(c5, "pos")
      }else{
        c5 <- append(c5, "neg")
      }
      ## c6 col ##
      if (c2[i] == c4[i]){
        c6 <- append(c6, c2[i])
      }else if (c2[i] == "neg"){
        c6 <- append(c6, "dud")
      }else{
        c6 <- append(c6, "pos")
      }
      ## c7 col ##
      if (c3[i] == c5[i]){
        c7 <- append(c7, c3[i])
      }else if (c3[i] == "neg"){
        c7 <- append(c7, "dud")
      }else{
        c7 <- append(c7, "pos")
      }
      ## c8 col ##
      if (c6[i] == c7[i]){
        c8 <- append(c8, c6[i])
      }else if (c6[i] == "pos"){
        c8 <- append(c8, "pos")
      }else if (c7[i] == "pos"){
        c8 <- append(c8, "pos")
      }else{
        c8 <- append(c8, "dud")
      }
      ## c9 col ##
      if(c1[i] == "OJO"){
        c9 <- append(c9, "CI")
      }else{
        c9 <- append(c9, c8[i])
      }
    }

    c1 <- as.data.frame(c1)
    c2 <- as.data.frame(c2)
    c3 <- as.data.frame(c3)
    c4 <- as.data.frame(c4)
    c5 <- as.data.frame(c5)
    c6 <- as.data.frame(c6)
    c7 <- as.data.frame(c7)
    c8 <- as.data.frame(c8)
    c9 <- as.data.frame(c9)
    c9[] <- lapply(c9, as.character)
    
    indv <- vector()
    for (i in c9[,1]){
      if (i == "CI"){
        indv <- append(indv, "Calidad Insuficiente")
      } else if (i == "pos"){
        indv <- append(indv, "Positive")
      } else if (i == "neg"){
        indv <- append(indv, "Negative")
      } else if (i == "dud"){
        indv <- append(indv, "Dudosa")
      } else if (is.na(i) == TRUE){
        indv <- append(indv, "-")
      }
    }
    
    indv <- as.data.frame(indv)
    indv[] <- lapply(indv, as.character)
    
    final <- as.data.frame(cbind(a, wells_n1_def, wells_n2_def, wells_rnasep_def, cq_n1, cq_n2, cq_rnasep, lgcop_n1, lgcop_n2, cop_n1, cop_n2, c1,c2,c3,c4,c5,c6,c7,c8,c9, indv))
    return(final)
  }
  
  AnalysisSamplesDT <- function(){
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
  }
  
  output$analysis <- renderDataTable({
    AnalysisSamplesDT()
  })
  
  ################## Interpretation Tab: Biorad #########################
  interpretation <- function(){
    an <- AnalysisSamples()
    ID <- unique(an$Real_ID)
    c10 <- vector()
    s1 <- seq(from=1, to=117, by=2)
    s2 <- seq(from=2, to=118, by=2)
    for (i in 1:length(s1)){
      if(an$indv[s1[i]] == an$indv[s2[i]]){
        c10 <- append(c10, an$indv[s1[i]])
      }else{
        c10 <- append(c10, "Repeat")
      }
    }
    c10 <- as.data.frame(c10)
    c10[] <- lapply(c10, as.character)
    ## Interpretation ##
    interp <- vector()
    for (i in c10[,1]){
      if (i == "Calidad Insuficiente"){
        interp = append(interp, "Repeat")
      } else if (i == "Dudosa"){
        interp = append(interp, "Repeat")
      } else {
        interp = append(interp, i)
      }
    }
    
    interp <- as.data.frame(interp)
    interp[] <- lapply(interp, as.character)
    final <- cbind(ID, c10, interp)
    colnames(final) <- c("ID","", "Interpretation")
    return(final)
  }
  
  interpretationDT <- function(){
    df <- interpretation()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 60))
    
    f <- defdt %>%
      formatStyle(
        columns = 3,
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat", "Calidad Insuficiente"), 
                                     c('seagreen', 'tomato', 'lightblue', "lemonchiffon"))
      )
    return(f)
  }
  
  ##### Render in App
  output$interpret <- renderDataTable(
    interpretationDT()
  )
  
  
  ################## ID Well Tab: Biorad ##############################
  SamplesToKeep <- function(){
    inp <- readBiorad()[[1]]
    g <- inp %>%
      group_by(ID)
    
    g_sp <- group_split(g)
    tokeep <- vector()
    for (g in g_sp){
      if(all(is.na(g$Cq)) == TRUE){
        next
      } else if (all(is.na(g$Cq)) == FALSE){
        tokeep <- append(tokeep, as.numeric(as.character(unique(g$ID))))
      } 
    }
    tokeep <-na.exclude(tokeep)
    return(tokeep)
  }

  IDWELLtab <- function(){
    inp <- readBiorad()[[1]]
    wellid <- cbind(as.character(inp$Well), as.character(inp$Target), as.numeric(as.character(inp$ID)))
    wellid <- as.data.frame(wellid)
    colnames(wellid) <- c("Well", "Target", "ID")
    
    ## Get list of samples to keep
    tokeep <- SamplesToKeep()
    ## Keep only samples from wellid df
    wellid <- wellid[wellid$ID %in% tokeep,]
    return(wellid)
  }
  
  output$IDWELL <- renderTable(
    IDWELLtab()
  )
  
  output$downIDWELL <- downloadHandler(
    filename = "ID_well.csv",
    content = function(fname){
      write.csv(IDWELLtab(), fname, quote = F, row.names = F)}
  )
  
    
  ################## ID Result Tab: Biorad ###############################
  idresult <- function(){
    an <- AnalysisSamples()
    int <- interpretation()
    realid <- unique(as.numeric(as.character(an$Real_ID)))
    idres <- as.data.frame(cbind(realid, int$Interpretation))
    colnames(idres)<- c("ID", "Interpretation")
    
    tokeep <- SamplesToKeep()
    idres <- idres[idres$ID %in% tokeep,]
    
    return(idres)
  }
  
  output$downidres <- downloadHandler(
    filename = "ID_result.csv",
    content = function(fname){
      write.csv(idresult(), fname, quote = F, row.names = F)}
  )  
  
  output$IDRESULT <- renderTable(
    idresult()
  )
  
  
  
  
  ################## APPLIED QUANT STUDIO ##################
  ##########################################################
  
  ################## Read Applied "Raw Data" Sheet: Applied ###################
    readApplied <- function(){
        res <- input$appl
        tbl_list <- list()
        if (!is.null(res)){
        ## Data
        samples <- read_xlsx(res$datapath, sheet = "Raw Data")
        df <- as.data.frame(samples)
        df1 <- df[-c(1:44),1:4]
        colnames(df1) <- c("Well", "Well_Position", "Cycle","Fluorescence")
          
        # Run Information
        df2 <- samples[1:41,1:2]
        df2 <- as.data.frame(df2)
          
        # List of genes
        genes <- read_xlsx(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df3 <- df[-c(1:44),5]
        df3 <- unique(df3)
          
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Run Information
        tbl_list[[3]] <- df3 # List of genes
          
        return(tbl_list)
        }
  }
  ### Tab Names
  output$gen1app <- renderText({
    readApplied()[[3]][1]
  })
  
  output$gen2app <- renderText({
    readApplied()[[3]][2]
  })
  
  output$gen3app <- renderText({
    readApplied()[[3]][3]
  })
  
  ### Render in App
  output$readapp <- renderTable({
    readApplied()[[1]]
  })
  
  ################## Transposed Tab: Applied #####################
  aggApplied <- function(){
    df <- readApplied()[[1]]
    final <- df %>%
      pivot_wider(names_from = Cycle, values_from = Fluorescence)
    
    final <- as.data.frame(t(as.matrix(final)))
    f <- final %>%
      row_to_names(row_number = 2)
    return(f)
  }
  
  ##### Render in App 
  output$trans <- renderTable(
    aggApplied()
  )
  
  ################## N1, N2 and RNAseP Tabs: Applied ######################
  getGenes <- function(){
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
  }
  
  ###### Render N1 tab 
  output$transgen1 <- renderTable(
    getGenes()[1]
  )
  
  ## Download CSV ##
  output$downtransgen1 <- downloadHandler(
    filename = "N1.csv",
    content = function(fname){
      write.csv(getGenes()[1], fname, quote = F, row.names = F)}
  )
  
  ###### Render N2 tab 
  output$transgen2 <- renderTable(
    getGenes()[2]
  )
  
  ## Download CSV ##
  output$downtransgen2 <- downloadHandler(
    filename = "N2.csv",
    content = function(fname){
      write.csv(getGenes()[2], fname, quote = F, row.names = F)}
  )
  
  ###### Render RNAseP tab 
  output$transgen3 <- renderTable(
    getGenes()[3]
  )
  
  ## Download CSV ##
  output$downtransgen3 <- downloadHandler(
    filename = "RNAseP.csv",
    content = function(fname){
      write.csv(getGenes()[3], fname, quote = F, row.names = F)}
  )
  
  ################## Read Run Info: Applied ####################
  'readAppliedRunInfo <- function(){
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Results")
    df <- data.frame(samples)
    df <- df[1:44,]
    return(df)
  }'
  
  #### Render in App
  output$appliedruninfo <- renderTable({
    readApplied()[[2]]
  })
  
  ################## Read Results: Applied #####################
  readAppliedResults <- function(){
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Results", na = "Undetermined")
    df <- data.frame(samples)
    df <- df[-c(1:44),c(1,2,4,5,7,9)]
    colnames(df) <- c("Well", "Well_Position","Sample", "Target","Fluor","Cq")
    return(df)
  }
  
  output$appliedres <- renderTable({
    readAppliedResults()
  })
  
  ################## Cq Plate Tab: Applied ##################
  cqPlateApp <- function(){
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
  }
  
  ## Convert to datatable to render nicely in shiny
  printCqPlateApp <- function(){
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
  }
  
  output$cqplateapp<- renderDataTable(
    printCqPlateApp()
  )
  
  ################## Sample Plate Tab: Applied ##########################
  samplePlateApp <- function(){
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
  }
  
  ## Convert to datatable to render nicely in shiny
  printSamplePlateApp <- function(){
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
  }
  
  output$sampleplateapp<- renderDataTable(
    printSamplePlateApp()
  )
  
  ################## Check Sample Order Tab: Applied #######################
  checkSamplesApp <- function(){
    s <- samplePlateApp()
    first <- data.frame(d = unlist(s[1:8], use.names = FALSE))
    second <- data.frame(d = unlist(s[9:16], use.names = FALSE))
    third <- data.frame(d = unlist(s[17:24], use.names = FALSE))
    third[ third == "HeLa"] <- NA
    all <- cbind(first, second, third)
    colnames(all) <- c("first", "second", "third")
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
  }
  
  ## Convert to datatable to render nicely in shiny
  checkSamplesAppDT <- function(){
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
        backgroundColor = styleEqual(c("Coinciden", "No coinciden"), c('seagreen', 'tomato'))
      )
    return(f)
  }
  
  ##### Render Table
  output$samplecheckapp <- renderDataTable(
    checkSamplesAppDT()
  )
  
  ################## Plate Setup MultiChanel for Applied #################### 
  setupMultiCApp <- function(){
    def <- checkSamplesApp()
    
    realid <- as.numeric(as.character(def$first))
    sample <- rep(paste("S",1:(length(realid)/2), sep = ""),each=2)
    replic <- rep(paste(sample, "-", 1:2, sep=""))
    
    multic <- as.data.frame(cbind(sample,replic, realid))
    colnames(multic)<- c("Sample", "Replicate","Real_ID")
    return(multic)
  }
  
  setupMultiCAppDT <- function(){
    a <- setupMultiCApp()
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
  }
  
  output$setupmulticapp <- renderDataTable(
    setupMultiCAppDT()
  )
  
  ################## Standard Curve Tab (Table): Applied #####################
  stCurveApp <- function(){
    a <- data.frame(matrix(0, nrow = 5, ncol = 19))
    colnames(a)<-c("Sample","Dilution","Copies","log(copies)","N1_dup1","N2_dup1","RNAseP_dup1","N1_dup2","N2_dup2","RNAseP_dup2","N1_avg","N2_avg","RNAseP_avg","N1_logcop","N2_logcop","RNAseP_logcop", "N1_copies", "N2_copies", "RNAseP_copies")
    a$Sample <- c("NTC","C(-)","C(+)10-2","C(+)10-4","C(+)10-5")
    a$Dilution <- c("-","-",100,1000,100000)
    a$Copies <- c("-","-",200000*2/as.numeric(a$Dilution[3]), 200000*2/as.numeric(a$Dilution[4]),200000*2/as.numeric(a$Dilution[5]))
    a$`log(copies)`<- c("-","-", 3.602, 1.602, 0.602)
    
    cq <- cqPlateApp()
    
    a$N1_dup1 <- as.numeric(c(cq[9,8],NA,cq[15,8],cq[13,8], cq[11,8]))
    a$N2_dup1 <- as.numeric(c(cq[9,16],NA,cq[15,16],cq[13,16], cq[11,16]))
    a$RNAseP_dup1 <- as.numeric(c(cq[9,24],cq[15,24],NA,NA,NA))
    a$N1_dup2 <- as.numeric(c(cq[10,8],NA,cq[16,8],cq[14,8], cq[12,8]))
    a$N2_dup2 <- as.numeric(c(cq[10,16],NA,cq[16,16],cq[14,16], cq[12,16]))
    a$RNAseP_dup2 <- as.numeric(c(cq[10,24],cq[16,24],NA,NA,NA))
    a$N1_avg <- c(NA,NA,mean(c(a$N1_dup1[3], a$N1_dup2[3]), na.rm = T), mean(c(a$N1_dup1[4], a$N1_dup2[4]), na.rm = T), mean(c(a$N1_dup1[5], a$N1_dup2[5]), na.rm = T))
    a$N2_avg <- c(NA,NA,mean(c(a$N2_dup1[3], a$N2_dup2[3]), na.rm = T), mean(c(a$N2_dup1[4], a$N2_dup2[4]), na.rm = T), mean(c(a$N2_dup1[5], a$N2_dup2[5]), na.rm = T))
    a$RNAseP_avg <- c(NA,mean(c(a$RNAseP_dup1[2], a$RNAseP_dup2[2]), na.rm = T),mean(c(a$RNAseP_dup1[3], a$RNAseP_dup2[3]), na.rm = T), mean(c(a$RNAseP_dup1[4], a$RNAseP_dup2[4]), na.rm = T), mean(c(a$RNAseP_dup1[5], a$RNAseP_dup2[5]), na.rm =T))
    
    ################### Coefficients ####################
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
    ###############################################################
    ## log(Copies) ##
    for (i in 1:length(a$N1_avg)){
      if(is.na(a$N1_avg[i])== TRUE){
        a$N1_logcop[i] <- NA
      }else{
        a$N1_logcop[i] <- a$N1_avg[i]*coeff[1,2]+coeff[1,1]
      }
      if(is.na(a$N2_avg[i])== TRUE){
        a$N2_logcop[i] <- NA
      }else{
        a$N2_logcop[i] <- a$N2_avg[i]*coeff[2,2]+coeff[2,1]
      }
      if(is.na(a$RNAseP_avg[i])== TRUE){
        a$RNAseP_logcop[i] <- NA
      }else{
        a$RNAseP_logcop[i] <- a$RNAseP_avg[i]*coeff[2,2]+coeff[2,1]
      }
    }
    ## Copies ##
    for (i in 1:length(a$N1_logcop)){
      if (is.na(a$N1_logcop[i])){
        a$N1_copies[i] <- NA
      }else{
        a$N1_copies[i] <-10^a$N1_logcop[i]
      }
      if (is.na(a$N2_logcop[i])){
        a$N2_copies[i] <- NA
      }else{
        a$N2_copies[i] <-10^a$N2_logcop[i]
      }
      if (is.na(a$RNAseP_logcop[i])){
        a$RNAseP_copies[i] <- NA
      }else{
        a$RNAseP_copies[i] <-10^a$RNAseP_logcop[i]
      }
    }
    return(a)
  }
  
  #### Convert to datatable 
  stCurveAppDT <- function(){
    a <- stCurveApp()
    ## Transform into data table
    defdt <- datatable(a, rownames = F)
    dt <- defdt %>%
      formatStyle(
        columns = c(5,8,11,14,17),
        backgroundColor = "yellow"
      )%>%
      formatStyle(
        columns = c(6,9,12,15,18),
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = c(7,10,13,16,19),
        backgroundColor = "lightblue"
      )
    
    return(dt)
  }
  
  ##### Render Table
  output$stdcurveapp <- renderDataTable(
    stCurveAppDT()
  )
  
  ####### Coefficients function #######
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
  
  ################## Standard Curve Tab (Plots): Applied ###########################
  stdPlotsApp <- function(){
    a <- stCurveApp()
    n11 <- a$N1_avg[3:5]
    n22 <- a$N2_avg[3:5]
    cpp <- a$`log(copies)`[3:5]
    pl1 <- as.data.frame(cbind(as.numeric(n11),as.numeric(cpp)))
    colnames(pl1) <- c("n1", "cp")
    pl1 <- pl1[complete.cases(pl1),]
    
    pl2 <- as.data.frame(cbind(as.numeric(n22),as.numeric(cpp)))
    colnames(pl2) <- c("n2", "cp")
    pl2 <- pl2[complete.cases(pl2),]
    
    form1 <- pl1$cp ~ pl1$n1
    
    p1 <- ggplot(pl1, aes(x=n1, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form1,
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl1)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl1)[1]))) +
      ylab("log(Copies)")
    
    form2 <- pl2$cp ~ pl2$n2
    
    p2 <- ggplot(pl2, aes(x=n2, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form2, 
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl2)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl2)[1]))) +
      ylab("log(Copies)")
    
    s <- grid.arrange(p1, p2, nrow = 1)
    return(s)
  }
  
  output$stdapp <- renderPlot(
    stdPlotsApp()
  )
  
  ################## Analysis Samples Tab: Applied ######################
  AnalysisSamplesApp <- function(){
    
    ## Real Id, Sample, Replicate columns
    a <- setupMultiCApp()
    
    ## Well columns
    wells_n1 <- vector()
    for (n in 1:8){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n1 <- append(wells_n1, p)
    }
    wells_n1_def <- as.data.frame(head(wells_n1, length(a[,1])))
    colnames(wells_n1_def) <- "Wells_N1"
    
    wells_n2 <- vector()
    for (n in 9:16){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_n2 <- append(wells_n2, p)
    }
    wells_n2_def <- as.data.frame(head(wells_n2, length(a[,1])))
    colnames(wells_n2_def) <- "Wells_N2"
    
    wells_rnasep <- vector()
    for (n in 17:24){
      p <- rep(paste(LETTERS[1:16], n, sep=""))
      wells_rnasep <- append(wells_rnasep, p)
    }
    wells_rnasep_def <- as.data.frame(head(wells_rnasep, length(a[,1])))
    colnames(wells_rnasep_def) <- "Wells_RNAseP"
    
    ###### Ct columns 
    cq <- cqPlateApp()
    cq_n1 <- stack(cq[,1:8])
    cq_n1 <- as.data.frame(head(cq_n1[,1], length(a[,1])))
    colnames(cq_n1) <- "Cq_N1"
    
    cq_n2 <- stack(cq[,9:16])
    cq_n2 <- as.data.frame(head(cq_n2[,1], length(a[,1])))
    colnames(cq_n2) <- "Cq_N2"
    
    cq_rnasep <- stack(cq[,17:24])
    cq_rnasep <- as.data.frame(head(cq_rnasep[,1], length(a[,1])))
    colnames(cq_rnasep) <- "Cq_RNAseP"
    
    #### log(copies) Column
    coeffs <- stdCoeffsApp()
    
    lgcop_n1 <- vector()
    for (i in 1:nrow(cq_n1)){
      if (is.na(cq_n1[i,1])== TRUE){
        lgcop_n1 <- append(lgcop_n1, NA)
      }else{
        t <- as.numeric(as.character(cq_n1[i,1]))*coeffs[1,2]+as.numeric(coeffs[1,1])
        lgcop_n1 <- append(lgcop_n1, t)
      }
    }
    lgcop_n1 <- as.data.frame(lgcop_n1)
    colnames(lgcop_n1)<- "LogCopies(N1)"
    
    lgcop_n2 <- vector()
    for (i in 1:nrow(cq_n2)){
      if (is.na(cq_n2[i,1])== TRUE){
        lgcop_n2 <- append(lgcop_n2, NA)
      }else{
        t <- as.numeric(as.character(cq_n2[i,1]))*coeffs[1,2]+coeffs[1,1]
        lgcop_n2 <- append(lgcop_n2, t)
      }
    }
    lgcop_n2 <- as.data.frame(lgcop_n2)
    colnames(lgcop_n2)<- "LogCopies(N2)"
    
    #### Copies column N1
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
    
    #### Copies column N2
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
    
    ### Hidden columns
    c1 <- vector()
    c2 <- vector()
    c3 <- vector()
    c4 <- vector()
    c5 <- vector()
    c6 <- vector()
    c7 <- vector()
    c8 <- vector()
    c9 <- vector()
    
    for (i in 1:length(cq_rnasep[,1])){
      ## c1 col ##
      if(is.na(as.numeric(as.character(cq_rnasep[i,1]))) == TRUE){
        c1 <- append(c1,"OJO")
      }else if (as.numeric(as.character(cq_rnasep[i,1])) > 35){
        c1 <- append(c1,"OJO")
      }else{
        c1 <- append(c1,"OK")
      }
      
      ## c2 col ##
      if(is.na(cop_n1[i,1]) == TRUE){
        c2 <- append(c2, "neg")
      } else if (as.numeric(cop_n1[i,1]) > 4){
        c2 <- append(c2, "pos")
      }else{
        c2 <- append(c2, "neg")
      }
      ## c3 col ##
      if(is.na(cop_n2[i,1]) == TRUE){
        c3 <- append(c3, "neg")
      }else if (as.numeric(cop_n2[i,1]) > 4){
        c3 <- append(c3, "pos")
      }else{
        c3 <- append(c3, "neg")
      }
      ## c4 col ##
      if(is.na(cq_n1[i,1]) == TRUE){
        c4 <- append(c4, "neg")
      }else if(as.numeric(cq_n1[i,1]) < 40){
        c4 <- append(c4, "pos")
      }else{
        c4 <- append(c4, "neg")
      }
      ## c5 col ##
      if(is.na(cq_n2[i,1]) == TRUE){
        c5 <- append(c5, "neg")
      }else if(as.numeric(cq_n2[i,1]) < 40){
        c5 <- append(c5, "pos")
      }else{
        c5 <- append(c5, "neg")
      }
      ## c6 col ##
      if (c2[i] == c4[i]){
        c6 <- append(c6, c2[i])
      }else if (c2[i,1] == "neg"){
        c6 <- append(c6, "dud")
      }else{
        c6 <- append(c6, "pos")
      }
      ## c7 col ##
      if (c3[i] == c5[i]){
        c7 <- append(c7, c2[i])
      }else if (c3[i] == "neg"){
        c7 <- append(c7, "dud")
      }else{
        c7 <- append(c7, "pos")
      }
      ## c8 col ##
      if (c6[i] == c7[i]){
        c8 <- append(c8, c6[i])
      }else if (c6[i] == "pos"){
        c8 <- append(c8, "pos")
      }else if (c7[i] == "pos"){
        c8 <- append(c8, "pos")
      }else{
        c8 <- append(c8, "dud")
      }
      ## c9 col ##
      if(c1[i] == "OJO"){
        c9 <- append(c9, "CI")
      }else{
        c9 <- append(c9, c8[i])
      }
      ##
    }
    c1 <- as.data.frame(c1)
    c2 <- as.data.frame(c2)
    c3 <- as.data.frame(c3)
    c4 <- as.data.frame(c4)
    c5 <- as.data.frame(c5)
    c6 <- as.data.frame(c6)
    c7 <- as.data.frame(c7)
    c8 <- as.data.frame(c8)
    c9 <- as.data.frame(c9)
    
    c9[] <- lapply(c9, as.character)
    
    indv <- vector()
    for (i in c9[,1]){
      if (i == "CI"){
        indv <- append(indv, "Calidad Insuficiente")
      } else if (i == "pos"){
        indv <- append(indv, "Positive")
      } else if (i == "neg"){
        indv <- append(indv, "Negative")
      } else if (i == "dud"){
        indv <- append(indv, "Dudosa")
      } else if (is.na(i) == TRUE){
        indv <- append(indv, "-")
      }
    }
    indv <- as.data.frame(indv)
    indv[] <- lapply(indv, as.character)
    
    final <- as.data.frame(cbind(a, wells_n1_def, wells_n2_def, wells_rnasep_def, cq_n1, cq_n2, cq_rnasep, lgcop_n1, lgcop_n2, cop_n1, cop_n2, c1,c2,c3,c4,c5,c6,c7,c8,c9,indv))
    return(final)
  }
  
  AnalysisSamplesAppDT <- function(){
    df <- AnalysisSamplesApp()
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
  }
  
  output$analysisapp <- renderDataTable({
    AnalysisSamplesAppDT()
  })
  
  ################## Interpretation Tab: Applied ####################
  interpretationApp <- function(){
    an <- AnalysisSamplesApp()
    c10 <- vector()
    s1 <- seq(from=1, to=length(an[,1]), by=2)
    s2 <- seq(from=2, to=length(an[,1]), by=2)
    for (i in 1:length(s1)){
      if(an$indv[s1[i]] == an$indv[s2[i]]){
        c10 <- append(c10, an$indv[s1[i]])
      }else{
        c10 <- append(c10, "Repeat")
      }
    }
    c10 <- as.data.frame(c10)
    c10[] <- lapply(c10, as.character)
    ## Interpretation ##
    interp <- vector()
    for (i in c10[,1]){
      if (i == "Calidad Insuficiente"){
        interp = append(interp, "Repeat")
      } else if (i == "Dudosa"){
        interp = append(interp, "Repeat")
      } else {
        interp = append(interp, i)
      }
    }
    
    interp <- as.data.frame(interp)
    interp[] <- lapply(interp, as.character)
    final <- cbind(c10, interp)
    colnames(final) <- c("", "Interpretation")
    return(final)
  }
  
  
  interpretationAppDT <- function(){
    df <- interpretationApp()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 60))
    
    f <- defdt %>%
      formatStyle(
        columns = 2,
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat", "Calidad Insuficiente"), 
                                     c('seagreen', 'tomato', 'lightblue', "lemonchiffon"))
      )
    return(f)
  }
  
  ### Render in App ###
  output$interpretapp <- renderDataTable(
    interpretationAppDT()
  )
  
  ################## ID Well Tab: Applied ###############
  IDWELLtabApp <- function(){
    inp <- readAppliedResults()
    wellid <- cbind(inp$Well_Position, as.character(inp$Target), inp$Sample)
    colnames(wellid) <- c("Well", "Target", "Sample")
    return(wellid)
  }
  
  output$IDWELLApp <- renderTable(
    IDWELLtabApp()
  )
  
  output$downidwellapp <- downloadHandler(
    filename = "ID_well.csv",
    content = function(fname){
      write.csv(IDWELLtabApp(), fname, quote = F, row.names = F)}
  )
  
  ################## ID Result Tab: Applied ###################
  idresultApp <- function(){
    an <- AnalysisSamplesApp()
    int <- interpretationApp()
    realid <- unique(as.numeric(as.character(an$Real_ID)))
    final <- as.data.frame(cbind(realid, int$Interpretation))
    colnames(final)<- c("ID", "Interpretation")
    return(final)
  }
  
  output$downidresapp <- downloadHandler(
    filename = "ID_result.csv",
    content = function(fname){
      write.csv(idresultApp(), fname, quote = F, row.names = F)}
  )  
  
  output$IDRESULTApp <- renderTable(
    idresultApp()
  )
  
  
  ################## TAQMAN ##################
  ############################################
  
  ################## Read CSVs: Taqman ############################
  cyclesInput <- function(){
    dat <- input$taqmancsv
    if (!is.null(dat)){
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=",", dec=".")
      tbl_list <- lapply(tbl_list, function(x){x[1]<-NULL;x})
      return(tbl_list)
    }
  }
  
  ################## Gene List from input CSVs: Taqman ###############
  geneList <- function(){
    dat <- input$taqmancsv
    lst <- lapply(dat$name, FUN = function(x) gsub(pattern = ".*[_]([^.]+)[.].*", replacement = "\\1",
                                                   basename(dat$name)))
    lst <- unlist(lst)
    out <- unique(lst)
    def <- lapply(out, FUN = function(x) gsub(pattern = ".csv", replacement = "", out))
    out <- unique(unlist(def))
    
    return(out)
  }
  
  ################## ID Well Tab: Taqman ###################
  taqwellID <- function(){
    raw <- input$taqwellid
    if (!is.null(raw)){
      data = read.csv(raw$datapath, sep = ',', header = TRUE, dec=',');
      data[is.na(data)]<-""
      colnames(data)[3] <- "ID"
      data$Well <- as.character(data$Well)
      data <- as.data.frame(data)
      return(data)
    }else{
      return (NULL);
    }  
  }
  ## b) Call function and display in app
  output$taqmanwell <- renderTable({
    newWell = taqwellID()
    if(is.null(newWell)){
      return();
    }else{
      newWell;
    }
  })
  
  ################## ID Result Tab: Taqman ###################
  taqIDRes <- function(){
    idres <- input$taqidres
    if (!is.null(idres)){
      data <- read.csv(idres$datapath, sep = ',', header = TRUE, dec=',');
      data <- data[,-3] 
      return(data)
    }
  }
  
  ## Display table 
  output$taqmanidres <- renderTable({
    res <- taqIDRes()
    if(is.null(res)){
      return();
    }else{
      res;
    }
  })
  
  ################## Combine ID well and ID results: Taqman ###################
  combIDwellIDres <- function(){
    wellid <- taqwellID()
    idres <- taqIDRes()
    info <-merge(idres,wellid,by="ID")
    info$Well2<-as.character(sub('(?<![0-9])0*(?=[0-9])', '', info$Well, perl=TRUE))
    info$Target <- as.character(info$Target)
    info <-as.data.frame(info)
    return(info)
  }
  
  ################## C Tab: Taqman #################
  constVarendoC <- function(){
    info <- combIDwellIDres()
    C<-info[which(info$Target==input$endoC),] # select target
    C<-C[order(C$Well2),]
    return(C)
  }
  
  output$c <- renderTable({
    constVarendoC()
  })
  
  ################## Info Tab: Taqman ########
  output$info <- renderTable({
    inf <- combIDwellIDres()
    if(is.null(inf)){
      return();
    }else{
      inf;
    }
  })
  
  ################## Plot dimensions: Taqman ############
  taqDim <- function(){
    C <- constVarendoC()
    #Search for control in CSV filenames
    control <- as.character(unique(C$Target))
    csvs <- as.data.frame(input$taqmancsv)
    out <- csvs[grep(control, csvs$name), ]
    
    # Read control data for 
    raw <- read.csv(out$datapath, header = TRUE, sep = ",", dec=".")
    raw<-raw[,-1]
    raw<-raw[,C$Well2]
    l<-round_any(rowMaxs(as.matrix(raw),value=T), 100, f = ceiling)
    l <- max(l)
    r<-dim(raw)
    r<-r[2]
    r<-as.integer(r)
    lr <- list(l,r)
    return(lr)
  }
  
  ################## General Plots Tab: Taqman ##################
  generalPlots <- function(){
    info <- combIDwellIDres()
    lr <- taqDim()
    genes <- geneList()
    
    pltList <- list()
    for (i in 1:length(genes)){
      print(genes[i])
      df<-info[which(info$Target==genes[i]),] # select target
      df<-df[order(df$Well2),]
      
      ## load raw data (fluorescence-RFU per temperature) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ## Read fluorescence data for target gene
      inp <- cyclesInput()
      
      raw <- inp[i][[1]]
      #raw<-raw[,-1]
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
  
  ################## Indetermined Plots: Taqman #################
  indetPlots <- function(){
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
      #raw <- raw[,-1]
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
  }
  
  ## Render plots in app
  output$indetgen1 <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[1]], ncol = 2)
  })
  
  ## Download plots
  output$downlgen1 <- downloadHandler(
    filename = ("N1_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[1]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
    }
  )
  # N2 tab
  output$indetgen2 <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[2]], ncol = 2)
  })
  
  ## Download plots
  output$downlgen2 <- downloadHandler(
    filename = ("N2_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[2]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
    }
  )
  
  #RNAsep tab
  output$indetgen3 <- renderPlot({
    ls <- indetPlots()
    plot_grid(plotlist = ls[[3]], ncol = 2)
  })
  
  ## Download plots
  output$downlgen3 <- downloadHandler(
    filename = ("RNAseP_IndetPlots.pdf"),
    content = function(file){
      plots <- indetPlots()
      p <- plot_grid(plotlist = plots[[3]], ncol = 2)
      save_plot(file, p, ncol = 2, paper="a4")
    }
  )
  
  output$gen1 <- renderText({
    geneList()[1]
  })
  
  output$gen2 <- renderText({
    geneList()[2]
  })
  
  output$gen3 <- renderText({
    geneList()[3]
  })

  
  
  
  ################## SYBR ##################
  ##########################################
  
  
  
  ################## Read "Run Information" and "Data": SYBR ##########
  # This function allows the user to upload either 1 Excel with 2 sheets (Data and Run Information)
  # or to upload 2 independent CSVs with the Run Information and Quantification Summary (with Sample,
  # Fluor, Target, Cq, etc columns)
  readSYBR <- function(){
      res <- input$sybr
      tbl_list <- list()
      if (!is.null(res)){
      if (grepl("csv",res$datapath[[1]][1]) == TRUE){
        tbl_list <- lapply(res$datapath, read.csv)
        tmp <- as.data.frame(tbl_list[1]) %>%
          select(Well, Fluor, Target, Content,Sample, Cq)
        tbl_list[[1]] <- tmp
        
      } else if (grepl("xlsx", res$datapath[[1]][1]) == TRUE){
        ## Read Data
        dat <- read_xlsx(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Content ,df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target", "Content","Sample", "Cq")
        
        ## Read Run Information
        run <- read_xlsx(res$datapath, sheet = "Run Information")
        df2 <- as.data.frame(run)
        
        ## Append to list
        tbl_list[[1]] <- df
        tbl_list[[2]] <- df2
      }
      return(tbl_list)
    }
  }
  
  output$sybrruninfo <- renderTable(
    readSYBR()[2]
  )
  
  output$sybrdata <- renderTable(
    readSYBR()[1]
  )
  
  ################## New_Cq Tab: Biorad ##################
  newDataSYBR <- function(){
    inp <- readSYBR()[[1]]
    dat <- inp
    dat$Sample <- as.character(dat$Sample)
    for (i in 1:length(dat$Sample)){
      ## Interm column
      if (grepl("o", dat$Sample[i]) == TRUE){
        dat$interm <- dat$Sample[i]
      } else{
        dat$interm <- ""
      }
      ## NewSample column
      if (dat$Sample[i] == dat$interm[i]){
        dat$news <- ""
      } else {
        dat$news <- dat$Sample[i]
      }
      ## Formula column
      form <- rep(c(0,16,32,48,64,80,96,112), times = length(dat$Sample)/8)
      dat$Formula <- as.character(form)
      ## NewCq column
      if (dat$Cq[i] == "Undetermined"){
        dat$newcq[i] <- as.character(0)
      } else if (is.na(dat$Cq[i]) == TRUE){
        dat$newcq[i] <- as.character(0)
      } else if (dat$Cq[i] == ""){
        dat$newcq[i] <- as.character(0)
      } else {
        dat$newcq[i] <- dat$Cq[i]
      }
    }
    return(dat)
  }
  
  output$newdatasybr <- renderTable({
    newDataSYBR()
  })
  
  ################## Cq Plate Tab: SYBR ##############
  cqPlateSYBR <- function(){
    inp <- readSYBR()[[1]]
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
  }
  
  ## Convert to datatable to render nicely in shiny
  printCqPlateSYBR <- function(){
    dt <- cqPlateSYBR()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:6,
        backgroundColor = "lightgreen"
      ) %>%
      formatStyle(
        columns = 7:12,
        backgroundColor = "pink"
      ) %>%
      formatStyle(
        columns = 13:18,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 19:24,
        backgroundColor = "lightblue"
      )
    return(f)
  }
  
  #### Render Table in App ### 
  output$cqplatesybr<- renderDataTable(
    printCqPlateSYBR()
  )
  
  ################## Sample Plate Tab: SYBR ############################
  samplePlateSYBR <- function(){
    inp <- readSYBR()[[1]]
    dt <- data.frame(matrix(nrow = 16, ncol = 24))
    s <- 1
    e <- 24
    for (i in 1:16){
      dt[i,] <- inp$Sample[s:e]
      s <- s + 24
      e <- e + 24
    }
    return(dt)
  }
  
  ## Convert to datatable to render nicely in shiny ##
  printSamplePlateSYBR <- function(){
    dt <- samplePlateSYBR()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))
    f <- defdt %>% 
      formatStyle(
        columns = 1:6,
        backgroundColor = "lightgreen"
      ) %>%
      formatStyle(
        columns = 7:12,
        backgroundColor = "pink"
      ) %>%
      formatStyle(
        columns = 13:18,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 19:24,
        backgroundColor = "lightblue"
      )
    return(f)
  }
  
  ###### Render Table in App 
  output$sampleplatesybr<- renderDataTable(
    printSamplePlateSYBR()
  )
  
  ################## Check Sample Order Tab: SYBR #######################
  checkSamplesSYBR <- function(){
    s <- samplePlateSYBR()
    check1 <- rbind(s[1:16,][1:6])
    check2 <- rbind(s[1:16,][7:12])
    check3 <- rbind(s[1:16,][13:18])
    check4 <- rbind(s[1:16,][19:24])
    odd<-seq(1,15,2)
    even<-seq(2,16,2)
    #colnames(check1) <- c("Dup_1.1", "Dup_1.2", "Dup_2.1", "Dup_2.2", "Dup_3.1", "Dup_3.2", "Dup_4.1", "Dup_4.2")
    df  <- data.frame()
    for (i in 1:6){
      for (j in 1:8){
        row <- as.numeric(cbind(check1[odd[j],i],check1[even[j],i],check2[odd[j],i],check2[even[j],i],check3[odd[j],i],check3[even[j],i],check4[odd[j],i],check4[even[j],i]))
        df <- rbind(df,row)
      }
    }
    print(df)
    colnames(df) <- c("Dup_1.1", "Dup_1.2", "Dup_2.1", "Dup_2.2", "Dup_3.1", "Dup_3.2", "Dup_4.1", "Dup_4.2")
    def <- df[complete.cases(df),]
    
    ## Add "coinciden" to script 
    for (i in 1:nrow(def)){
      if (def[i,"Dup_1.1"] == def[i, "Dup_1.2"] && def[i,"Dup_1.1"] == def[i, "Dup_2.1"] && def[i, "Dup_1.1"] == def[i, "Dup_2.2"] && def[i,"Dup_1.1"] == def[i, "Dup_3.1"] && def[i,"Dup_1.1"] == def[i, "Dup_3.2"] && def[i,"Dup_1.1"] == def[i, "Dup_4.1"] && def[i,"Dup_1.1"] == def[i, "Dup_4.2"]){
        def$Check[i] <- "Coinciden"
      } else {
        def$Check[i] <- "No coinciden"
      }
    }
    return(def)
  }
  
  ## Convert to datatable to render nicely in shiny
  checkSamplesSYBRDT <- function(){
    def <- checkSamplesSYBR()
    defdt <- datatable(def, rownames = F,  options = list(pageLength = 60))
    f <- defdt %>% 
      formatStyle(
        columns = 1:2,
        backgroundColor = "lightgreen"
      )%>%
      formatStyle(
        columns = 3:4,
        backgroundColor = "pink"
      ) %>%
      formatStyle(
        columns = 5:6,
        backgroundColor = "orange"
      ) %>%
      formatStyle(
        columns = 7:8,
        backgroundColor = "lightblue"
      ) %>%
      formatStyle(
        columns = 9,
        backgroundColor = styleEqual(c("Coinciden", "No coinciden"), c('seagreen', 'tomato'))
      )
    return(f)
  }
  
  ####### Render Table in App
  output$samplechecksybr <- renderDataTable(
    checkSamplesSYBRDT()
  )
  
  
  ################## Plate Setup MultiChanel: SYBR ########################
  setupMultiCSYBR <- function(){
    def <- checkSamplesSYBR()
    col1 <- rep(paste("S",1:43, sep = ""),each=1)
    col2 <- unique(def$Dup_1.1)
    multic <- as.data.frame(cbind(col1, col2))
    colnames(multic)<- c("Sample", "Real_ID")
    print(multic)
    return(multic)
  }
  
  setupMultiCSYBRDT <- function(){
    a <- setupMultiCSYBR()
    defdt <- datatable(a, rownames = F, 
                       options = list(pageLength = 50))
    f <- defdt %>%
      formatStyle(
        columns = 1,
        backgroundColor = "seagreen"
      )%>%
      formatStyle(
        columns = 2,
        backgroundColor = "lightgreen"
      )
    return(f)
  }
  
  output$setupmulticsybr <- renderDataTable(
    setupMultiCSYBRDT()
  )
  
  ################## Standard Curve Tab (Table): SYBR #####################
  stCurveSYBR <- function(){
    a <- data.frame(matrix(0, nrow = 6, ncol = 24))
    colnames(a)<-c("Sample","Dilution","Copies","log(copies)","N_dup1","S_dup1","RdRP_dup1","RPP30_dup1","N_dup2","S_dup2","RdRP_dup2","RPP30_dup2","N_avg","S_avg","RdRP_avg","RPP30_avg","N_logcop","S_logcop","RdRP_logcop", "RPP30_logcop","N_copies", "S_copies", "RdRP_copies","RPP30_copies")
    a$Sample <- c("NTC","C(-)","C(+)10-2", "C(+)10-3","C(+)10-4","C(+)10-5")
    a$Dilution <- c("-","-",100,1000,10000,100000)
    a$Copies <- c("-","-",200000*2/as.numeric(a$Dilution[3]), 200000*2/as.numeric(a$Dilution[4]),200000*2/as.numeric(a$Dilution[5]),200000*2/as.numeric(a$Dilution[6]))
    a$`log(copies)`<- c("-","-", 3.602, 2.602,1.602, 0.602)
    
    cq <- cqPlateSYBR()
    
    a$N_dup1 <- as.numeric(c(cq[8,6],cq[16,24],cq[16,6],cq[14,6], cq[12,6], cq[10,6]))
    a$S_dup1 <- as.numeric(c(cq[8,12],cq[16,11],NA,NA,NA,NA))
    a$RdRP_dup1 <- as.numeric(c(cq[8,18],cq[16,13],cq[16,18],cq[14,18],cq[12,18],cq[10,18]))
    a$RPP30_dup1 <- as.numeric(c(cq[8,24],cq[16,15],NA,NA,NA,NA))
    a$N_dup2 <- as.numeric(c(cq[7,6],cq[15,24],cq[15,6],cq[13,6], cq[11,6], cq[9,6]))
    a$S_dup2 <- as.numeric(c(cq[7,12],cq[16,12],NA,NA,NA,NA))
    a$RdRP_dup2 <- as.numeric(c(cq[7,18],cq[16,14],cq[15,18], cq[13,18], cq[11,18],cq[9,18]))
    a$RPP30_dup2 <- as.numeric(c(cq[7,24],cq[16,16],NA,NA,NA,NA))
    a$N_avg <- c(mean(c(a$N_dup1[1], a$N_dup2[1]), na.rm = T), mean(c(a$N_dup1[2], a$N_dup2[2]), na.rm = T), mean(c(a$N_dup1[3], a$N_dup2[3]), na.rm = T), mean(c(a$N_dup1[4], a$N_dup2[4]), na.rm = T), mean(c(a$N_dup1[5], a$N_dup2[5]), na.rm = T), mean(c(a$N_dup1[6], a$N_dup2[6]), na.rm = T))
    a$S_avg <- c(mean(c(a$S_dup1[1], a$S_dup2[1]), na.rm = T), mean(c(a$S_dup1[2], a$S_dup2[2]), na.rm = T), NA, NA,NA,NA)
    a$RdRP_avg <- c(mean(c(a$RdRP_dup1[1], a$RdRP_dup2[1]), na.rm = T), mean(c(a$RdRP_dup1[2], a$RdRP_dup2[2]), na.rm = T), mean(c(a$RdRP_dup1[3], a$RdRP_dup2[3]), na.rm = T), mean(c(a$RdRP_dup1[4], a$N_dup2[4]), na.rm = T), mean(c(a$RdRP_dup1[5], a$RdRP_dup2[5]), na.rm = T), mean(c(a$RdRP_dup1[6], a$RdRP_dup2[6]), na.rm = T))
    a$RPP30_avg <- c(mean(c(a$RPP30_dup1[1], a$RPP30_dup2[1]), na.rm = T), mean(c(a$RPP30_dup1[2], a$RPP30_dup2[2]), na.rm = T), NA, NA,NA,NA)
    
    ################### Coefficients ####################
    n <- a$N_avg[3:6]
    s <- a$S_avg[3:6]
    rdrp <- a$RdRP_avg[3:6]
    rpp30 <- a$RPP30[3:6]
    cp <- a$`log(copies)`[3:6]
    p <- as.data.frame(cbind(as.numeric(n),as.numeric(s),as.numeric(rdrp),as.numeric(rpp30),as.numeric(cp)))
    #colnames(p) <- c("N", "S","RdRP", "RPP30","cp")
    
    model_n <- lm(cp ~ n, p)
    n_coeff1 <- as.numeric(model_n$coefficients[1])
    n_coeff2 <- as.numeric(model_n$coefficients[2])
    n_coeff <- as.data.frame(cbind(n_coeff1, n_coeff2))
    colnames(n_coeff) <- c("Intercept", "Slope")
    
    model_rdrp <- lm(cp ~ rdrp, p)
    rdrp_coeff1 <- as.numeric(model_rdrp$coefficients[1])
    rdrp_coeff2 <- as.numeric(model_rdrp$coefficients[2])
    rdrp_coeff <- as.data.frame(cbind(rdrp_coeff1, rdrp_coeff2))
    colnames(rdrp_coeff) <- c("Intercept", "Slope")
    
    
    coeff <- as.data.frame(rbind(n_coeff, rdrp_coeff))
    rownames(coeff) <- c("N", "RdRP")
    ###############################################################
    ## log(Copies) ##
    for (i in 1:length(a$N_avg)){
      if(is.na(a$N_avg[i])== TRUE){
        a$N_logcop[i] <- NA
      }else{
        a$N_logcop[i] <- a$N_avg[i]*coeff[1,2]+coeff[1,1]
      }
      if(is.na(a$RdRP_avg[i])== TRUE){
        a$RdRP_logcop[i] <- NA
      }else{
        a$RdRP_logcop[i] <- a$RdRP_avg[i]*coeff[2,2]+coeff[2,1]
      }
      
      a$S_logcop[i] <- NA
      a$RPP30_logcop[i] <- NA
    }
    ## Copies ##
    for (i in 1:length(a$N_logcop)){
      if (is.na(a$N_logcop[i])){
        a$N_copies[i] <- NA
      }else{
        a$N_copies[i] <-10^a$N_logcop[i]
      }
      if (is.na(a$RdRP_logcop[i])){
        a$RdRP_copies[i] <- NA
      }else{
        a$RdRP_copies[i] <-10^a$RdRP_logcop[i]
      }
      a$S_copies[i] <- NA
      a$RPP30_copies[i] <- NA
    }
    return(a)
  }
  
  # Convert to datatable 
  stCurveSYBRDT <- function(){
    a <- stCurveSYBR()
    ## Transform into data table
    defdt <- datatable(a, rownames = F)
    dt <- defdt %>%
      formatStyle(
        columns = c(5,9,13,17,21),
        backgroundColor = "lightgreen"
      )%>%
      formatStyle(
        columns = c(6,10,14,18,22),
        backgroundColor = "pink"
      ) %>%
      formatStyle(
        columns = c(7,11,15,19,23),
        backgroundColor = "lightblue"
      ) %>%
      formatStyle(
        columns = c(8,12,16,20,24),
        backgroundColor = "orange"
      )
    
    return(dt)
  }
  
  ##### Render Table
  output$stdcurvesybr <- renderDataTable(
    stCurveSYBRDT()
  )
  
  ################## Standard Curve Tab (Plots): SYBR ###########################
  stdPlotsSYBR <- function(){
    a <- stCurveSYBR()
    nn <- a$N_avg[3:6]
    rdrpp <- a$RdRP_avg[3:6]
    cpp <- a$`log(copies)`[3:6]
    pl1 <- as.data.frame(cbind(as.numeric(nn),as.numeric(cpp)))
    colnames(pl1) <- c("n", "cp")
    pl1 <- pl1[complete.cases(pl1),]
    
    pl2 <- as.data.frame(cbind(as.numeric(rdrpp),as.numeric(cpp)))
    colnames(pl2) <- c("rdrp", "cp")
    pl2 <- pl2[complete.cases(pl2),]
    
    form1 <- pl1$cp ~ pl1$n
    
    p1 <- ggplot(pl1, aes(x=n, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form1,
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                   parse = TRUE, na.rm = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl1)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl1)[1]))) +
      ylab("log(Copies)")
    
    form2 <- pl2$cp ~ pl2$rdrp
    
    p2 <- ggplot(pl2, aes(x=rdrp, y = cp)) + 
      geom_point()+
      geom_smooth(method = lm, se = F) +
      stat_poly_eq(formula = form2, 
                   label.x.npc = "right", label.y.npc = "top",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                   parse = TRUE, na.rm = TRUE) +
      ggtitle(paste("Standard curve for",toupper(colnames(pl2)[1])))+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
      xlab(paste("Dilutions",toupper(colnames(pl2)[1]))) +
      ylab("log(Copies)")
    
    s <- grid.arrange(p1, p2, nrow = 1)
    return(s)
  }
  
  output$stdsybr <- renderPlot(
    stdPlotsSYBR()
  )
  
  ################## Analysis Samples Tab: SYBR ######################
  AnalysisSamplesSYBR <- function(){
    
    ## Real Id, Sample, Replicate columns
    a <- setupMultiCSYBR()
    sample <- rep(paste("S",1:length(a$Real_ID), sep=""), each=1)
    dup1 <- paste(sample,"-1", sep="")
    dup2 <- paste(sample,"-2", sep="")
    
    cq <- cqPlateSYBR()
    n <- as.numeric(stack(cq[,1:6])[,1])

    odd<-seq(1,length(n)-1,2)
    even<-seq(2,length(n),2)
    
    n_dup1 <- vector()
    for (i in odd){
      n_dup1 <- append(n_dup1,n[i])
    }
    n_dup2 <- vector()
    for (i in even){
      n_dup2 <- append(n_dup2, n[i])
    }
    
    s <- as.numeric(stack(cq[,7:12])[,1])
    s_dup1 <- vector()
    for (i in odd){
      s_dup1 <- append(s_dup1,s[i])
    }
    s_dup2 <- vector()
    for (i in even){
      s_dup2 <- append(s_dup2,s[i])
    }
    
    rdrp <- as.numeric(stack(cq[,13:18])[,1])
    rdrp_dup1 <- vector()
    for (i in odd){
      rdrp_dup1 <- append(rdrp_dup1,rdrp[i])
    }
    rdrp_dup2 <- vector()
    for (i in even){
      rdrp_dup2 <- append(rdrp_dup2,rdrp[i])
    }
    
    rpp30 <- as.numeric(stack(cq[,19:24])[,1])
    rpp30_dup1 <- vector()
    for (i in odd){
      rpp30_dup1 <- append(rpp30_dup1,rpp30[i])
    }
    rpp30_dup2 <- vector()
    for (i in even){
      rpp30_dup2 <- append(rpp30_dup2,rpp30[i])
    }
    
    #### Controls ####
    h1 <- vector()
    h2 <- vector()
    n1 <- vector()
    n2 <- vector()
    s1 <- vector()
    s2 <- vector()
    r1 <- vector()
    r2 <- vector()
    n1n2 <- vector()
    s1s2 <- vector()
    r1r2 <- vector()
    spos <- vector()
    npos <- vector()
    rpos <- vector()
    score <- vector()
    assignment <- vector()
    comments <- vector()
    for (i in 1:nrow(a)){
      #### h1 ####
      if (is.na(rpp30_dup1[i]) == TRUE){
        h1 <- append(h1, "OJO")
      } else if (as.numeric(rpp30_dup1[i]) > 35){
        h1 <- append(h1, "OJO")
      } else {
        h1 <- append(h1, "OK")
      }
      #### h2 ####
      if (is.na(rpp30_dup2[i]) == TRUE){
        h2 <- append(h2, "OJO")
        }else if (as.numeric(rpp30_dup2[i]) > 35){
          h2 <- append(h2, "OJO")
        }else {
          h2 <- append(h2, "OK")
        }
     
      ### n1 ####
      if (is.na(n_dup1[i]) == TRUE){
        n1 <- append(n1, "null")
      } else if (as.numeric(n_dup1[i]) == 0){
        n1 <- append(n1, "neg")
      } else if (as.numeric(n_dup1[i]) < 40){
        n1 <- append(n1, "pos")
      } else {
        n1 <- append(n1, "neg")
      }
      
      ## n2 ##
      if (is.na(n_dup2[i]) == TRUE){
        n2 <- append(n2, "null")
      } else if (as.numeric(n_dup2[i]) == 0){
        n2 <- append(n2, "neg")
      } else if (as.numeric(n_dup2[i]) < 40){
        n2 <- append(n2, "pos")
      } else {
        n2 <- append(n2, "neg")
      }
      
      ## s1 ##
      if (is.na(s_dup1[i]) == TRUE){
        s1 <- append(s1, "null")
      } else if (as.numeric(s_dup1[i]) == 0){
        s1 <- append(s1, "neg")
      } else if (as.numeric(s_dup1[i]) < 40){
        s1 <- append(s1, "pos")
      } else {
        s1 <- append(s1, "neg")
      }
      
      ## s2 ##
      if (is.na(s_dup2[i]) == TRUE){
        s2 <- append(s2, "null")
      } else if (as.numeric(s_dup2[i]) == 0){
        s2 <- append(s2, "neg")
      } else if (as.numeric(s_dup2[i]) < 40){
        s2 <- append(s2, "pos")
      } else {
        s2 <- append(s2, "neg")
      }
      
      ## rdrp1 ##
      if (is.na(rdrp_dup1[i]) == TRUE){
        r1 <- append(r1, "null")
      } else if (as.numeric(rdrp_dup1[i]) == 0){
        r1 <- append(r1, "neg")
      } else if (as.numeric(rdrp_dup1[i]) < 40){
        r1 <- append(r1, "pos")
      } else {
        r1 <- append(r1, "neg")
      }
      
      ## rdrp2 ##
      if (is.na(rdrp_dup2[i]) == TRUE){
        r2 <- append(r2, "null")
      } else if (as.numeric(rdrp_dup2[i]) == 0){
        r2 <- append(r2, "neg")
      } else if (as.numeric(rdrp_dup2[i]) < 40){
        r2 <- append(r2, "pos")
      } else {
        r2 <- append(r2, "neg")
      }
      
      n1n2 <- append(n1n2, paste(n1[i],"N",n2[i], sep=""))
      s1s2 <- append(s1s2, paste(s1[i],"S",s2[i], sep=""))
      r1r2 <- append(r1r2, paste(r1[i],"R",r2[i], sep=""))
      
      ### spos ###
      if (s1s2[i] == "negSneg"){
        spos <- append(spos, 0)
      }else if(s1s2[i] == "posSpos"){
        spos <- append(spos, 1)
      } else {
        spos <- append(spos, 9)
      }
      
      ### npos ###
      if (n1n2[i] == "negNneg"){
        npos <- append(npos, 0)
      }else if(n1n2[i] == "posNpos"){
        npos <- append(npos, 1)
      } else {
        npos <- append(npos, 9)
      }
      
      ### npos ###
      if (r1r2[i] == "negNneg"){
        rpos <- append(rpos, 0)
      }else if(r1r2[i] == "posNpos"){
        rpos <- append(rpos, 1)
      } else {
        rpos <- append(rpos, 9)
      }
      
      ### Score ###
      score <- append(score, sum(npos[i],spos[i],rpos[i]))
      
      ### Assignment ###
      if (is.na(a$Real_ID[i]) == TRUE){
        assignment <- append(assignment, "")
      } else if(h1[i] == "OJO"){
        assignment <- append(assignment, "Null")
      } else if(h2[i] == "OJO"){
        assignment <- append(assignment, "Null")
      } else if (as.numeric(score[i]) < 2){
        assignment <- append(assignment, "Negative")
      } else if (as.numeric(score[i]) < 4){
        assignment <- append(assignment, "Positive")
      } else if (as.numeric(score[i]) == 10){
        assignment <- append(assignment, "Null")
      } else if (as.numeric(score[i]) < 10){
        assignment <- append(assignment, "Negative")
      } else if (as.numeric(score[i]) < 12){
        assignment <- append(assignment, "Positive")
      } else {
        assignment <- append(assignment, "Null")
      }
      
      if (is.na(a$Real_ID[i]) == TRUE){
        comments <- append(comments, "")
      } else if (h1[i] == "OJO") {
        comments <- append(comments, "Housekeeping failed")
      } else if (h2[i] == "OJO") {
        comments <- append(comments, "Housekeeping failed")
      } else {
        comments <- append(comments, "")
      }
      
    }
    
    final <- cbind(a$Real_ID, sample, dup1, n_dup1, s_dup1, rdrp_dup1, rpp30_dup1, dup2, n_dup2, s_dup2, rdrp_dup2, rpp30_dup2,h1,h2,n1,n2,r1,r2, n1n2,s1s2,r1r2, npos, spos, rpos, score, assignment, comments)
    final <- as.data.frame(final)
    colnames(final) <- c("Real_ID","Sample", "Dup1", "N_Ct_1", "S_Ct_1", "RdRP_Ct_1", "RPP30_Ct_1","Dup2","N_Ct_2", "S_Ct_2", "RdRP_Ct_2", "RPP30_Ct_2", "H1", "H2", "N1","N2","R1","R2", "N1N2", "S1S2", "R1R2", "Npos", "Spos", "Rpos", "Score", "Assignment", "Comments")
    return(final)
  }
  
  AnalysisSamplesSYBRDT <- function(){
    df <- AnalysisSamplesSYBR()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    f <- defdt %>%
      formatStyle(
        columns = 1:2,
        backgroundColor = "lemonchiffon"
      )
    return(f)
  }
  
  output$analysissybr <- renderDataTable({
    AnalysisSamplesSYBRDT()
  })
  
  ################## ID Well Tab: SYBR ##############################
  IDWELLtabSYBR <- function(){
    inp <- readSYBR()[[1]]
    wellid <- cbind(as.character(inp$Well), as.character(inp$Target), as.character(inp$Sample))
    wellid <- as.data.frame(wellid)
    colnames(wellid) <- c("Well", "Target", "Sample")
    out <- wellid[order(as.numeric(as.character(wellid$Sample))),]
    out[out == ""] <- NA
    def <- out[complete.cases(out),]
    def <- def[!grepl("HeLa", def$Sample),]
    def <- def[!grepl("10-*", def$Sample),]
    return(def)
  }
  
  output$downIDWELLsybr <- downloadHandler(
    filename = "ID_well.csv",
    content = function(fname){
      write.csv(IDWELLtabSYBR(), fname, quote = F, row.names = F)}
  )
  output$IDWELLsybr <- renderTable(
    IDWELLtabSYBR()
  )
  
  ################## ID Result Tab: SYBR ###################
  idresultSYBR <- function(){
    an <- AnalysisSamplesSYBR()
    realid <- as.numeric(as.character(an$Real_ID))
    final <- as.data.frame(cbind(realid, as.character(an$Assignment)))
    colnames(final)<- c("ID", "Interpretation")
    return(final)
  }
  
  output$downIDRESsybr <- downloadHandler(
    filename = "ID_result.csv",
    content = function(fname){
      write.csv(idresultSYBR(), fname, quote = F, row.names = F)}
  )  
  
  output$IDRESsybr <- renderTable(
    idresultSYBR()
  )
  
  
  ################### SYBR MELTING CURVE ANALYSIS ###################
  
  ################## Read CSVs with Tm: SYBR ######### 
  TMinput <- function(){
    dat <- input$sybrcsv
    if (!is.null(dat)){
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=",", dec=",")
      tbl_list <- lapply(tbl_list, function(x){x[1]<-NULL;x})
      return(tbl_list)
    }
  }
  
  ################## Gene List: SYBR ###########
  SybrGeneList <- function(){
    dat <- input$sybrcsv
    lst <- lapply(dat$name, FUN = function(x) gsub(pattern = ".*[_]([^.]+)[.].*", replacement = "\\1",
                                                   basename(dat$name)))
    lst <- unlist(lst)
    out <- unique(lst)
    return(out)
  }
  
  ################## Well ID: SYBR ###########
  wellID <- function(){
    raw <- input$wellid
    if (!is.null(raw)){
      data = read.csv(raw$datapath, sep = ',', header = TRUE, dec=',');
      ## Correct Well names, adds new column called Well2 to data 
      data$Well2<-sub('(^[A-Z]+)0', "\\1", data$Well, perl=TRUE) #replace f.i. A01 by A1
      return(data)
    }else{
      return (NULL);
    }  
  }
  
  ## b) Call function and display in app
  output$well <- renderTable({
    newWell = wellID()
    if(is.null(newWell)){
      return();
    }else{
      newWell;
    }
  })

  
  ################## Fluos_Gene Tab: SYBR ###############
  ##### Select target gene and match columns in fluorescence file
  ## matchTarget por archivo independiente
  matchTarget <- function(tm, gene){
    Well <- wellID()
    #print(Well[which(Well$Target==gene),])
    well_gene <- Well[which(Well$Target==gene),] # select target
    #Select fluorescence matching columns
    fluos_gene<-tm[,well_gene$Well2]
    
    return(fluos_gene)
  }
  
  ## Misma función que matchTarget pero con un loop para guardar todos los archivos en la misma tabla
  ## Bastante guarro, pero funciona
  matchAllTarget <- function(){
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
  }
  
  output$fluosgene <- renderTable({
    matchAllTarget()
  })
  
  '############ Fluo Temp Tab: SYBR ###############
  fluoTemp <- function(){
    fluos_gene <- matchAllTarget()
    melt_gene <- TMinput()
    genes <- SybrGeneList()
    LsforPlot <- list()
    for (i in 1:length(genes)){
      fluo_temp <- cbind(melt_gene[[i]]$Temperature, fluos_gene[[i]])
      #fluo_temp <- fluos_gene[[i]]
      #colnames(fluo_temp)[1] <- "Temperature"
      LsforPlot[[i]] <- fluo_temp
    }
    return(LsforPlot)
  }
  
  output$fluotemp <-renderTable({
    fluoTemp()
  })'
  
  ################## Melting Curve Plots: SYBR ##################
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
    c<-as.integer(b+1)
    d<-as.integer(2*b)
    data_gene <- cbind(fluos_gene,rep.col(tm$Temperature,b))
    f <- function(df){
      out <- as.character(as.numeric(df))
    }
    #data_gene <- sapply(data_gene, f)
    #print(str(data_gene))
    data_gene <- as.data.frame(data_gene)
    
    if (input$isderiv == TRUE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5),is.deriv=T)
    } else if (input$isderiv == FALSE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=10,Tm.border=c(0.5,0.5))
    }
  }
  
  ##### b) Render plots for each gene in a different tab
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
  
  ################## TM Results: SYBR #######################
  ## 7) Load data and generate output table
  ## a) Generate table for all genes
  TMTable <- function(){
    #Combine all input genes
    melt_gene <- TMinput()
    fluos_gene <- matchAllTarget()
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
  }
  
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
  
  
}
