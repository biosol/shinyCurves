library("shiny")
library("data.table")
library("qpcR")
library("reshape2")
library("cowplot")
library("rlist")
library("DT")
library("readxl")
library("plyr")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("matrixStats")
library("ggpmisc")
library("outviz")
library("tidyverse")
library("tidyr")
library("janitor")
library("shinyWidgets")
library("ggpubr")
library("stringi")
################################# SHINY APP (COMPLETE) FOR QPCR ANALYSIS ############################################
########################### Sonia Olaechea-Lázaro (UPV/EHU, May 2020) ############################################  

server <- function(input, output) {
  
  ################ TOY DATASET #################
  output$downtoy <- downloadHandler(
    filename <- function() {
      paste("toy", "zip", sep=".")
    },
    
    content <- function(file) {
      file.copy("toy.zip", file)
    },
    contentType = "toy.zip"
  )
  
  ############### MANUAL ################
  output$downmanual <- downloadHandler(
    filename <- function() {
      paste("Tutorial_App", "docx", sep=".")
    },
    
    content <- function(file) {
      file.copy("Tutorial_App.docx", file)
    },
    contentType = "Tutorial_App.docx"
  )
  
  
  ################## BIORAD - TAQMAN ##################
  #####################################################
  
  ################## Empty Plate Tab: Biorad #####################
  'platesTable <- function(cnames, controls){
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
  )'
  
  ################## Read "Run Information" and "Data" Sheet: Biorad ##########################
  # This function allows the user to upload either 1 Excel with 2 sheets (Data and Run Information)
  # or to upload 2 independent CSVs with the Run Information and Quantification Summary (with Sample,
  # Fluor, Target, Cq, etc columns)
  readBiorad <- function(){
    res <- input$biorad
    tbl_list <- list()
    ## Read data from CSVs
    if (!is.null(res)){
      if (grepl("csv",res$datapath[1]) == TRUE){
        tbl_list <- lapply(res$datapath, read.csv)
        tmp <- as.data.frame(tbl_list[1]) %>%
          select(Well, Fluor, Target, Sample, Cq)
        colnames(tmp) <- c("Well", "Fluor", "Target", "ID", "Cq")
        
        ## Correct A1 to A01, etc
        newwell <- vector()
        for (i in tmp$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        tmp$Well <- newwell
        colnames(tmp) <- c("Well", "Fluor", "Target", "ID", "Cq")
        tmp$Cq <- format(tmp$Cq, digits = 4)
        
        tbl_list[[1]] <- tmp
        
        ## Read Data from 1 Excel
      } else if (grepl("xlsx", res$datapath[1]) == TRUE){
        dat <- read_xlsx(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target",  "ID", "Cq")
        
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        
        df$Well <- newwell
        df$Cq <- format(df$Cq, digits = 4)
        
        
        ## Read Run Information
        run <- read_xlsx(res$datapath, sheet = "Run Information")
        df2 <- as.data.frame(run)
        ## Append to list
        tbl_list[[1]] <- df
        tbl_list[[2]] <- df2
        
      } else if (grepl("xls$", res$datapath[1]) == TRUE){
        dat <- read_xls(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target", "ID", "Cq")
        
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        
        df$Well <- newwell
        df$Cq <- format(df$Cq, digits = 4)
        
        
        ## Read Run Information
        run <- read_xls(res$datapath, sheet = "Run Information")
        df2 <- as.data.frame(run)
        ## Append to list
        tbl_list[[1]] <- df
        tbl_list[[2]] <- df2
      }
      
      ## Gene List
      inp <- tbl_list[[1]]
      genes <- unique(inp$Target)
      genes <- as.character(genes)
      genes[genes == ""] <- NA
      genes <- as.character(na.exclude(genes))
      
      # Add to list
      tbl_list[[3]] <- genes
      return(tbl_list)
    }
    
  }
  
  output$BRruninfo <- renderTable(
    if (!is.null(readBiorad()[[2]])){
      readBiorad()[[2]]
    } else{
      "No run information to show"
    }
    
  )
  
  output$inputdf <- renderTable(
    readBiorad()[[1]]
  )
  
  ################## Ct Plate Tab: Biorad ######################
  ctPlate <- function(){
    inp <- readBiorad()[[1]]
    r <- vector()
    c <- vector()
    newwell <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    ncols <- length(c)
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(inp$Well)){
      ro <- str_sub(inp$Well[i],1,1)
      col <- as.numeric(str_sub(inp$Well[i],-2))
      df[ro,col] <- inp$Cq[i]
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printCtPlate <- function(){
    dt <- ctPlate()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20)) %>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  #### Render Table in App ### 
  output$ctplate<- renderDataTable(
    printCtPlate()
  )
  
  ################## Sample Plate Tab: Biorad ############################
  samplePlate <- function(){
    inp <- readBiorad()[[1]]
    inp$ID <- as.character(inp$ID)
    r <- vector()
    c <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    ncols <- length(c)
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(inp$Well)){
      ro <- str_sub(inp$Well[i],1,1)
      col <- as.numeric(str_sub(inp$Well[i],-2))
      df[ro,col] <- as.character(inp$ID[i])
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny ##
  printSamplePlate <- function(){
    dt <- samplePlate()
    defdt <- datatable(dt, rownames = T,  options = list(pageLength = 20))
    return(defdt)
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
    if (input$posctrlbiorad != 0){
      inp <- readBiorad()[[1]]
      genes <- readBiorad()[[3]]
      # Get control rows
      ctrls <- lapply(genes, function(x){
        inp[grepl(paste(x,"_",sep=""), inp$ID),]
      })
      ctrls <- bind_rows(ctrls)
      
      ntcs <- lapply(genes, function(x){
        inp[grepl(paste("NTC_",x,sep=""), inp$ID),]
      })
      ntcs <- bind_rows(ntcs)
      
      ctrls <- as.data.frame(rbind(ctrls, ntcs))
      
      # Prepare data frame depending on nb of controls and duplicates
      nbctrls <- as.character(unique(ctrls$ID))
      
      ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
      if(input$dupsbiorad == TRUE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlbiorad), ncol = 3+(length(genes)*2)+length(genes)*3))
        
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(rep(paste(x,"(Dup",1:2,")",sep=""),each=1),paste(x,"(Avg)", sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies <- lapply(a$Dilution, function(x){as.numeric(input$concstdbiorad)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Gene duplicates and avg
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")] <- tmp$Cq[1]
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")] <- tmp$Cq[2]
              a[tmp_s[2],paste(unique(tmp$Target),"(Avg)",sep="")] <- mean(c(a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")],a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")]), na.rm = TRUE)
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocbiorad, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","Fluor","Target","ID","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")] <- negc$Cq[1]
        a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")] <- negc$Cq[2]
        a["C(-)", paste(unique(negc$Target),"(Avg)",sep="")] <- mean(c(a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")],a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")]), na.rm = TRUE)
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        ntc <- bind_rows(lapply(ntc, rbind))
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc),]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Fluor","Target","ID","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")] <- df$Cq[1]
          a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")] <- df$Cq[2]
          a["NTC", paste(unique(df$Target),"(Avg)", sep = "")] <- mean(a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")],a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")], na.rm = TRUE)
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        coefficients <- as.data.frame(coefficients)
        names(coefficients)[names(coefficients) == "cp"] <- paste(input$endocapplied,"(Avg)",sep="")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        return(a)
        
        ## If user has no duplicates
      } else if (input$dupsbiorad == FALSE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlbiorad), ncol = 3+length(genes)*3))
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(paste(x,"(Ct)",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdbiorad)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Genes control Cq
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Ct)",sep="")] <- tmp$Cq
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocbiorad, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","Fluor","Target","ID","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Ct)", sep="")] <- negc$Cq
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc)]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Fluor","Target","ID","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Ct)",sep = "")] <- df$Cq
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        
        return(a)
      }
      
    } else if (input$posctrlbiorad == 0){
      return(NULL)
    }
    
  }
  
  ##### Render Table
  output$stdcurve <- renderTable(
    if (!is.null(stCurve())){
      stCurve()
    } else {
      "No standard curve to show"
    },
    rownames = TRUE
  )
  
  ####### Coefficients function ########
  stdCoeffs <- function(){
    a <- stCurve()
    if (input$dupsbiorad == TRUE){
      cp <- a$logCopies
      avgs <- a[,grep("Avg", names(a))]
      avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      
      return(coefficients2)
    } else if (input$dupsbiorad == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Ct", names(a))]
      avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    }
    
  }
  
  ################## Standard Curve Tab (Plots): Biorad ##########################
  stdPlots <- function(){
    if (input$posctrlbiorad != 0){
      a <- stCurve()
      g <- readBiorad()[[3]]
      genes <- grep(input$endoC, g, value = TRUE,invert = TRUE)
      
      if(input$dupsbiorad == TRUE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
        
      } else if (input$dupsbiorad == FALSE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocbiorad, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
      }
    } else if (input$posctrlbiorad == 0){
      return(NULL)
    }
    
  }
  
  ##### Render Plots in App
  output$std <- renderPlot(
    if (!(is.null(stdPlots()))){
      stdPlots()
    } else {
      NULL
    }
    
  )
  
  ################## Analysis Samples Tab: Biorad ###################
  AnalysisSamples <- function(){
    
    ## Real Id, Sample, Replicate columns
    d <- readBiorad()[[1]]
    genes <- readBiorad()[[3]]
    studygenes <- grep(input$endoC, genes, value = TRUE, invert = TRUE)
    d$Cq <- as.numeric(d$Cq)
    d$ID <- as.character(d$ID)
    for (i in 1:length(d$Cq)){
      if (is.na(d[i, "Target"]) == TRUE){
        next
      } else if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } 
    }
    
    ## If user has duplicates, calculate mean Cq
    if(input$dupsbiorad == TRUE & input$copiesforassig == TRUE){
      
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div >= 1.5){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctbiorad & as.numeric(paste(x)) <= input$cyclesbiorad){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endoC,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocbiorad){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocbiorad){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endoC,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks 
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgenbiorad & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgenbiorad & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## Log Copies
      coeffs <- stdCoeffs()
      
      studygenes <- genes[grep(input$endoC, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvbiorad){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvbiorad & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgenbiorad & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgenbiorad){
            "Positive"
          }else if (pos < input$numposgenbiorad & neg > pos & neg > rep){
            "Negative"
          }else{
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupsbiorad == FALSE & input$copiesforassig == TRUE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Ct")
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, value = TRUE ,invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctbiorad & as.numeric(paste(x)) <= input$cyclesbiorad){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$cyclesbiorad){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      
      controlcq <- sapply(mean_cq_merged[,paste(input$endoC,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocbiorad){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocbiorad){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endoC,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgenbiorad & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgenbiorad & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      coeffs <- stdCoeffs()
      studygenes <- genes[grep(input$endoC, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (mean_cq_merged[j,paste(studygenes[i],"(Ct)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(Ct)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Ct)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Ct)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvbiorad){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvbiorad & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgenbiorad & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgenbiorad){
            "Positive"
          }else if (pos < input$numposgenbiorad & neg > pos & neg > rep){
            "Negative"
          } else{
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupsbiorad == TRUE & input$copiesforassig == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div >= 1.5){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctbiorad){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endoC,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocbiorad){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocbiorad){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endoC,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      mean_cq_merged[["Assignment"]] <- ctassig
      
      return(mean_cq_merged)
    } else if (input$dupsbiorad == FALSE & input$copiesforassig == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, value = TRUE ,invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctbiorad){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endoC,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocbiorad){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocbiorad){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endoC,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenbiorad){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
    
  }
  
  output$downanbiorad <- downloadHandler(
    filename = "Analysis.csv",
    content = function(fname){
      write.csv(AnalysisSamples(), fname, quote = F, row.names = F)}
  )
  
  AnalysisSamplesDT <- function(){
    df <- AnalysisSamples()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    t <- defdt %>%
      formatStyle(
        columns = grep("Final", names(df)),
        backgroundColor = "khaki"
      ) %>%
      formatStyle(
        columns = "Assignment",
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat"), 
                                     c("limegreen", "tomato", "lightblue")),
        fontWeight = "bold"
      ) 
    
    return(t)
  }
  
  output$analysis <- renderDataTable({
    AnalysisSamplesDT()
  })
  
  ################## ID Well Tab: Biorad ##############################
  'SamplesToKeep <- function(){
    inp <- readBiorad()[[1]]
    g <- inp %>%
      group_by(ID)
    
    g_sp <- group_split(g)
    tokeep <- vector()
    for (g in g_sp){
      if(all(is.na(as.numeric(g$Cq))) == TRUE){
        next
      } else if (all(is.na(g$Cq)) == FALSE){
        tokeep <- append(tokeep, as.numeric(as.character(unique(g$ID))))
      } 
    }
    tokeep <- na.exclude(tokeep)
    return(tokeep)
  }'
  
  IDWELLtab <- function(){
    inp <- readBiorad()[[1]]
    wellid <- cbind(as.character(inp$Well), as.character(inp$Target), as.numeric(as.character(inp$ID)))
    wellid <- as.data.frame(wellid)
    colnames(wellid) <- c("Well", "Target", "ID")
    
    ## Get list of samples to keep
    #tokeep <- SamplesToKeep()
    ## Keep only samples from wellid df
    #wellid <- wellid[wellid$ID %in% tokeep,]
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
  
  
  ################## ID Result Tab: Biorad    ###############################
  idresult <- function(){
    an <- AnalysisSamples()
    
    id <-unique(as.numeric(as.character(an$ID)))
    idres <- as.data.frame(cbind(id, an$Assignment))
    colnames(idres)<- c("ID", "Interpretation")
    idres <- idres[which(is.na(idres$Interpretation) == FALSE),]
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ################## APPLIED QUANT STUDIO - TAQMAN ##################
  ###################################################################
  
  ################## Read Applied "Raw Data" Sheet: Applied ###################
  readApplied <- function(){
    res <- input$appl
    tbl_list <- list()
    if (!is.null(res)){
      if (grepl("xlsx",res$datapath[1]) == TRUE){
        ## Data
        samples <- read_xlsx(res$datapath, sheet = "Raw Data")
        df <- as.data.frame(samples)
        start <- as.numeric(as.character(grep("Well Position", df[[2]])))
        df1 <- df[-c(1:start),1:4]
        colnames(df1) <- c("Well", "Well_Position", "Cycle","Fluorescence")
        
        ## A1 to A01
        newwell <- vector()
        for(i in df1$Well_Position){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        
        df1$Well_Position <- newwell
        
        # Run Information
        df2 <- samples[1:(start-1),1:2]
        df2 <- as.data.frame(df2)
        
        # List of genes
        genes <- read_xlsx(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df3 <- df[-c(1:start),5]
        df3 <- unique(df3)
        df3 <- as.character(na.exclude(df3))
        
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Run Information
        tbl_list[[3]] <- df3 # List of genes
        return(tbl_list)
        
      } else if (grepl("xls$",res$datapath[1]) == TRUE){
        ## Data
        samples <- read_xls(res$datapath, sheet = "Raw Data")
        df <- as.data.frame(samples)
        start <- as.numeric(as.character(grep("Well Position", df[[2]])))
        df1 <- df[-c(1:start),1:4]
        colnames(df1) <- c("Well", "Well_Position", "Cycle","Fluorescence")
        
        ## A1 to A01
        newwell <- vector()
        for(i in df1$Well_Position){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df1$Well_Position <- newwell
        
        # Run Information
        df2 <- samples[1:(start-1),1:2]
        df2 <- as.data.frame(df2)
        
        # List of genes
        genes <- read_xls(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df3 <- df[-c(1:start),5]
        df3 <- unique(df3)
        df3 <- as.character(na.exclude(df3))
        
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Run Information
        tbl_list[[3]] <- df3 # List of genes
        
        return(tbl_list)
      }
    }
  }
  
  output$readapp <- renderTable(
    readApplied()[[1]]
  )
  
  ################## Read Results: Applied #####################
  readAppliedResults <- function(){
    res <- input$appl
    if (!is.null(res)){
      if (grepl("xlsx", res$datapath[1]) == TRUE){
        samples <- read_xlsx(res$datapath, sheet = "Results", na = "Undetermined")
        df <- data.frame(samples)
        start <- grep("Well Position", df[[2]])
        df <- df[-c(1:start),c(2,4,5,7,9)]
        colnames(df) <- c("Well","ID", "Target","Fluor","Cq")
        df$Cq <- format(df$Cq, digits = 4)
        
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df$Well <- newwell
        return(df)
        
      } else if (grepl("xls$", res$datapath[1]) == TRUE){
        samples <- read_xls(res$datapath, sheet = "Results", na = "Undetermined")
        df <- data.frame(samples)
        start <- grep("Well Position", df[[2]])
        df <- df[-c(1:start),c(2,4,5,7,9)]
        colnames(df) <- c("Well","ID", "Target","Fluor","Cq")
        df$Cq <- format(df$Cq, digits = 4)
        
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df$Well <- newwell
        return(df)
      }
    }
  }
  
  output$appliedres <- renderTable({
    readAppliedResults()
  })
  ################## Conversion Tab: Applied #####################
  conversionApplied <- function(){
    df <- readApplied()[[1]]
    res <- readAppliedResults()
    colnames(res)[1] <- "Well_Position"
    m <- merge(df, res, by = "Well_Position")
    m <- m[,c("Well_Position", "Cycle", "Fluorescence", "Target")]
    m_s <- m[order(m$Cycle, m$Target), ]
    
    final <- m_s %>%
      group_by(Target) %>%
      pivot_wider(names_from = Cycle, values_from = Fluorescence)
    final <- as.data.frame(final)
    rownames(final) <- final$Well_Position
    final$Well_Position <- NULL
    final_t<-t(final)
    final_t_s <- final_t[order(as.numeric(rownames(final_t))),]
    
    genes <- readApplied()[[3]]
    
    l <- lapply(genes, function(x){
      df <- final_t_s[, final_t_s["Target", ] == x]
      df <- as.data.frame(df)
      df <- df[!row.names(df) %in% "Target",]
      df <- df[!is.na(names(df))]
      rownames_to_column(df, "Cycle")
    })
    
    return(l)
  }
  
  output$conversionapp <- renderUI({
    genes <- readApplied()[[3]]
    nbgenes <- c(1:length(genes))
    tabs <- lapply(genes, function(x){
      tabPanel(
        title = uiOutput(x),
        uiOutput(paste("down",x,sep = "")),
        tableOutput(paste("trans",x,sep=""))
      )
    })
    
    ## Tab Names
    lapply(genes, function(x){
      output[[x]] <- renderText({
        x 
      })
    })
    
    ## Download Plots (Button)
    lapply(genes, function(x){
      output[[paste("down",x,sep="")]] <- renderUI({
        ls <- conversionApplied()
        if (is.list(ls) == TRUE){
          downloadButton(paste("downl",x,sep=""), "Download CSV")
        } else{
          NULL
        }
      })
    })
    
    ## Download Plots (Handler)
    lapply(nbgenes, function(x){
      output[[paste("downl",genes[x],sep="")]] <- downloadHandler(
        filename = paste(genes[x],".csv",sep=""),
        content = function(fname){
          write.csv(conversionApplied()[x], fname, quote = F, row.names = F, sep = ",", dec = ".")
        }
      )
    })
    
    ## Render Tables
    lapply(nbgenes, function(x){
      output[[paste("trans",genes[x],sep="")]] <- renderTable({
        conversionApplied()[x]
      })
    })
    
    do.call(tabsetPanel,c(tabs))
  })
  
  ################## Read Run Info: Applied ####################
  output$appliedruninfo <- renderTable({
    if (!is.null(readApplied()[[2]])){
      readApplied()[[2]]
    } else {
      "No run information to show"
    }
  })
  
  ################## Ct Plate Tab: Applied ##################
  ctPlateApp <- function(){
    inp <- readAppliedResults()
    r <- vector()
    c <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    ncols <- length(c)
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(inp$Well)){
      ro <- str_sub(inp$Well[i],1,1)
      col <- as.numeric(str_sub(inp$Well[i],-2))
      df[ro,col] <- inp$Cq[i]
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printCtPlateApp <- function(){
    dt <- ctPlateApp()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20)) %>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  output$ctplateapp<- renderDataTable(
    printCtPlateApp()
  )
  
  ################## Sample Plate Tab: Applied ##########################
  samplePlateApp <- function(){
    inp <- readAppliedResults()
    inp$ID <- as.character(inp$ID)
    r <- vector()
    c <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    ncols <- length(c)
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(inp$Well)){
      ro <- str_sub(inp$Well[i],1,1)
      col <- as.numeric(str_sub(inp$Well[i],-2))
      df[ro,col] <- as.character(inp$ID[i])
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printSamplePlateApp <- function(){
    dt <- samplePlateApp()
    defdt <- datatable(dt, rownames = T,  options = list(pageLength = 20))
    return(defdt)
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
    if (input$posctrlapplied != 0){
      inp <- readAppliedResults()
      genes <- readApplied()[[3]]
      
      # Get control rows
      ctrls <- lapply(genes, function(x){
        inp[grepl(paste(x,"_",sep=""), inp$ID),]
      })
      ctrls <- bind_rows(ctrls)
      
      ntcs <- lapply(genes, function(x){
        inp[grepl(paste("NTC_",x,sep=""), inp$ID),]
      })
      ntcs <- bind_rows(ntcs)
      
      ctrls <- as.data.frame(rbind(ctrls, ntcs))
      
      # Prepare data frame depending on nb of controls and duplicates
      nbctrls <- as.character(unique(ctrls$ID))
      
      ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
      if(input$dupsapplied == TRUE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlapplied), ncol = 3+(length(genes)*2)+length(genes)*3))
        
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(rep(paste(x,"(Dup",1:2,")",sep=""),each=1),paste(x,"(Avg)", sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdapplied)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        
        ## Gene duplicates and avg
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")] <- tmp$Cq[1]
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")] <- tmp$Cq[2]
              a[tmp_s[2],paste(unique(tmp$Target),"(Avg)",sep="")] <- mean(c(a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")],a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")]), na.rm = TRUE)
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocapplied, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","ID","Target","Fluor","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")] <- negc$Cq[1]
        a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")] <- negc$Cq[2]
        a["C(-)", paste(unique(negc$Target),"(Avg)",sep="")] <- mean(c(a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")],a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")]), na.rm = TRUE)
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        ntc <- bind_rows(lapply(ntc, rbind))
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc),]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","ID","Target","Fluor","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")] <- df$Cq[1]
          a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")] <- df$Cq[2]
          a["NTC", paste(unique(df$Target),"(Avg)", sep = "")] <- mean(a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")],a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")], na.rm = TRUE)
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        coefficients <- as.data.frame(coefficients)
        names(coefficients)[names(coefficients) == "cp"] <- paste(input$endocapplied,"(Avg)",sep="")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        
        return(a)
        
        ## If user has no duplicates
      } else if (input$dupsapplied == FALSE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlapplied), ncol = 3+length(genes)*3))
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(paste(x,"(Ct)",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdapplied)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Genes control Cq
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Ct)",sep="")] <- tmp$Cq
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocapplied, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","ID","Target","Fluor","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Ct)", sep="")] <- negc$Cq
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc)]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","ID","Target","Fluor","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Ct)",sep = "")] <- df$Cq
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        
        return(a)
      }
    } else if (input$posctrlapplied){
      return(NULL)
    }
    
    
  }
  
  ##### Render Table
  output$stdcurveapp <- renderTable(
    if (!is.null(stCurveApp())){
      stCurveApp()
    } else {
      "No standard curve to show" 
    }, rownames = TRUE
  )
  
  ####### Coefficients function #######
  stdCoeffsApp <- function(){
    a <- stCurveApp()
    if (input$dupsapplied == TRUE){
      cp <- a$logCopies
      avgs <- a[,grep("Avg", names(a))]
      avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    } else if (input$dupsapplied == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Ct", names(a))]
      avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    }
    
  }
  
  ################## Standard Curve Tab (Plots): Applied ###########################
  stdPlotsApp <- function(){
    if (input$posctrlapplied != 0){
      a <- stCurveApp()
      g <- readAppliedResults()
      genes <- readApplied()[[3]]
      genes <- grep(input$endocapplied, genes, value = TRUE, invert = TRUE)
      
      if(input$dupsapplied == TRUE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
        
      } else if (input$dupsapplied == FALSE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocapplied, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
      }
    } else if (input$posctrlapplied == 0){
      return(NULL)
    }
  }
  
  ##### Render Plots in App
  output$stdapp <- renderPlot(
    stdPlotsApp()
  )
  
  ################## Analysis Samples Tab: Applied ######################
  AnalysisSamplesApp <- function(){
    
    ## Real Id, Sample, Replicate columns
    d <- readAppliedResults()
    genes <- readApplied()[[3]]
    studygenes <- grep(input$endocapplied, genes, value = TRUE, invert = TRUE)
    d$Cq <- as.numeric(d$Cq)
    d$ID <- as.character(d$ID)
    for (i in 1:length(d$Cq)){
      if (is.na(d[i, "Target"]) == TRUE){
        next
      } else if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } 
    }
    
    ## If user has duplicates, calculate mean Cq
    if(input$dupsapplied == TRUE & input$copiesforassigapplied == TRUE){
      
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div >= 1.5){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctapplied & as.numeric(paste(x)) <= input$cyclesapplied){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$cyclesapplied){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocapplied){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocapplied){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgenapplied & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgenapplied & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## LogCopies
      coeffs <- stdCoeffsApp()
      studygenes <- genes[grep(input$endocapplied, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvapplied){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvapplied & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgenapplied & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgenapplied){
            "Positive"
          }else if (pos < input$numposgenapplied & neg > pos & neg > rep){
            "Negative"
          } else{
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupsapplied == FALSE & input$copiesforassigapplied == TRUE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, value = TRUE ,invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctapplied & as.numeric(paste(x)) <= input$cyclesapplied){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$cyclesapplied){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocapplied){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocapplied){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgenapplied & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgenapplied & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      coeffs <- stdCoeffsApp()
      studygenes <- genes[grep(input$endocapplied, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvapplied){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvapplied & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgenapplied & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgenapplied){
            "Positive"
          }else if (pos < input$numposgenapplied & neg > pos & neg > rep){
            "Negative"
          } else {
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupsapplied == TRUE & input$copiesforassigapplied == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div >= 1.5){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) >= input$minctapplied){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocapplied){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocapplied){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      
      return(mean_cq_merged)
    } else if (input$dupsapplied == FALSE & input$copiesforassigapplied == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, value = TRUE ,invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctapplied){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocapplied){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocapplied){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocapplied,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgenapplied){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
    
  }
  
  output$downanapplied <- downloadHandler(
    filename = "Analysis.csv",
    content = function(fname){
      write.csv(AnalysisSamplesApp(), fname, quote = F, row.names = F)}
  )
  
  AnalysisSamplesDTApp <- function(){
    df <- AnalysisSamplesApp()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    t <- defdt %>%
      formatStyle(
        columns = grep("Final", names(df)),
        backgroundColor = "khaki"
      ) %>%
      formatStyle(
        columns = "Assignment",
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat"), 
                                     c("limegreen", "tomato", "lightblue")),
        fontWeight = "bold"
      ) 
    
    return(t)
  }
  
  output$analysisapp <- renderDataTable({
    AnalysisSamplesDTApp()
  })
  
  ################## ID Well Tab: Applied ###############
  IDWELLtabApp <- function(){
    inp <- readAppliedResults()
    wellid <- cbind(inp$Well, as.character(inp$Target), inp$ID)
    wellid <- as.data.frame(wellid)
    colnames(wellid) <- c("Well", "Target", "ID")
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
    realid <- unique(as.numeric(as.character(an$ID)))
    final <- as.data.frame(cbind(realid, an$Assignment))
    colnames(final)<- c("ID", "Interpretation")
    return(final)
  }
  
  output$IDRESULTApp <- renderTable(
    idresultApp()
  )
  
  output$downidresapp <- downloadHandler(
    filename = "ID_result.csv",
    content = function(fname){
      write.csv(idresultApp(), fname, quote = F, row.names = F)}
  )  
  
  
  
  
  ################## CT CURVES - TAQMAN ##################
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
      df<-info[which(info$Target==genes[i]),] # select target
      df<-df[order(df$Well2),]
      
      ## load raw data (fluorescence-RFU per temperature) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ## Read fluorescence data for target gene
      inp <- cyclesInput()
      
      raw <- inp[i][[1]]
      colnames(raw) <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', colnames(raw), perl=TRUE))
      raw<-raw[,df$Well2]
      raw_s <- stack(raw)
      
      raw_s$cycles<-rep(1:input$ct,lr[[2]]) # rep(1:cycle number,sample number x2 replicates)
      colnames(raw_s)<-c("RFU","Well2","Cycles")
      merge<-merge(raw_s,df,by="Well2")
      
      ###plot
      pltName <- paste(genes[i])
      merge$cols <- ifelse(merge$Interpretation == "Negative", 'tomato', ifelse(merge$Interpretation == "Positive", "green3", "blue"))
      pltList[[ pltName ]] = ggplot(data = merge)+
        geom_line(aes(x=Cycles, y=RFU, by=Well2, color=Interpretation))+
        scale_color_manual(values = c("tomato", "green3", "blue"))+
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
      colnames(raw) <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', colnames(raw), perl=TRUE))
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
      
      if (length(merge_r_persample) == 0){
        mes <- "No indetermined samples to plot"
        return(mes)
      } else {
        ## Plots
        pltList <- list()
        for (z in 1:length(merge_r_persample)){
          pltName <- paste(genes[i],"_",names(merge_r_persample)[[z]], sep = "")
          merge_pn$cols <- ifelse(merge_pn$Interpretation == "Negative", 'tomato', ifelse(merge_pn$Interpretation == "Positive", "green3", "blue"))
          pltList[[pltName]] <- ggplot(data = merge_pn, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
            geom_line()+
            geom_line(data = merge_r_persample[[z]], aes(x=Cycles, y=RFU, by=Well2, color=Interpretation))+
            scale_color_manual(values = c("tomato", "green3", "blue"))+
            ylim(-100, lr[[1]])+ 
            theme(legend.position = "none")+ 
            ggtitle(paste(genes[i],"_",names(merge_r_persample)[[z]]))
        }
        defPltLs[[i]] <- pltList
      }
    }
    return(defPltLs) 
  }
  
  ### Generate # tabs with X names depeding on user input files
  output$indetplots <- renderUI({
    genes <- geneList()
    nbgenes <- c(1:length(geneList()))
    tabs <- lapply(genes, function(x){
      tabPanel(
        title=uiOutput(x), 
        textOutput(paste("text",x,sep="")), 
        uiOutput(paste("down",x,sep="")), 
        plotOutput(paste("indet",x,sep=""), height = "500px")
      )
    })
    
    ## Tab Names
    lapply(genes, function(x){
      output[[x]] <- renderText({
        x 
      })
    })
    
    ## Message for no "Repeat" samples
    lapply(genes, function(x){
      output[[paste("text",x,sep="")]] <- renderText({
        ls <- indetPlots()
        if (is.list(ls) == FALSE){
          print(ls)
        } else {
          print(NULL)
        } 
      })
    })
    
    ## Download Plots (Button)
    lapply(genes, function(x){
      output[[paste("down",x,sep="")]] <- renderUI({
        ls <- indetPlots()
        if (is.list(ls) == TRUE){
          downloadButton(paste("downl",x,sep=""), "Download PDF")
        } else{
          NULL
        }
      })
    })
    
    ## Download Plots (Handler)
    lapply(nbgenes, function(x){
      output[[paste("downl",genes[x],sep="")]] <- downloadHandler(
        filename = paste(genes[x],"_IndetPlots.pdf",sep=""),
        content = function(file){
          plots <- indetPlots()
          p <- plot_grid(plotlist = plots[[x]], ncol = 2)
          save_plot(file, p, ncol = 1, base_height = 6, base_width = 8)
        }
      )
    })
    
    ## Render Plots
    lapply(nbgenes, function(x){
      output[[paste("indet",genes[x],sep="")]] <- renderPlot({
        plots <- indetPlots()
        if (is.list(plots) == TRUE){
          plot_grid(plotlist = plots[[x]], ncol = 2)
        } else {
        }
      })
    })
    
    do.call(tabsetPanel,c(tabs))
    
  })
  
  
  
  
  
  
  
  
  ################### MELTING CURVE - SYBR #########################
  ###############################################
  
  ################## Read Applied "Raw Data" Sheet: SYBR-Applied ###################
  'readAppliedSYBRRFU <- function(){
    res <- input$sybrapp
    tbl_list <- list()
    if (!is.null(res)){
      if (grepl("xls$",res$datapath[1]) == TRUE){
        ## Data
        samples <- read_xls(res$datapath, sheet = "Raw Data")
        df <- as.data.frame(samples)
        start <- as.numeric(as.character(grep("Well Position", df[[2]])))
        df1 <- df[-c(1:start),1:4]
        colnames(df1) <- c("Well", "Well_Position", "Cycle","Fluorescence")
        
        
        newwell <- vector()
        for(i in df1$Well_Position){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df1$Well_Position <- newwell
        
        # Samples with melting OK
        idwell <- IDWELLtabSYBRApp()
        idwell$Well <- as.character(sub("(?<![0-9])0*(?=[0-9])", "", idwell$Well, perl=TRUE))
        idwell$ID <- as.character(idwell$ID)
        ## Add controls, dils
        ctrls <- tmp[grep("negC|10-|NTC", df1$ID),]
        m <- tmp[df1$Well %in% idwell$Well,]
        df2 <- as.data.frame(rbind(m, ctrls))
        
        # Run Information
        df3 <- samples[1:(start-1),1:2]
        df3 <- as.data.frame(df3)
        
        # List of genes
        genes <- read_xlsx(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df4 <- df[-c(1:start),5]
        df4 <- unique(df4)
        
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Samples with melting OK
        tbl_list[[3]] <- df3 # Run Information
        tbl_list[[4]] <- df4 # List of genes
        
        return(tbl_list)
      }
    }
  }
  
  output$readappsybr <- renderTable(
    readAppliedSYBR()[[2]]
  )'
  
  
  
  ################## Conversion Tab: SYBR-Applied #####################
  'conversionAppliedSYBR <- function(){
    df <- readAppliedSYBR()[[2]]
    res <- readAppliedResultsSYBR()
    colnames(res)[1] <- "Well_Position"
    
    m <- merge(df, res, by = "Well_Position")
    m <- m[,c("Well_Position", "Cycle", "Fluorescence", "Target")]
    m_s <- m[order(m$Cycle, m$Target), ]
    
    
    final <- m_s %>%
      group_by(Target) %>%
      pivot_wider(names_from = Cycle, values_from = Fluorescence)
    
    final <- as.data.frame(final)
    rownames(final) <- final$Well_Position
    final$Well_Position <- NULL
    
    final_t<-t(final)
    final_t_s <- final_t[order(as.numeric(rownames(final_t))),]
    
    genes <- readAppliedSYBR()[[4]]
    
    l <- lapply(genes, function(x){
      df <- final_t_s[, final_t_s["Target", ] == x]
      df <- as.data.frame(df)
      df <- df[!row.names(df) %in% "Target",]
      rownames_to_column(df, "Cycle")
    })
    
    return(l)
  }
  
  output$conversionsybrapp <- renderUI({
    genes <- readAppliedSYBR()[[4]]
    nbgenes <- c(1:length(genes))
    tabs <- lapply(genes, function(x){
      tabPanel(
        title = uiOutput(x),
        uiOutput(paste("down",x,sep = "")),
        tableOutput(paste("trans",x,sep=""))
      )
    })
    
    ## Tab Names
    lapply(genes, function(x){
      output[[x]] <- renderText({
        x 
      })
    })
    
    ## Download Plots (Button)
    lapply(genes, function(x){
      output[[paste("down",x,sep="")]] <- renderUI({
        ls <- conversionAppliedSYBR()
        if (is.list(ls) == TRUE){
          downloadButton(paste("downl",x,sep=""), "Download CSV")
        } else{
          NULL
        }
      })
    })
    
    ## Download Plots (Handler)
    lapply(nbgenes, function(x){
      output[[paste("downl",genes[x],sep="")]] <- downloadHandler(
        filename = paste(genes[x],".csv",sep=""),
        content = function(fname){
          write.csv(conversionAppliedSYBR()[x], fname, quote = F, row.names = F)
        }
      )
    })
    
    ## Render Tables
    lapply(nbgenes, function(x){
      output[[paste("trans",genes[x],sep="")]] <- renderTable({
        conversionAppliedSYBR()[x]
      })
    })
    
    do.call(tabsetPanel,c(tabs))
  })'
  
  
  
  ################## Read Excel to generate ID_well #############
  readSYBRIDWell <- function(){
    res <- input$sybrinp
    if (!is.null(res)){
      if (grepl("csv",res$datapath[[1]][1]) == TRUE){
        tbl_list <- lapply(res$datapath, read.csv)
        tmp <- as.data.frame(tbl_list[1]) %>%
          select(Well, Target, Sample)
        colnames(tmp) <- c("Well", "Target", "ID")
        
        newwell <- vector()
        for(i in tmp$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        tmp$Well <- newwell
        colnames(tmp) <- c("Well", "Target", "ID")
        
        return(tmp)
        
      } else if (grepl("xlsx", res$datapath[[1]][1]) == TRUE){
        ## Read Data from Biorad
        if ("^Data$" %in% excel_sheets(res$datapath)){
          dat <- read_xlsx(res$datapath, sheet = "Data")
          df <- as.data.frame(dat)
          df <- cbind(df$Well, df$Target, df$Sample)
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Target","ID")
          
          ## A01 to A1
          newwell <- vector()
          for(i in df$Well){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df$Well <- newwell
          colnames(df) <- c("Well", "Target","ID")
          
          ## Append to list
          return(df)
          
          # Read "Results" from Applied 
        } else if ("Results" %in% excel_sheets(res$datapath)){
          dat <- read_xlsx(res$datapath, sheet = "Results")
          df <- as.data.frame(dat)
          start <- as.numeric(as.character(grep("Well Position", df[[2]])))
          names(df) <- df[start,]
          dfi <- df[-c(1:start),]
          df1 <- dfi %>%
            select("Well Position", "Target Name", "Sample Name")
          colnames(df1) <- c("Well","Target","ID")
          ## A01 to A1
          newwell <- vector()
          for(i in df1$Well){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df1$Well <- newwell
          
          ## Append to list
          return(df1)
        }
        
      } else if (grepl("xls$", res$datapath[[1]][1]) == TRUE){
        if ("^Data$" %in% excel_sheets(res$datapath)){
          dat <- read_xls(res$datapath, sheet = "Data")
          df <- as.data.frame(dat)
          df <- cbind(df$Well, df$Target, df$Sample)
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Target","ID")
          
          ## A01 to A1
          newwell <- vector()
          for(i in df$Well){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df$Well <- newwell
          colnames(df) <- c("Well", "Target","ID")
          
          ## Append to list
          return(df)
        } else if ("Results" %in% excel_sheets(res$datapath)){
          dat <- read_xls(res$datapath, sheet = "Results")
          df <- as.data.frame(dat)
          start <- as.numeric(as.character(grep("Well Position", df[[2]])))
          names(df) <- df[start,]
          dfi <- df[-c(1:start),]
          df1 <- dfi %>%
            select("Well Position", "Target Name", "Sample Name")
          colnames(df1) <- c("Well","Target","ID")
          ## A01 to A1
          newwell <- vector()
          for(i in df1$Well){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df1$Well <- newwell
          
          ## Append to list
          return(df1)
        }
      }
    }
  }
  
  output$rawidwell <- renderTable({
    if(is.null(readSYBRIDWell())){
      return();
    }else{
      readSYBRIDWell();
    }
  })
  
  ################## Read CSVs with Tm: SYBR ######### 
  TMinput <- function(){
    dat <- input$sybrcsv
    if (!is.null(dat)){
      # Read independent CSVs from Taqman, else xlsx from SYBR
      if (length(dat$datapath) != 1){
        tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=",", dec=".")
        tbl_list <- lapply(tbl_list, function(x){x[1]<-NULL;x})
        return(tbl_list)
        
      } else {
        # Read data from Applied   
        if (grepl("xlsx",dat$datapath[1]) == TRUE){
          ## Data
          samples <- read_xlsx(dat$datapath, sheet = "Melt Curve Raw Data")
          df <- as.data.frame(samples)
          start <- as.numeric(as.character(grep("Well Position", df[[2]])))
          names(df) <- df[start,]
          dfi <- df[-c(1:start),]
          df1 <- dfi %>%
            select("Well", "Well Position","Temperature", "Fluorescence")
          colnames(df1) <- c("Well", "Well_Position", "Temperature","Fluorescence")
          
          newwell <- vector()
          for(i in as.character(df1$Well_Position)){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df1$Well_Position <- newwell
          
          idwell <- readSYBRIDWell()
          colnames(idwell) <- c("Well_Position", "Target", "ID")
          m <- merge(df1, idwell, by = "Well_Position")
          m <- m[,c("Well_Position", "Temperature", "Fluorescence", "Target")]
          m_s <- m[order(m$Temperature, m$Target), ]
          
          final <- m_s %>%
            group_by(Target) %>%
            pivot_wider(names_from = Temperature, values_from = Fluorescence)
          
          final <- as.data.frame(final)
          rownames(final) <- final$Well_Position
          final$Well_Position <- NULL
          
          final_t<-t(final)
          final_t_s <- final_t[order(as.numeric(rownames(final_t))),]
          
          genes <- unique(idwell$Target)
          
          l <- lapply(genes, function(x){
            df <- final_t_s[, final_t_s["Target", ] == x]
            df <- as.data.frame(df)
            df <- df[!row.names(df) %in% "Target",]
            rownames_to_column(df, "Temperature")
          })
          
          return(l)
          
        } else if (grepl("xls$",dat$datapath[1]) == TRUE){
          samples <- read_xls(dat$datapath, sheet = "Melt Curve Raw Data")
          df <- as.data.frame(samples)
          start <- as.numeric(as.character(grep("Well Position", df[[2]])))
          names(df) <- df[start,]
          dfi <- df[-c(1:start),]
          df1 <- dfi %>%
            select("Well", "Well Position","Temperature", "Fluorescence")
          colnames(df1) <- c("Well", "Well_Position", "Temperature","Fluorescence")
          
          newwell <- vector()
          for(i in as.character(df1$Well_Position)){
            if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
              stri_sub(i, 2, 1) <- 0
              newwell <- append(newwell, i)
            } else {
              newwell <- append(newwell, i)
            }
          }
          df1$Well_Position <- newwell
          
          idwell <- readSYBRIDWell()
          colnames(idwell) <- c("Well_Position", "Target", "ID")
          m <- merge(df1, idwell, by = "Well_Position")
          m <- m[,c("Well_Position", "Temperature", "Fluorescence", "Target")]
          m_s <- m[order(m$Temperature, m$Target), ]
          
          final <- m_s %>%
            group_by(Target) %>%
            pivot_wider(names_from = Temperature, values_from = Fluorescence)
          
          final <- as.data.frame(final)
          rownames(final) <- final$Well_Position
          final$Well_Position <- NULL
          
          final_t<-t(final)
          final_t_s <- final_t[order(as.numeric(rownames(final_t))),]
          
          genes <- unique(idwell$Target)
          
          l <- lapply(genes, function(x){
            df <- final_t_s[, final_t_s["Target", ] == x]
            df <- as.data.frame(df)
            df <- df[!row.names(df) %in% "Target",]
            rownames_to_column(df, "Temperature")
          })
          print(str(l))
          return(l)
        }
      }
    }
  }
  
  output$rawdata <- renderTable({
    TMinput()
  })
  
  ################## Gene List: SYBR ###########
  SybrGeneList <- function(){
    dat <- input$sybrcsv
    if (length(dat$datapath) != 1){
      lst <- lapply(dat$name, FUN = function(x) gsub(pattern = ".*[_]([^.]+)[.].*", replacement = "\\1",
                                                     basename(dat$name)))
      lst <- unlist(lst)
      out <- unique(lst)
      return(out)
    } else {
      if (grepl("xlsx",dat$datapath[1]) == TRUE){
        dat <- read_xlsx(dat$datapath, sheet = "Results")
        df <- as.data.frame(dat)
        start <- as.numeric(as.character(grep("Well Position", df[[2]])))
        names(df) <- df[start,]
        dfi <- df[-c(1:start),]
        df1 <- dfi %>%
          select("Target Name")
        colnames(df1) <- "Target"
        genes <- as.vector(unique(df1$Target))
        return(genes)
      } else if (grepl("xls$",dat$datapath[1]) == TRUE){
        dat <- read_xls(dat$datapath, sheet = "Results")
        df <- as.data.frame(dat)
        start <- as.numeric(as.character(grep("Well Position", df[[2]])))
        names(df) <- df[start,]
        dfi <- df[-c(1:start),]
        df1 <- dfi %>%
          select("Target Name")
        colnames(df1) <- "Target"
        genes <- as.vector(unique(df1$Target))
        return(genes)
      }
    }
    
  }
  
  ################## Fluos_Gene Tab: SYBR ###############
  ##### Select target gene and match columns in fluorescence file
  ## matchTarget por archivo independiente
  matchTarget <- function(tm, gene){
    # Modify colnames of tm
    Well <- readSYBRIDWell()
    # Remove ctrls and dils, we don't want to plot them
    Well <- Well[grep("negC|10-|NTC", Well$ID, invert = TRUE),]
    
    Well$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', Well$Well, perl=TRUE))
    well_gene <- Well[which(Well$Target==gene),] # select target
    
    colnames(well_gene) <- c("Well", "Target", "ID")
    
    #Select fluorescence matching columns
    fluos_gene <- tm[,well_gene$Well]
    
    return(fluos_gene)
  }
  
  ## Misma función que matchTarget pero con un loop para guardar todos los archivos en la misma tabla
  ## Bastante guarro, pero funciona
  matchAllTarget <- function(){
    Well <- readSYBRIDWell()
    # Remove ctrls and dils, we don't want to plot them
    Well <- Well[grep("negC|10-|NTC", Well$ID, invert = TRUE),]
    
    genes <- SybrGeneList()
    melt_gene <- TMinput()
    
    matchLs <- list()
    for (i in 1:length(genes)){
      #Select matching target rows from well_id
      Well$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', Well$Well, perl=TRUE))
      well_gene <- Well[which(Well$Target==genes[i]),] # select target
      #Select fluorescence matching columns
      fluos_gene<-melt_gene[[i]][,well_gene$Well]
      matchLs[[i]] <- fluos_gene
    }
    return(matchLs)
  }
  
  output$fluosgene <- renderTable({
    matchAllTarget()
  })
  
  ################## Melting Curve Plots: SYBR ##################
  # Function to Repeat temperature columns
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
    data_gene <- as.data.frame(data_gene)
    if (input$isderiv == TRUE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=as.numeric(input$cutarea),Tm.border=c(as.numeric(input$lowertmborder),as.numeric(input$uppertmborder)),is.deriv=T)
    } else if (input$isderiv == FALSE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=as.numeric(input$cutarea),Tm.border=c(as.numeric(input$lowertmborder),as.numeric(input$uppertmborder)))
    }
  }
  
  ##### b) Render plots for each gene in a different tab
  output$tmplots <- renderUI({
    genes <- SybrGeneList()
    nbgenes <- c(1:length(SybrGeneList()))
    tabs <- lapply(genes, function(x){
      tabPanel(
        title=uiOutput(x), 
        textOutput(paste("text",x,sep="")), 
        uiOutput(paste("down",x,sep="")), 
        plotOutput(paste("tm",x,sep=""), height = "500px")
      )
    })
    
    ## Tab Names
    lapply(genes, function(x){
      output[[x]] <- renderText({
        x 
      })
    })
    
    ## Download Plots (Button)
    lapply(genes, function(x){
      output[[paste("down",x,sep="")]] <- renderUI({
        downloadButton(paste("downl",x,sep=""), "Download PDF")
      })
    })
    
    ## Download Plots (Handler)
    lapply(nbgenes, function(x){
      output[[paste("downl",genes[x],sep="")]] <- downloadHandler(
        filename = (paste(genes[x],"_TMPlots.pdf")),
        content = function(file){
          pdf(file)
          tm <- TMinput()
          genes <- SybrGeneList()
          p <- TMPlots(tm[[x]], genes[[x]])
          dev.off()
        })
    })
    
    ## Render Plots
    lapply(nbgenes, function(x){
      output[[paste("tm",genes[x],sep="")]] <- renderPlot({
        tm <- TMinput()
        TMPlots(tm[[x]], genes[x])
      })
    })
    
    do.call(tabsetPanel,c(tabs))
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
    newidwell_all <- list()
    for (i in 1:length(genes)){
      a<-dim(fluos_gene[[i]])
      b<-a[2]
      c<-as.integer(b+1)
      d<-as.integer(2*b)
      data_gene <- cbind(fluos_gene[[i]],rep.col(melt_gene[[i]]$Temperature,b))
      data_gene <- as.data.frame(data_gene)
      
      if (input$isderiv == TRUE){
        res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=as.numeric(input$cutarea),Tm.border=c(as.numeric(input$lowertmborder),as.numeric(input$uppertmborder)),is.deriv=T)
      } else if (input$isderiv == FALSE){
        res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=as.numeric(input$cutarea),Tm.border=c(as.numeric(input$lowertmborder),as.numeric(input$uppertmborder)))
      }
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
      results_gene$Well <- rownames(results_gene)
      results_gene$Target <- rep(genes[[i]],dim(results_gene)[1])
      
      idwell <- readSYBRIDWell()
      idwell$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', idwell$Well, perl=TRUE))
      results_gene <- merge(results_gene, idwell)
      ord <- c("ID","Target","Well", "Tm", "Area")
      results_gene2 <- results_gene[,ord]
      results_all[[i]] <- results_gene2
      
      ## New IDWell
      newid<- results_gene2[,c("ID", "Target", "Well")]
      newidwell_all[[i]] <- newid
    }
    
    ls <- list()
    results <- list.stack(results_all)
    results <- results[order(as.numeric(as.character(results$ID))),]
    newidwell <- list.stack(newidwell_all)
    newidwell <- newidwell[order(as.numeric(as.character(newidwell$ID))),]
    
    ls[[1]] <- results
    ls[[2]] <- newidwell
    
    return(ls)
  }
  
  ## b) Render table
  output$tmtable <-renderTable({
    TMTable()[[1]]
  })
  
  ## c) Download table
  output$downloadtable <- downloadHandler(
    filename = "TMtable.csv",
    content = function(fname){
      write.csv(TMTable()[[1]], fname, quote = F, row.names = F)}
  )
  
  ################## New ID_Well: SYBR ####################
  output$newidwell <- renderTable({
    TMTable()[[2]]
  })
  
  ##Download table
  output$downnewidwell <- downloadHandler(
    filename = "ID_well.csv",
    content = function(fname){
      write.csv(TMTable()[[2]], fname, quote = F, row.names = F)}
  )
  
  
  'fullpath<-paste("/drive/my-drive", "results_SYBR", sep="/")
  drive_mkdir(fullpath, overwrite = FALSE)#the new folder will be called "results"
  
  output$send2drive <-observeEvent(
    t <- TMTable(),
    drive_upload(t, path = fullpath)
  )'
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ################## BIORAD - SYBR ##################
  ###################################################
  
  ################## Read ID Well Tab: SYBR ##############################
  IDWELLtabSYBR <- function(){
    inp <- input$idwellsybr
    if (!is.null(inp)){
      idwell <- read_csv(inp$datapath)
      idwell <- as.data.frame(idwell)
      idwell$ID <- as.character(idwell$ID)
      idwell <- idwell[order(as.numeric(as.character(idwell$ID))),]
      return(idwell)
    } else{
      return(NULL)
    }
  }
  
  output$IDWELLsybr <- renderTable(
    IDWELLtabSYBR()
  )
  
  ################## Read "Run Information" and "Data": SYBR ##########
  # This function allows the user to upload either 1 Excel with 2 sheets (Data and Run Information)
  # or to upload 2 independent CSVs with the Run Information and Quantification Summary (with Sample,
  # Fluor, Target, Cq, etc columns)
  ## tbl_list will contain:
  ## 1) Df with all data
  ## 2) Df only with IDs to analyse (taken from ID_well from melting curves)
  ## 3) Run Information
  ## 4) Gene List
  readSYBR <- function(){
    res <- input$sybr
    tbl_list <- list()
    if (!is.null(res)){
      if (grepl("csv",res$datapath[[1]][1]) == TRUE){
        tbl_list <- lapply(res$datapath, read.csv)
        tmp <- as.data.frame(tbl_list[1]) %>%
          select(Well, Fluor, Target, Sample, Cq)
        colnames(tmp) <- c("Well", "Fluor", "Target", "ID", "Cq")
        tmp$Cq <- format(tmp$Cq, digits = 4)
        
        newwell <- vector()
        for(i in tmp$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        
        tmp$Well <- newwell
        colnames(tmp) <- c("Well", "Fluor", "Target", "ID", "Cq")
        tmp$ID <- as.character(tmp$ID)
        
        tbl_list[[1]] <- tmp
        
        idwell <- IDWELLtabSYBR()
        new <- vector()
        for(i in idwell$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            new <- append(new, i)
          } else {
            new <- append(new, i)
          }
        }
        idwell$Well <- new
        idwell$ID <- as.character(idwell$ID)
        
        ## Add controls, dils
        ctrls <- tmp[grep("negC|10-|NTC", tmp$ID),]
        m <- tmp[tmp$Well %in% idwell$Well,]
        def <- as.data.frame(rbind(m, ctrls))
        def <- def[order(as.numeric(as.character(def$ID)), def$Target),]
        
        tbl_list[[2]] <- def
        
      } else if (grepl("xlsx", res$datapath[[1]][1]) == TRUE){
        ## Read Data
        dat <- read_xlsx(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target", "ID", "Cq")
        df$Cq <- format(df$Cq, digits = 4)
        
        ## A01 to A1
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df$Well <- newwell
        colnames(df) <- c("Well", "Fluor", "Target", "ID", "Cq")
        df$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', df$Well, perl=TRUE))
        df$ID <- as.character(df$ID)
        
        tbl_list[[1]] <- df
        
        idwell <- IDWELLtabSYBR()
        idwell$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', idwell$Well, perl=TRUE))
        idwell$ID <- as.character(idwell$ID)
        
        ## Add controls, dils
        ctrls <- tmp[grep("negC|10-|NTC", df$ID),]
        m <- tmp[df$Well %in% idwell$Well,]
        def <- as.data.frame(rbind(m, ctrls))
        
        tbl_list[[2]] <- def
        
        
        ## Read Run Information
        run <- read_xlsx(res$datapath, sheet = "Run Information")
        df2 <- as.data.frame(run)
        
        ## Append to list
        tbl_list[[3]] <- df2
        
      } else if (grepl("xls$", res$datapath[[1]][1]) == TRUE){
        ## Read Data
        dat <- read_xls(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target", "ID", "Cq")
        df$Cq <- format(df$Cq, digits = 4)
        
        ## A01 to A1
        newwell <- vector()
        for(i in df$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df$Well <- newwell
        colnames(df) <- c("Well", "Fluor", "Target", "ID", "Cq")
        df$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', df$Well, perl=TRUE))
        df$ID <- as.character(df$ID)
        
        tbl_list[[1]] <- df
        
        idwell <- IDWELLtabSYBR()
        idwell$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', idwell$Well, perl=TRUE))
        idwell$ID <- as.character(idwell$ID)
        
        ## Add controls, dils
        ctrls <- tmp[grep("negC|10-|NTC", df$ID),]
        m <- tmp[df$Well %in% idwell$Well,]
        def <- as.data.frame(rbind(m, ctrls))
        
        tbl_list[[2]] <- def
        
        
        ## Read Run Information
        run <- read_xls(res$datapath, sheet = "Run Information")
        df2 <- as.data.frame(run)
        
        ## Append to list
        tbl_list[[3]] <- df2
        
      }
      
      ## Gene List
      inp <- tbl_list[[1]]
      genes <- unique(inp$Target)
      genes <- as.character(genes)
      genes[genes == ""] <- NA
      genes <- as.character(na.exclude(genes))
      
      # Add to list
      tbl_list[[4]] <- genes
      return(tbl_list)
    }
  }
  
  output$sybrruninfo <- renderTable(
    if (!is.null(readSYBR()[[3]])){
      print(readSYBR()[3])
      readSYBR()[[3]]
    } else if (is.null(readSYBR()[[3]]) == TRUE) {
      "No run information to show"
    }
  )
  
  output$sybrdata <- renderTable(
    readSYBR()[2]
  )
  
  ################## Ct Plate Tab: SYBR ##############
  ctPlateSYBR <- function(){
    inp <- readSYBR()[[1]]
    tofill <- readSYBR()[[2]]
    
    r <- vector()
    c <- vector()
    newwell <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    #ncols <- length(c)
    ncols <- max(as.numeric(c))
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    
    
    for (i in 1:length(tofill$Well)){
      ro <- str_sub(tofill$Well[i],1,1)
      col <- as.numeric(str_sub(tofill$Well[i],-2))
      df[ro,col] <- as.numeric(as.character(tofill$Cq[i]))
    }
    
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printCtPlateSYBR <- function(){
    dt <- ctPlateSYBR()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20))%>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  #### Render Table in App ### 
  output$ctplatesybr<- renderDataTable(
    printCtPlateSYBR()
  )
  
  ################## Sample Plate Tab: SYBR ############################
  samplePlateSYBR <- function(){
    inp <- readSYBR()[[1]]
    tofill <- readSYBR()[[2]]
    
    inp$ID <- as.character(inp$ID)
    r <- vector()
    c <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    #ncols <- length(c)
    ncols <- max(as.numeric(c))
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(tofill$Well)){
      ro <- str_sub(tofill$Well[i],1,1)
      col <- as.numeric(str_sub(tofill$Well[i],-2))
      df[ro,col] <- as.character(tofill$ID[i])
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny ##
  printSamplePlateSYBR <- function(){
    dt <- samplePlateSYBR()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20))
    return(defdt)
  }
  
  ###### Render Table in App 
  output$sampleplatesybr<- renderDataTable(
    printSamplePlateSYBR()
  )
  
  ################## Standard Curve Tab (Table): SYBR #####################
  stCurveSYBR <- function(){
    if (input$posctrlsybr != 0){
      inp <- readSYBR()[[2]]
      genes <- readSYBR()[[4]]
      
      # Get control rows
      ctrls <- lapply(genes, function(x){
        inp[grepl(paste(x,"_",sep=""), inp$ID),]
      })
      ctrls <- bind_rows(ctrls)
      
      ntcs <- lapply(genes, function(x){
        inp[grepl(paste("NTC_",x,sep=""), inp$ID),]
      })
      ntcs <- bind_rows(ntcs)
      
      ctrls <- as.data.frame(rbind(ctrls, ntcs))
      
      # Prepare data frame depending on nb of controls and duplicates
      nbctrls <- as.character(unique(ctrls$ID))
      
      ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
      if(input$dupssybr == TRUE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlsybr), ncol = 3+(length(genes)*2)+length(genes)*3))
        
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(rep(paste(x,"(Dup",1:2,")",sep=""),each=1),paste(x,"(Avg)", sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdsybr)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Gene duplicates and avg
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")] <- tmp$Cq[1]
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")] <- tmp$Cq[2]
              a[tmp_s[2],paste(unique(tmp$Target),"(Avg)",sep="")] <- mean(c(a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")],a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")]), na.rm = TRUE)
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocsybr, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","Fluor","Target","ID","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")] <- negc$Cq[1]
        a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")] <- negc$Cq[2]
        a["C(-)", paste(unique(negc$Target),"(Avg)",sep="")] <- mean(c(a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")],a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")]), na.rm = TRUE)
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        ntc <- bind_rows(lapply(ntc, rbind))
        
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc$Target),]
          
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Fluor","Target","ID","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")] <- df$Cq[1]
          a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")] <- df$Cq[2]
          formean <- na.omit(c(a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")],a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")]))
          a["NTC", paste(unique(df$Target),"(Avg)", sep = "")] <- mean(formean)
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        names(avgs_cp)[names(avgs_cp) == "cp"] <- paste(input$endocsybr,"(Avg)", sep="")
        
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        return(a)
        
        ## If user has no duplicates
      } else if (input$dupssybr == FALSE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlsybr), ncol = 3+length(genes)*3))
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(paste(x,"(Ct)",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdsybr)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Genes control Cq
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Ct)",sep="")] <- tmp$Cq
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocsybr, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","Fluor","Target","ID","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Ct)", sep="")] <- negc$Cq
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc),]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Fluor","Target","ID","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Ct)",sep = "")] <- df$Cq
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("(Ct)", names(a))]
        avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        return(a)
      }
    } else if (input$posctrlsybr == 0){
      return(NULL)
    }
  }
  
  ##### Render Table
  output$stdcurvesybr <- renderTable(
    if (!is.null(stCurveSYBR())){
      stCurveSYBR()
    } else{
      "No standard curve to show" 
    }, rownames = TRUE
  )
  
  ####### Coefficients function #######
  stdCoeffsSYBR <- function(){
    a <- stCurveSYBR()
    if (input$dupssybr == TRUE){
      cp <- a$logCopies
      avgs <- a[,grep("Avg", names(a))]
      avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    } else if (input$dupssybr == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Ct", names(a))]
      avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    }
    
  }
  
  ################## Standard Curve Tab (Plots): SYBR ###########################
  stdPlotsSYBR <- function(){
    if (input$posctrlsybr != 0){
      a <- stCurveSYBR()
      g <- readSYBR()[[2]]
      genes <- readSYBR()[[4]]
      genes <- grep(input$endocsybr, genes, value = TRUE, invert = TRUE)
      
      if(input$dupssybr == TRUE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
        
      } else if (input$dupssybr == FALSE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocsybr, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
      }
    } else if (input$posctrlsybr){
      return(NULL)
    }
    
    
  }
  
  ##### Render Plots in App
  output$stdsybr <- renderPlot(
    if (!is.null(stdPlotsSYBR())){
      stdPlotsSYBR()
    } else {
      NULL
    }
    
  )
  ################## Analysis Samples Tab: SYBR ######################
  AnalysisSamplesSYBR <- function(){
    
    ## Real Id, Sample, Replicate columns
    d <- readSYBR()[[2]]
    genes <- readSYBR()[[4]]
    studygenes <- grep(input$endocsybr, genes, value = TRUE, invert = TRUE)
    d$Cq <- as.numeric(as.character(d$Cq))
    d$ID <- as.character(d$ID)
    for (i in 1:length(d$Cq)){
      if (is.na(d[i, "Target"]) == TRUE){
        next
      } else if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } 
    }
    
    ## If user has duplicates, calculate mean Cq
    if(input$dupssybr == TRUE & input$copiesforassigsybr == TRUE){
      
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          # Check that duplicates exist (2 rows)
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1){
            not0 <- as.numeric(grep("^0$", as.character(y$Cq), invert = TRUE, value = TRUE))
            div <- abs(not0-not0) 
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          # Check if difference between dups is > 1.5
          if (length(y$Cq) == 1){
            if (is.na(y$Cq) == TRUE){
              byidemean <- 0
            } else {
              byidemean <- y$Cq
            }
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1 & div == 0){
            byidemean <- grep("^0", as.character(y$Cq), invert = TRUE, value = TRUE)
            byidemean <- as.numeric(byidemean)
          } else if (div == 0){
            byidmean <- 0
          } else if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
        
        # Remove controls and modify colnames
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove empty df if there is
      mean_cq <- Filter(function(x) dim(x)[1] > 0, mean_cq)
      # Combine dataframes by ID
      mean_cq_merged <- Reduce(function(x, y) merge(x, y, all = TRUE), mean_cq)
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybr)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies
      for (i in studygenes){
        if (all(grepl(i, colnames(mean_cq_merged)) == FALSE)){
          next
        } else if (any(grepl(i, colnames(mean_cq_merged)) == TRUE)){
          l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
            if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
              "Repeat"
            } else if (is.na(as.numeric(paste(x))) == TRUE){
              NA
            } else if (as.numeric(paste(x)) > input$minctsybr & as.numeric(paste(x)) <= input$cyclessybr){
              "Check copy number"
            } else if (as.numeric(paste(x)) <= input$minctsybr & as.numeric(paste(x)) != 0){
              "Positive"
            } else if (as.numeric(paste(x)) > input$cyclessybr){
              "Negative"
            } else if (as.numeric(paste(x)) == 0){
              "Negative"
            }
          })
          mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
        }
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybr,"(MeanCt)", sep = "")], function(x){
        if (is.na(x) == TRUE){
          NA
        } else if (as.character(x) == "Repeat"){
          "Repeat"
        } else if (as.numeric(as.character(x)) >= input$maxendocsybr){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybr){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybr,"(MeanCt)", sep = "")] <- controlcq
      
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgensybr & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgensybr & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgensybr){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## Log Copies
      coeffs <- stdCoeffsSYBR()
      studygenes <- genes[grep(input$endocsybr, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        if (all(grepl(studygenes[i], colnames(mean_cq_merged)) == FALSE)){
          next
        } else {
          for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
            if (is.na(mean_cq_merged[j,"FinalCtCheck"]) == TRUE){
              logcop <- append(logcop, NA)
              cop <- append(cop, NA)
            } else if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
              if (any(grepl(studygenes[i], colnames(mean_cq_merged)) == TRUE)){
                if (is.na(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")]) == TRUE){
                  logcop <- append(logcop, NA)
                  cop <- append(cop, NA)
                } else if (mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")] == 0){
                  logcop <- append(logcop, 0)
                  cop <- append(cop, 0)
                } else {
                  lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
                  if (is.na(lg) == TRUE){
                    logcop <- append(logcop, "-")
                    cop <- append(cop, "-")
                  } else if (is.na(lg) == FALSE){
                    logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                    cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
                  }
                }
              }
            } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
              logcop <- append(logcop, "-")
              cop <- append(cop, "-")
            }
          }
          mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
          mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
        }
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        if (all(grepl(i, colnames(mean_cq_merged)) == FALSE)){
          next
        } else if(any(grepl(i, colnames(mean_cq_merged)) == TRUE)){
          l <- apply(p, MARGIN = 1, function(x){
            if (is.na(x[["FinalCtCheck"]]) == TRUE){
              NA
            } else if (x[["FinalCtCheck"]] == "Check copy number"){
              if (is.na(x[[paste(i,"(Copies)",sep="")]]) == TRUE){
                NA
              } else if (x[[paste(i,"(Copies)",sep="")]] == "-"){
                "-"
              } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
                "Repeat"
              } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvsybr){
                "Positive"
              } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvsybr & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
                "Repeat"
              } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
                "Negative"
              } 
            } else if (x[["FinalCtCheck"]] != "Check copy number"){
              "-"
            }
          })
          mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
        }
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (is.na(x[["FinalCtCheck"]]) == TRUE){
          NA
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgensybr & rep > pos){
            "Repeat"
          } else if (pos >= input$numposgensybr){
            "Positive"
          } else if (pos < input$numposgensybr & neg > pos & neg > rep){
            "Negative"
          } else {
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if (is.na(x[["FinalCtCheck"]]) == TRUE){
          "Repeat"
        } else if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupssybr == FALSE & input$copiesforassigsybr == TRUE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
        #Remove controls
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybr)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctsybr & as.numeric(paste(x)) <= input$cyclessybr){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctsybr & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$cyclessybr){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybr,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybr){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybr){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybr,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgensybr & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgensybr & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgensybr){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## Log Copies
      coeffs <- stdCoeffsSYBR()
      studygenes <- genes[grep(input$endocsybr, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (is.na(mean_cq_merged[j,paste(studygenes[i],"(Ct)", sep="")]) == TRUE){
              logcop <- append(logcop, NA)
              cop <- append(cop, NA)
            } else if (mean_cq_merged[j,paste(studygenes[i],"(Ct)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(Ct)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Ct)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Ct)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } else if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvsybr){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvsybr & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgensybr & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgensybr){
            "Positive"
          }else if (pos < input$numposgensybr & neg > pos & neg > rep){
            "Negative"
          } else {
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupssybr == TRUE & input$copiesforassigsybr == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          # Check that duplicates exist (2 rows)
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1){
            not0 <- as.numeric(grep("^0$", as.character(y$Cq), invert = TRUE, value = TRUE))
            div <- abs(not0-not0) 
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          # Check if difference between dups is > 1.5
          if (length(y$Cq) == 1){
            if (is.na(y$Cq) == TRUE){
              byidemean <- 0
            } else {
              byidemean <- y$Cq
            }
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1 & div == 0){
            byidemean <- grep("^0", as.character(y$Cq), invert = TRUE, value = TRUE)
            byidemean <- as.numeric(byidemean)
          } else if (div == 0){
            byidmean <- 0
          } else if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
        
        # Remove controls and modify colnames
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove empty df if there is
      mean_cq <- Filter(function(x) dim(x)[1] > 0, mean_cq)
      # Combine dataframes by ID
      mean_cq_merged <- Reduce(function(x, y) merge(x, y, all = TRUE), mean_cq)
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybr)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctsybr & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctsybr){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybr,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybr){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybr){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybr,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgensybr){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
      
    } else if (input$dupssybr == FALSE & input$copiesforassigsybr == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
        #Remove controls
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybr)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctsybr & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctsybr){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybr,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybr){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybr){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybr,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgensybr){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybr){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
  }
  
  output$downansybr <- downloadHandler(
    filename = "Analysis.csv",
    content = function(fname){
      write.csv(AnalysisSamplesSYBR(), fname, quote = F, row.names = F)}
  )
  
  AnalysisSamplesDTSYBR <- function(){
    df <- AnalysisSamplesSYBR()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    t <- defdt %>%
      formatStyle(
        columns = grep("Final", names(df)),
        backgroundColor = "khaki"
      ) %>%
      formatStyle(
        columns = "Assignment",
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat"), 
                                     c("limegreen", "tomato", "lightblue")),
        fontWeight = "bold"
      )
    
    return(t)
  }
  
  output$analysissybr <- renderDataTable({
    AnalysisSamplesDTSYBR()
  })
  
  ################## ID Result Tab: SYBR ###################
  idresultSYBR <- function(){
    an <- AnalysisSamplesSYBR()
    realid <- as.numeric(as.character(an$ID))
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
  
  
  
  
  
  
  
  
  
  ################ APPLIED - SYBR ##############
  ##############################################
  
  ################## Read ID Well Tab: SYBR-Applied ##############################
  IDWELLtabSYBRApp <- function(){
    inp <- input$idwellsybrapp
    if (!is.null(inp)){
      idwell <- read_csv(inp$datapath)
      idwell <- as.data.frame(idwell)
      idwell$ID <- as.character(idwell$ID)
      idwell <- idwell[order(as.numeric(as.character(idwell$ID))),]
      return(idwell)
    } else{
      return(NULL)
    }
  }
  
  output$IDWELLsybrapp <- renderTable(
    IDWELLtabSYBRApp()
  )
  
  ################## Read Run Info and Raw Data: SYBR-Applied #####################
  ## tbl_list will contain:
  ## 1) Df with all data
  ## 2) Df only with IDs to analyse (taken from ID_well from melting curves)
  ## 3) Run Information
  ## 4) Gene List
  readAppliedResultsSYBR <- function(){
    res <- input$sybrapp
    tbl_list <- list()
    if (!is.null(res)){
      if (grepl("xlsx", res$datapath[1]) == TRUE){
        samples <- read_xlsx(res$datapath, sheet = "Results", na = "Undetermined")
        df1 <- data.frame(samples)
        start <- grep("Well Position", df1[[2]])
        df1 <- df1[-c(1:start),c(2,4,5,7,9)]
        colnames(df1) <- c("Well","ID", "Target","Fluor","Cq")
        df1$Cq <- format(df1$Cq, digits = 4)
        
        newwell <- vector()
        for(i in df1$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df1$Well <- newwell
        
        # Samples with melting OK
        idwell <- IDWELLtabSYBRApp()
        new <- vector()
        for(i in idwell$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            new <- append(new, i)
          } else {
            new <- append(new, i)
          }
        }
        idwell$Well <- new
        idwell$ID <- as.character(idwell$ID)
        
        ## Add controls, dils
        ctrls <- df1[grep("negC|10-|NTC", df1$ID),]
        m <- df1[df1$Well %in% idwell$Well,]
        
        df2 <- as.data.frame(rbind(m, ctrls))
        
        # Run Information
        df3 <- samples[1:(start-1),1:2]
        df3 <- as.data.frame(df3)
        
        # List of genes
        genes <- read_xlsx(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df4 <- df[-c(1:start),5]
        df4 <- unique(df4)
        
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Samples with melting OK
        tbl_list[[3]] <- df3 # Run Information
        tbl_list[[4]] <- df4 # List of genes
        
        return(tbl_list)
        
      } else if (grepl("xls$", res$datapath[1]) == TRUE){
        samples <- read_xls(res$datapath, sheet = "Results", na = "Undetermined")
        df1 <- data.frame(samples)
        start <- grep("Well Position", df1[[2]])
        df1 <- df1[-c(1:start),c(2,4,5,7,9)]
        colnames(df1) <- c("Well","ID", "Target","Fluor","Cq")
        df1$Cq <- format(df1$Cq, digits = 4)
        
        newwell <- vector()
        for(i in df1$Well){
          if (grepl("^[[:upper:]][[:digit:]]$",i) == TRUE){
            stri_sub(i, 2, 1) <- 0
            newwell <- append(newwell, i)
          } else {
            newwell <- append(newwell, i)
          }
        }
        df1$Well <- newwell
        
        # Samples with melting OK
        idwell <- IDWELLtabSYBRApp()
        idwell$Well <- as.character(sub('(?<![0-9])0*(?=[0-9])', '', idwell$Well, perl=TRUE))
        idwell$ID <- as.character(idwell$ID)
        ## Add controls, dils
        ctrls <- df1[grep("negC|10-|NTC", df1$ID),]
        m <- df1[df1$Well %in% idwell$Well,]
        df2 <- as.data.frame(rbind(m, ctrls))
        
        # Run Information
        df3 <- samples[1:(start-1),1:2]
        df3 <- as.data.frame(df3)
        
        # List of genes
        genes <- read_xls(res$datapath, sheet = "Results")
        df <- as.data.frame(genes)
        df4 <- df[-c(1:start),5]
        df4 <- unique(df4)
        
        ## Append to list
        tbl_list[[1]] <- df1 # Data
        tbl_list[[2]] <- df2 # Samples with melting OK
        tbl_list[[3]] <- df3 # Run Information
        tbl_list[[4]] <- df4 # List of genes
        return(tbl_list)
      }
    }
    
  }
  
  output$sybrdataapp <- renderTable({
    readAppliedResultsSYBR()[[2]]
  })
  
  output$sybrruninfoapp <- renderTable({
    if (!is.null(readAppliedResultsSYBR()[[3]])){
      readAppliedResultsSYBR()[[3]]
    } else {
      "No run information to show"
    }
  })
  ################## Ct Plate Tab: SYBR-Applied ##################
  ctPlateAppSYBR <- function(){
    inp <- readAppliedResultsSYBR()[[1]]
    tofill <- readAppliedResultsSYBR()[[2]]
    
    r <- vector()
    c <- vector()
    newwell <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    #ncols <- length(c)
    ncols <- max(as.numeric(c))
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    
    for (i in 1:length(tofill$Well)){
      ro <- str_sub(tofill$Well[i],1,1)
      col <- as.numeric(str_sub(tofill$Well[i],-2))
      df[ro,col] <- as.numeric(as.character(tofill$Cq[i]))
    }
    
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printCtPlateAppSYBR <- function(){
    dt <- ctPlateAppSYBR()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20)) %>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  output$ctplatesybrapp <- renderDataTable(
    printCtPlateAppSYBR()
  )
  
  ################## Sample Plate Tab: SYBR-Applied ##########################
  samplePlateAppSYBR <- function(){
    inp <- readAppliedResultsSYBR()[[1]]
    tofill <- readAppliedResultsSYBR()[[2]]
    
    inp$ID <- as.character(inp$ID)
    r <- vector()
    c <- vector()
    for (i in inp$Well){
      r <- append(r,str_sub(i,1,1))
      c <- append(c, str_sub(i,-2))
    }
    r <-unique(r)
    c <- unique(c)
    nrows <- length(r)
    #ncols <- length(c)
    ncols <- max(as.numeric(c))
    
    df <- data.frame(matrix(nrow = nrows, ncol = ncols))
    colnames(df) <- c(1:ncols)
    rownames(df) <- r
    for (i in 1:length(tofill$Well)){
      ro <- str_sub(tofill$Well[i],1,1)
      col <- as.numeric(str_sub(tofill$Well[i],-2))
      df[ro,col] <- as.character(tofill$ID[i])
    }
    return(df)
  }
  
  ## Convert to datatable to render nicely in shiny
  printSamplePlateAppSYBR <- function(){
    dt <- samplePlateAppSYBR()
    defdt <- datatable(dt, rownames = T,  options = list(pageLength = 20))
    return(defdt)
  }
  
  output$sampleplatesybrapp<- renderDataTable(
    printSamplePlateAppSYBR()
  )
  
  ################## Standard Curve Tab (Table): SYBR-Applied #####################
  stCurveAppSYBR <- function(){
    if (input$posctrlsybrapp != 0){
      inp <- readAppliedResultsSYBR()[[2]]
      genes <- readAppliedResultsSYBR()[[4]]
      
      # Get control rows
      ctrls <- lapply(genes, function(x){
        inp[grepl(paste(x,"_",sep=""), inp$ID),]
      })
      ctrls <- bind_rows(ctrls)
      
      ntcs <- lapply(genes, function(x){
        inp[grepl(paste("NTC_",x,sep=""), inp$ID),]
      })
      ntcs <- bind_rows(ntcs)
      
      ctrls <- as.data.frame(rbind(ctrls, ntcs))
      
      # Prepare data frame depending on nb of controls and duplicates
      nbctrls <- as.character(unique(ctrls$ID))
      
      ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
      if(input$dupssybrapp == TRUE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlsybrapp), ncol = 3+(length(genes)*2)+length(genes)*3))
        
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(rep(paste(x,"(Dup",1:2,")",sep=""),each=1),paste(x,"(Avg)", sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdsybrapp)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Gene duplicates and avg
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")] <- tmp$Cq[1]
              a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")] <- tmp$Cq[2]
              a[tmp_s[2],paste(unique(tmp$Target),"(Avg)",sep="")] <- mean(c(a[tmp_s[2],paste(unique(tmp$Target),"(Dup1)",sep="")],a[tmp_s[2],paste(unique(tmp$Target),"(Dup2)",sep="")]), na.rm = TRUE)
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocsybrapp, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","ID","Target","Fluor","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")] <- negc$Cq[1]
        a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")] <- negc$Cq[2]
        a["C(-)", paste(unique(negc$Target),"(Avg)",sep="")] <- mean(c(a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")],a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")]), na.rm = TRUE)
        
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        ntc <- bind_rows(lapply(ntc, rbind))
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc$Target),]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","Fluor","Target","ID","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")] <- df$Cq[1]
          a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")] <- df$Cq[2]
          formean <- na.omit(c(a["NTC", paste(unique(df$Target),"(Dup1)", sep = "")],a["NTC", paste(unique(df$Target),"(Dup2)", sep = "")]))
          a["NTC", paste(unique(df$Target),"(Avg)", sep = "")] <- mean(formean)
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        return(a)
        
        ## If user has no duplicates
      } else if (input$dupssybrapp == FALSE){
        a <- data.frame(matrix(NA, nrow = 2+as.numeric(input$posctrlsybrapp), ncol = 3+length(genes)*3))
        #Prepare colnames
        colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
        l <- lapply(genes, function(x){
          list(paste(x,"(Ct)",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
        })
        l <- unlist(l)
        colnames(a)[4:ncol(a)] <- l
        
        ## Sample
        smp <- lapply(nbctrls, function(x){str_split(x,"_")})
        smp <- unlist(smp)
        smp <- as.character(unique(grep("10-", smp,value = TRUE)))
        smp <- smp[order(nchar(smp), smp)]
        rownames(a) <- c("NTC", "C(-)", smp)
        
        ## Dilution
        dils <- lapply(smp, function(x){str_split(x,"-")})
        dilsnum <- lapply(dils, function(x){as.numeric(x[[1]])})
        dilsfin <- lapply(dilsnum, function(x){x[[1]]^x[[2]]})
        dilsfin <- unlist(dilsfin)
        dilsfin <- sort(dilsfin)
        a$Dilution <- c(NA, NA, dilsfin)
        
        ##Copies
        copies<- lapply(a$Dilution, function(x){as.numeric(input$concstdsybrapp)/x})
        copies <- unlist(copies)
        a$Copies <- as.numeric(copies)
        
        ##log(Copies)
        logcopies <- lapply(a$Copies, function(x){log(x, base = 10)})
        logcopies <- unlist(logcopies)
        a$logCopies <- as.numeric(logcopies)
        
        ## Duplicates
        ctrls$Target <- as.character(ctrls$Target)
        ctrls$ID <- as.character(ctrls$ID)
        ctrls$Cq <- as.character(ctrls$Cq)
        ctrls_sp <- split(ctrls, ctrls$Target)
        
        ## Genes control Cq
        for (dil in rownames(a)[3:length(rownames(a))]){
          for (df in ctrls_sp){
            tmp <- df[grep(dil,df$ID),]
            tmp$Cq <- as.numeric(tmp$Cq)
            tmp_s <- str_split(tmp$ID[1], "_")
            tmp_s <- unlist(tmp_s)
            if (!is.na(tmp_s[1])){
              a[tmp_s[2],paste(unique(tmp$Target),"(Ct)",sep="")] <- tmp$Cq
            }
          }
        }
        
        ## Negative control
        control <- ctrls_sp[grep(input$endocsybrapp, ctrls_sp)]
        control <- as.data.frame(control)
        colnames(control) <- c("Well","Fluor","Target","ID","Cq")
        negc <- control[grep("negC", control$ID),]
        negc$Cq <- as.numeric(negc$Cq)
        a["C(-)", paste(unique(negc$Target),"(Ct)", sep="")] <- negc$Cq
        
        ## NTC
        ntc <- lapply(ctrls_sp, function(x){
          x[grep("NTC", x$ID),]
        })
        
        for (i in 1:length(genes)){
          df <- ntc[grep(genes[i], ntc)]
          df <- as.data.frame(df)
          colnames(df) <- c("Well","ID","Target","Fluor","Cq")
          df$Cq <- as.numeric(df$Cq)
          a["NTC", paste(unique(df$Target),"(Ct)",sep = "")] <- df$Cq
        }
        
        ## Coefficients
        cp <- a$logCopies
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        
        avgs_cp <- cbind(avgs_def, cp)
        d <- data.frame()
        coefficients <- sapply(avgs_cp, function(x){
          model <- lm(cp ~ x, avgs_cp)
          coeff1 <- as.numeric(model$coefficients[1])
          coeff2 <- as.numeric(model$coefficients[2])
          coeffs <- as.data.frame(cbind(coeff1, coeff2))
        })
        rownames(coefficients) <- c("Intercept", "Slope")
        neworder <- order(as.character(colnames(coefficients)))
        coefficients2 <- coefficients[,neworder]
        
        ##### LogCopies and Copies
        for (i in 1:nrow(a)){
          for (j in 1:length(avgs)){
            if(is.na(avgs[i,j]) == TRUE){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else if (avgs[i,j] == 0){
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
            } else {
              a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients2[2,j])+unlist(coefficients2[1,j])
            }
            if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
              a[i, paste(genes[j],"(Copies)", sep="")] <- NA
            } else{
              a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
            }
          }
        }
        
        return(a)
      }
    } else if (input$posctrlsybrapp == 0){
      return(NULL)
    }
  }
  
  ##### Render Table
  output$stdcurvesybrapp <- renderTable(
    if (!is.null(stCurveAppSYBR())){
      stCurveAppSYBR()
    } else {
      "No standard curve to show" 
    }, rownames = TRUE
  )
  
  ####### Coefficients function #######
  stdCoeffsAppSYBR <- function(){
    a <- stCurveAppSYBR()
    if (input$dupssybrapp == TRUE){
      cp <- a$logCopies
      avgs <- a[,grep("Avg", names(a))]
      avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    } else if (input$dupssybrapp == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Ct", names(a))]
      avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
      avgs_def <- avgs[,avgs_genes]
      avgs_cp <- cbind(avgs_def, cp)
      d <- data.frame()
      coefficients <- sapply(avgs_cp, function(x){
        model <- lm(cp ~ x, avgs_cp)
        coeff1 <- as.numeric(model$coefficients[1])
        coeff2 <- as.numeric(model$coefficients[2])
        coeffs <- as.data.frame(cbind(coeff1, coeff2))
      })
      rownames(coefficients) <- c("Intercept", "Slope")
      neworder <- order(as.character(colnames(coefficients)))
      coefficients2 <- coefficients[,neworder]
      return(coefficients2)
    }
    
  }
  
  ################## Standard Curve Tab (Plots): SYBR-Applied ###########################
  stdPlotsAppSYBR <- function(){
    if (input$posctrlsybrapp != 0){
      a <- stCurveAppSYBR()
      g <- readAppliedResultsSYBR()[[2]]
      genes <- readAppliedResultsSYBR()[[4]]
      genes <- grep(input$endocsybrapp, genes, value = TRUE, invert = TRUE)
      
      if(input$dupssybrapp == TRUE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Avg", names(a))]
        avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
        
      } else if (input$dupssybrapp == FALSE){
        ## Get copies column
        cp <- a$logCopies
        ## Get avgs column (only for not control genes)
        avgs <- a[,grep("Ct", names(a))]
        avgs_genes <- grep(input$endocsybrapp, names(avgs), value = TRUE, invert = TRUE)
        avgs_def <- avgs[,avgs_genes]
        ## Combine copies with avgs
        avgs_cp <- cbind(avgs_def, cp)
        avgs_cp <- avgs_cp[3:nrow(avgs_cp),]
        
        doStdPlots <- function(gene, nam){
          ggplot(avgs_cp, aes(x=gene, y = cp)) + 
            geom_point()+
            geom_smooth(method = lm, se = F) +
            stat_poly_eq(formula = cp ~ gene,
                         label.x.npc = "right", label.y.npc = "top",
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                         parse = TRUE) +
            ggtitle(paste("Standard curve for",nam))+
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
            xlab(paste("Dilutions",nam)) +
            ylab("log(Copies)")
        }
        
        p <- map2(avgs_cp[grep("cp", names(avgs_cp), value = TRUE, invert = TRUE)], genes, doStdPlots)
        ggarrange(plotlist = p, nrow = 1)
      }
    } else if (input$posctrlsybrapp == 0){
      return(NULL)
    }
  }
  
  ##### Render Plots in App
  output$stdsybrapp <- renderPlot(
    stdPlotsAppSYBR()
  )
  
  ################## Analysis Samples Tab: SYBR-Applied ######################
  AnalysisSamplesAppSYBR <- function(){
    
    ## Real Id, Sample, Replicate columns
    d <- readAppliedResultsSYBR()[[2]]
    genes <- readAppliedResultsSYBR()[[4]]
    studygenes <- grep(input$endocsybrapp, genes, value = TRUE, invert = TRUE)
    d$Cq <- as.numeric(d$Cq)
    d$ID <- as.character(d$ID)
    for (i in 1:length(d$Cq)){
      if (is.na(d[i, "Target"]) == TRUE){
        next
      } else if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } 
    }
    
    ## If user has duplicates, calculate mean Cq
    if(input$dupssybrapp == TRUE & input$copiesforassigsybrapp == TRUE){
      
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          # Check that duplicates exist (2 rows)
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1){
            not0 <- as.numeric(grep("^0$", as.character(y$Cq), invert = TRUE, value = TRUE))
            div <- abs(not0-not0) 
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          # Check if difference between dups is > 1.5
          if (length(y$Cq) == 1){
            if (is.na(y$Cq) == TRUE){
              byidemean <- 0
            } else {
              byidemean <- y$Cq
            }
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1 & div == 0){
            byidemean <- grep("^0", as.character(y$Cq), invert = TRUE, value = TRUE)
            byidemean <- as.numeric(byidemean)
          } else if (div == 0){
            byidmean <- 0
          } else if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
        
        # Remove controls and modify colnames
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove empty df if there is
      mean_cq <- Filter(function(x) dim(x)[1] > 0, mean_cq)
      # Combine dataframes by ID
      mean_cq_merged <- Reduce(function(x, y) merge(x, y, all = TRUE), mean_cq)
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybrapp)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies
      for (i in studygenes){
        if (all(grepl(i, colnames(mean_cq_merged)) == FALSE)){
          next
        } else if (any(grepl(i, colnames(mean_cq_merged)) == TRUE)){
          l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
            if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
              "Repeat"
            } else if (is.na(as.numeric(paste(x))) == TRUE){
              NA
            } else if (as.numeric(paste(x)) > input$minctsybrapp & as.numeric(paste(x)) <= input$cyclessybrapp){
              "Check copy number"
            } else if (as.numeric(paste(x)) <= input$minctsybrapp & as.numeric(paste(x)) != 0){
              "Positive"
            } else if (as.numeric(paste(x)) > input$cyclessybrapp){
              "Negative"
            } else if (as.numeric(paste(x)) == 0){
              "Negative"
            }
          })
          mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
        }
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybrapp,"(MeanCt)", sep = "")], function(x){
        if (is.na(x) == TRUE){
          NA
        } else if (as.character(x) == "Repeat"){
          "Repeat"
        } else if (as.numeric(as.character(x)) >= input$maxendocsybrapp){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybrapp){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybrapp,"(MeanCt)", sep = "")] <- controlcq
      
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgensybrapp & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgensybrapp & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## Log Copies
      coeffs <- stdCoeffsAppSYBR()
      studygenes <- genes[grep(input$endocsybrapp, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        if (all(grepl(studygenes[i], colnames(mean_cq_merged)) == FALSE)){
          next
        } else {
          for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
            if (is.na(mean_cq_merged[j,"FinalCtCheck"]) == TRUE){
              logcop <- append(logcop, NA)
              cop <- append(cop, NA)
            } else if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
              if (any(grepl(studygenes[i], colnames(mean_cq_merged)) == TRUE)){
                if (is.na(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")]) == TRUE){
                  logcop <- append(logcop, NA)
                  cop <- append(cop, NA)
                } else if (mean_cq_merged[j,paste(studygenes[i],"(MeanCt)", sep="")] == 0){
                  logcop <- append(logcop, 0)
                  cop <- append(cop, 0)
                } else {
                  lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(MeanCt)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
                  if (is.na(lg) == TRUE){
                    logcop <- append(logcop, "-")
                    cop <- append(cop, "-")
                  } else if (is.na(lg) == FALSE){
                    logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                    cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
                  }
                }
              }
            } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
              logcop <- append(logcop, "-")
              cop <- append(cop, "-")
            }
          }
          mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
          mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
        }
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        if (all(grepl(i, colnames(mean_cq_merged)) == FALSE)){
          next
        } else if(any(grepl(i, colnames(mean_cq_merged)) == TRUE)){
          l <- apply(p, MARGIN = 1, function(x){
            if (is.na(x[["FinalCtCheck"]]) == TRUE){
              NA
            } else if (x[["FinalCtCheck"]] == "Check copy number"){
              if (is.na(x[[paste(i,"(Copies)",sep="")]]) == TRUE){
                NA
              } else if (x[[paste(i,"(Copies)",sep="")]] == "-"){
                "-"
              } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
                "Repeat"
              }  else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvsybrapp){
                "Positive"
              } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvsybrapp & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
                "Repeat"
              } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
                "Negative"
              } 
            } else if (x[["FinalCtCheck"]] != "Check copy number"){
              "-"
            }
          })
          mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
        }
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (is.na(x[["FinalCtCheck"]]) == TRUE){
          NA
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgensybrapp & rep > pos){
            "Repeat"
          } else if (pos >= input$numposgensybrapp){
            "Positive"
          } else if (pos < input$numposgensybrapp & neg > pos & neg > rep){
            "Negative"
          } else {
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if (is.na(x[["FinalCtCheck"]]) == TRUE){
          "Repeat"
        } else if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupssybrapp == FALSE & input$copiesforassigsybrapp == TRUE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
        #Remove controls
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybrapp)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) > input$minctsybrapp & as.numeric(paste(x)) <= input$cyclessybrapp){
            "Check copy number"
          } else if (as.numeric(paste(x)) <= input$minctsybrapp & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$cyclessybrapp){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybrapp,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybr){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybrapp){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybrapp,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          check <- sum(str_count(ctchecks[i,], "Check copy number"), na.rm = TRUE)
          if (check >= input$numposgensybrapp & check >= pos){
            ctassig <- append(ctassig, "Check copy number")
          } else if (check < input$numposgensybrapp & check >= pos & check >= neg & check != 0){
            ctassig <- append(ctassig, "Check copy number")
          } else if(pos >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > check){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## Log Copies
      coeffs <- stdCoeffsAppSYBR()
      studygenes <- genes[grep(input$endocsybrapp, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in 1:length(mean_cq_merged[,"FinalCtCheck"])){
          if (mean_cq_merged[j,"FinalCtCheck"] == "Check copy number"){
            if (is.na(mean_cq_merged[j,paste(studygenes[i],"(Ct)", sep="")]) == TRUE){
              logcop <- append(logcop, NA)
              cop <- append(cop, NA)
            } else if (mean_cq_merged[j,paste(studygenes[i],"(Ct)", sep="")] == 0){
              logcop <- append(logcop, 0)
              cop <- append(cop, 0)
            } else {
              lg <- as.numeric(as.character(mean_cq_merged[j,paste(studygenes[i],"(Ct)",sep = "")]))*unlist(coeffs[2,paste(studygenes[i],"(Ct)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Ct)",sep="")])
              if (is.na(lg) == TRUE){
                logcop <- append(logcop, "-")
                cop <- append(cop, "-")
              } else if (is.na(lg) == FALSE){
                logcop <- append(logcop, format(round(lg, 3), nsmall = 3))
                cop <- append(cop, format(round(as.numeric(10^lg), 3), nsmall = 3))
              }
            }
          } else if (mean_cq_merged[j,"FinalCtCheck"] != "Check copy number"){
            logcop <- append(logcop, "-")
            cop <- append(cop, "-")
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (is.na(as.numeric(x[[paste(i,"(Copies)",sep="")]])) == TRUE){
              NA
            } else if (x[[paste(i,"(Copies)",sep="")]] == "-"){
              "-"
            } else if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) >= input$mincnvsybrapp){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvsybrapp & as.numeric(x[[paste(i,"(Copies)",sep="")]]) != 0){
              "Repeat"
            } 
          } else if (x[["FinalCtCheck"]] != "Check copy number"){
            "-"
          }
        })
        mean_cq_merged[[paste("CopyCheck:",i, sep="")]] <- l
      }
      
      ## Final assignation according to copy number
      copychecks <- mean_cq_merged[,grep("CopyCheck|FinalCtCheck", names(mean_cq_merged))]
      
      copyassig <- apply(copychecks, MARGIN = 1, function(x){
        if (x[["FinalCtCheck"]] == "Check copy number"){
          cc <- x[grep("CopyCheck", names(x))]
          pos <- sum(str_count(cc, "Positive"), na.rm = TRUE)
          neg <- sum(str_count(cc, "Negative"), na.rm = TRUE)
          rep <- sum(str_count(cc, "Repeat"), na.rm = TRUE)
          if (rep >= input$numposgensybrapp & rep > pos){
            "Repeat"
          }else if (pos >= input$numposgensybrapp){
            "Positive"
          }else if (pos < input$numposgensybrapp & neg > pos & neg > rep){
            "Negative"
          } else {
            "Repeat"
          }
        } else if (x[["FinalCtCheck"]] != "Check copy number"){
          "-"
        }
      })
      mean_cq_merged[["FinalCopyCheck"]] <- copyassig
      
      ## Assignment
      final <- mean_cq_merged[,grep("Final", names(mean_cq_merged))]
      finalassig <- apply(final, MARGIN = 1, function(x){
        if(x[["FinalCtCheck"]] != "Check copy number"){
          x[["FinalCtCheck"]]
        } else if (x[["FinalCtCheck"]] == "Check copy number"){
          x[["FinalCopyCheck"]]
        }
      })
      
      mean_cq_merged[["Assignment"]] <- finalassig
      return(mean_cq_merged)
      
    } else if (input$dupssybrapp == TRUE & input$copiesforassigsybrapp == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      # Split by ID
      # Check if Cq1/Cq2 is bigger than 1.5
      # If not, calculate mean, else Repeat sample
      mean_cq <- lapply(d_g_clean, function(x){
        id <- split(x, x$ID)
        byidmean <- lapply(id, function(y){
          # Check that duplicates exist (2 rows)
          if (is.na(y$Cq[1]) == TRUE & is.na(y$Cq[2]) == TRUE){
            div <- NA
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1){
            not0 <- as.numeric(grep("^0$", as.character(y$Cq), invert = TRUE, value = TRUE))
            div <- abs(not0-not0) 
          } else {
            div <- abs(y$Cq[1]-y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          
          # Check if difference between dups is > 1.5
          if (length(y$Cq) == 1){
            if (is.na(y$Cq) == TRUE){
              byidemean <- 0
            } else {
              byidemean <- y$Cq
            }
          } else if (sum(str_count(as.character(y$Cq), pattern = "^0$")) == 1 & div == 0){
            byidemean <- grep("^0", as.character(y$Cq), invert = TRUE, value = TRUE)
            byidemean <- as.numeric(byidemean)
          } else if (div == 0){
            byidmean <- 0
          } else if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else {
            byidemean <- format(round(mean(c(y$Cq[1], y$Cq[2])), 3), nsmall = 3)
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
        
        # Remove controls and modify colnames
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCt)", sep=""))
      }
      
      # Remove empty df if there is
      mean_cq <- Filter(function(x) dim(x)[1] > 0, mean_cq)
      # Combine dataframes by ID
      mean_cq_merged <- Reduce(function(x, y) merge(x, y, all = TRUE), mean_cq)
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybrapp)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCt)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctsybrapp & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctsybrapp){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybrapp,"(MeanCt)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybrapp){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybrapp){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybrapp,"(MeanCt)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
      
    } else if (input$dupssybrapp == FALSE & input$copiesforassigsybrapp == FALSE){
      # Split by gene
      d_g <- split(d, d$Target)
      # Keep only df with genes names (if no value is given in ID, a new "empty" df is built)
      d_g_clean <- list()
      for (i in 1:length(genes)){
        d_g_clean[i] <- d_g[grep(genes[i], names(d_g))]
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- as.data.frame(cbind(d_g_clean[[i]]$ID, d_g_clean[[i]]$Cq))
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
        #Remove controls
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Ct)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Remove samples with more NAs than the min number to assign a sample as POS
      toremove <- vector()
      for (i in 1:nrow(mean_cq_merged)){
        na <- as.numeric(sum(is.na(mean_cq_merged[i,])))
        if (na >= as.numeric(input$numposgensybrapp)){
          toremove <- append(toremove, as.numeric(as.character(mean_cq_merged$ID[i])))
        } else {
          next
        }
      }
      mean_cq_merged <- mean_cq_merged[!mean_cq_merged$ID %in% toremove,]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Ct)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) <= input$minctsybrapp & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctsybrapp){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      #### Repeat for wrong controls (> maxctendoc)
      controlcq <- sapply(mean_cq_merged[,paste(input$endocsybrapp,"(Ct)", sep = "")], function(x){
        if (as.numeric(as.character(x)) >= input$maxendocsybrapp){
          "Repeat"
        } else if (as.numeric(as.character(x)) < input$maxendocsybrapp){
          as.numeric(as.character(x))
        }
      })
      mean_cq_merged[,paste(input$endocsybrapp,"(Ct)", sep = "")] <- controlcq
      
      ## Final assignation according to Ct checks: FINAL ASSIGNMENT
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
          neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
          unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
          if(pos >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Positive")
          } else if (as.numeric(unc) >= as.numeric(pos) & as.numeric(unc) >= as.numeric(neg) & as.numeric(unc) >= input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == unc & unc < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg == pos & neg < input$numposgensybrapp){
            ctassig <- append(ctassig, "Repeat")
          } else if (neg > pos & neg > unc){
            ctassig <- append(ctassig, "Negative")
          } 
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
  }
  
  output$downansybrapp <- downloadHandler(
    filename = "Analysis.csv",
    content = function(fname){
      write.csv(AnalysisSamplesAppSYBR(), fname, quote = F, row.names = F)}
  )
  
  AnalysisSamplesDTAppSYBR <- function(){
    df <- AnalysisSamplesAppSYBR()
    defdt <- datatable(df, rownames = F, options = list(pageLength = 120))
    t <- defdt %>%
      formatStyle(
        columns = grep("Final", names(df)),
        backgroundColor = "khaki"
      ) %>%
      formatStyle(
        columns = "Assignment",
        backgroundColor = styleEqual(c("Positive", "Negative", "Repeat"), 
                                     c("limegreen", "tomato", "lightblue")),
        fontWeight = "bold"
      )
    
    return(t)
  }
  
  output$analysissybrapp <- renderDataTable({
    AnalysisSamplesDTAppSYBR()
  })
  
  ################## ID Result Tab: SYBR-Applied ###################
  idresultAppSYBR <- function(){
    an <- AnalysisSamplesAppSYBR()
    realid <- unique(as.numeric(as.character(an$ID)))
    final <- as.data.frame(cbind(realid, an$Assignment))
    colnames(final)<- c("ID", "Interpretation")
    return(final)
  }
  
  output$downIDRESsybrapp <- downloadHandler(
    filename = "ID_result.csv",
    content = function(fname){
      write.csv(idresultAppSYBR(), fname, quote = F, row.names = F)}
  )  
  
  output$IDRESsybrapp <- renderTable(
    idresultAppSYBR()
  )
  
}