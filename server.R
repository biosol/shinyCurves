################################# SHINY APP (COMPLETE) FOR QPCR ANALYSIS ############################################
########################### Sonia Olaechea-Lázaro (UPV/EHU, May 2020) ############################################  

server <- function(input, output) {
  
  ################## BIORAD ##################
  ############################################
  
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
          select(Well, Fluor, Target, Content,Sample, Cq)
        colnames(tmp) <- c("Well", "Fluor", "Target", "Content", "ID", "Cq")
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
        colnames(tmp) <- c("Well", "Fluor", "Target", "Content", "ID", "Cq")
        tmp$Cq <- format(tmp$Cq, digits = 4)
        
        tbl_list[[1]] <- tmp
        
        ## Read Data from 1 Excel
      } else if (grepl("xlsx", res$datapath[1]) == TRUE){
        dat <- read_xlsx(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Content ,df$Sample, df$Cq)
        df <- as.data.frame(df)
        
        ## Correct A1 to A01, etc
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
        colnames(df) <- c("Well", "Fluor", "Target", "Content","Sample", "Cq")
        df$Cq <- format(df$Cq, digits = 4)
        
        
        ## Read Run Information
        run <- read_xlsx(res$datapath, sheet = "Run Information")
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
  
  ################## Cq Plate Tab: Biorad ######################
  cqPlate <- function(){
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
  printCqPlate <- function(){
    dt <- cqPlate()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20)) %>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  #### Render Table in App ### 
  output$cqplate<- renderDataTable(
    printCqPlate()
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
    inp <- readBiorad()[[1]]
    genes <- readBiorad()[[3]]
    
    # Get control rows
    ctrls <- lapply(genes, function(x){
      inp[grepl(x, inp$ID),]
      })
    ctrls <- bind_rows(ctrls)
    
    # Prepare data frame depending on nb of controls and duplicates
    nbctrls <- as.character(unique(ctrls$ID))
    
    ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
    if(input$dupsbiorad == TRUE){
      a <- data.frame(matrix(0, nrow = 2+as.numeric(input$posctrlbiorad), ncol = 3+(length(genes)*2)+length(genes)*3))
      
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
      copies<- lapply(a$Dilution, function(x){200000*2/x})
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
      colnames(control) <- c("Well","Fluor","Target","Content","ID","Cq")
      negc <- control[grep("negC", control$ID),]
      negc$Cq <- as.numeric(negc$Cq)
      a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")] <- negc$Cq[1]
      a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")] <- negc$Cq[2]
      a["C(-)", paste(unique(negc$Target),"(Avg)",sep="")] <- mean(c(a["C(-)", paste(unique(negc$Target),"(Dup1)",sep="")],a["C(-)", paste(unique(negc$Target),"(Dup2)",sep="")]), na.rm = TRUE)
      
      ## NTC
      ntc <- lapply(ctrls_sp, function(x){
        x[grep("NTC", x$ID),]
      })
      
      for (i in 1:length(genes)){
        df <- ntc[grep(genes[i], ntc)]
        df <- as.data.frame(df)
        colnames(df) <- c("Well","Fluor","Target","Content","ID","Cq")
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
      
      ##### LogCopies and Copies
      for (i in 1:nrow(a)){
        for (j in 1:length(avgs)){
          if(is.na(avgs[i,j]) == TRUE){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
          } else if (avgs[i,j] == 0){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- 0
          } else {
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients[2,j])+unlist(coefficients[1,j])
          }
          if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
            a[i, paste(genes[j],"(Copies)", sep="")] <- NA
          } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
            a[i, paste(genes[j],"(Copies)", sep="")] <- 0
          } else{
            a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
          }
        }
      }
      return(a)
      
      ## If user has no duplicates
    } else if (input$dupsbiorad == FALSE){
      a <- data.frame(matrix(0, nrow = 2+as.numeric(input$posctrlbiorad), ncol = 3+length(genes)*3))
      #Prepare colnames
      colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
      l <- lapply(genes, function(x){
        list(paste(x,"Cq",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
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
      copies<- lapply(a$Dilution, function(x){200000*2/x})
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
            a[tmp_s[2],paste(unique(tmp$Target),"Cq",sep="")] <- tmp$Cq
          }
        }
      }
      
      ## Negative control
      control <- ctrls_sp[grep(input$endocbiorad, ctrls_sp)]
      control <- as.data.frame(control)
      colnames(control) <- c("Well","Fluor","Target","Content","ID","Cq")
      negc <- control[grep("negC", control$ID),]
      negc$Cq <- as.numeric(negc$Cq)
      a["C(-)", paste(unique(negc$Target),"Cq", sep="")] <- negc$Cq
      
      ## NTC
      ntc <- lapply(ctrls_sp, function(x){
        x[grep("NTC", x$ID),]
      })
      
      for (i in 1:length(genes)){
        df <- ntc[grep(genes[i], ntc)]
        df <- as.data.frame(df)
        colnames(df) <- c("Well","Fluor","Target","Content","ID","Cq")
        df$Cq <- as.numeric(df$Cq)
        a["NTC", paste(unique(df$Target),"Cq",sep = "")] <- df$Cq
      }
      
      ## Coefficients
      cp <- a$logCopies
      avgs <- a[,grep("Cq", names(a))]
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
      
      ##### LogCopies and Copies
      for (i in 1:nrow(a)){
        for (j in 1:length(avgs)){
          if(is.na(avgs[i,j]) == TRUE){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
          } else if (avgs[i,j] == 0){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- 0
          } else {
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients[2,j])+unlist(coefficients[1,j])
          }
          if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
            a[i, paste(genes[j],"(Copies)", sep="")] <- NA
          } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
            a[i, paste(genes[j],"(Copies)", sep="")] <- 0
          } else{
            a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
          }
        }
      }
      
      return(a)
    }
    
  }
  
  #### Convert to datatable
  stCurveDT <- function(){
    a <- stCurve()
    ## Transform into data table
    defdt <- datatable(a, rownames = TRUE)
    dt <- defdt %>%
      formatRound(
        columns = c(1:ncol(a)),
        digits = 3
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
      return(coefficients)
    } else if (input$dupsbiorad == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Cq", names(a))]
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
      return(coefficients)
    }

  }
  
  ################## Standard Curve Tab (Plots): Biorad ##########################
  stdPlots <- function(){
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
      avgs <- a[,grep("Cq", names(a))]
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

  }
  
  ##### Render Plots in App
  output$std <- renderPlot(
    stdPlots()
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
      if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } else{
        next
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
            div <- abs(y$Cq[1]/y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else if (div < 0.5 & div > 0){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- mean(c(y$Cq[1], y$Cq[2]))
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
 
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
         mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Log Copies
      coeffs <- stdCoeffs()
      studygenes <- genes[grep(input$endoC, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in mean_cq_merged[,paste(studygenes[i],"(MeanCq)",sep = "")]){
          if (j == "Repeat"){
            logcop <- append(logcop, "Repeat")
            cop <- append(cop, "Repeat")
          } else if(is.na(j) == TRUE){
            logcop <- append(logcop, NA)
            cop <- append(cop, NA)
          } else if(j == 0){
            logcop <- append(logcop, 0)
            cop <- append(cop, 0)
          } else {
            lg <- as.numeric(j)*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
            logcop <- append(logcop, lg)
            cop <- append(cop, 10^lg)
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$rangebiorad[1] & as.numeric(paste(x)) < input$rangebiorad[2]){
            "Check copy number"
          } else if (as.numeric(paste(x)) > input$rangebiorad[2]){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }

      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          if (any(grepl("Check copy number", ctchecks[i,])) == TRUE){
            ctassig <- append(ctassig, "Check copy number")
          } else if (any(grepl("Check copy number", ctchecks[i,]) == FALSE)){
            pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
            neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
            unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
            if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
              ctassig <- append(ctassig, "Repeat")
            } else if(pos >= input$numposgenbiorad){
              ctassig <- append(ctassig, "Positive")
            } else if (pos < input$numposgenbiorad){
              ctassig <- append(ctassig, "Negative")
            } else if (neg > pos){
              ctassig <- append(ctassig, "Negative")
            }
          }
        }
      }
      
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
                "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) > input$mincnvbiorad){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvbiorad){
              "Repeat"
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
          if (pos >= input$numposgenbiorad){
            "Positive"
          } else if (pos < input$numposgenbiorad){
            "Negative"
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
        colnames(d_g_clean[[i]]) <- c("ID", "Cq")
      }
      
      for (i in 1:length(d_g_clean)){
        d_g_clean[[i]] <- d_g_clean[[i]][grep("10-|NTC|negC|PosCtrl", d_g_clean[[i]]$ID, value = TRUE ,invert = TRUE),]
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Cq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Log Copies
      coeffs <- stdCoeffs()
      studygenes <- genes[grep(input$endoC, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in mean_cq_merged[,paste(studygenes[i],"(Cq)",sep = "")]){
          if (j == "Repeat"){
            logcop <- append(logcop, "Repeat")
            cop <- append(cop, "Repeat")
          } else if(is.na(j) == TRUE){
            logcop <- append(logcop, NA)
            cop <- append(cop, NA)
          } else if(j == 0){
            logcop <- append(logcop, 0)
            cop <- append(cop, 0)
          } else {
            lg <- as.numeric(j)*unlist(coeffs[2,paste(studygenes[i],"Cq",sep="")])+unlist(coeffs[1,paste(studygenes[i],"Cq",sep="")])
            logcop <- append(logcop, lg)
            cop <- append(cop, 10^lg)
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Cq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$rangebiorad[1] & as.numeric(paste(x)) < input$rangebiorad[2]){
            "Check copy number"
          } else if (as.numeric(paste(x)) > input$rangebiorad[2]){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          if (any(grepl("Check copy number", ctchecks[i,])) == TRUE){
            ctassig <- append(ctassig, "Check copy number")
          } else if (any(grepl("Check copy number", ctchecks[i,]) == FALSE)){
            pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
            neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
            unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
            if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
              ctassig <- append(ctassig, "Repeat")
            } else if(pos >= input$numposgenbiorad){
              ctassig <- append(ctassig, "Positive")
            } else if (pos < input$numposgenbiorad){
              ctassig <- append(ctassig, "Negative")
            } else if (neg > pos){
              ctassig <- append(ctassig, "Negative")
            }
          }
        }
      }
      
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) > input$mincnvbiorad){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvbiorad){
              "Repeat"
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
          if (pos >= input$numposgenbiorad){
            "Positive"
          } else if (pos < input$numposgenbiorad){
            "Negative"
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
            div <- abs(y$Cq[1]/y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else if (div < 0.5 & div > 0){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- mean(c(y$Cq[1], y$Cq[2]))
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctbiorad){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
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
            if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
              ctassig <- append(ctassig, "Repeat")
            } else if(pos >= input$numposgenbiorad){
              ctassig <- append(ctassig, "Positive")
            } else if (pos < input$numposgenbiorad){
              ctassig <- append(ctassig, "Negative")
            } else if (neg > pos){
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
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Cq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Cq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctbiorad & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctbiorad){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
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
          if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
            ctassig <- append(ctassig, "Repeat")
          } else if(pos >= input$numposgenbiorad){
            ctassig <- append(ctassig, "Positive")
          } else if (pos < input$numposgenbiorad){
            ctassig <- append(ctassig, "Negative")
          } else if (neg > pos){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
    
  }
      
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
      ) %>%
      formatRound(
        columns = grep("Copies", names(df)),
        digits = 3
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
  
    
  ################## ID Result Tab: Biorad ###############################
  idresult <- function(){
    an <- AnalysisSamples()
    
    id <-unique(as.numeric(as.character(an$ID)))
    idres <- as.data.frame(cbind(id, an$Assignment))
    colnames(idres)<- c("ID", "Interpretation")
    idres <- idres[which(is.na(idres$Interpretation) == FALSE),]
    #Remove Repeat from ID_Result (might give problems with the curves)
    #idres <- idres[grep("Repeat", idres$Interpretation, invert = TRUE),]
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
  
  ################## Genes Tabs: Applied ######################
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
        ls <- getGenes()
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
          write.csv(getGenes()[x], fname, quote = F, row.names = F)
        }
      )
    })
    
    ## Render Tables
    lapply(nbgenes, function(x){
      output[[paste("trans",genes[x],sep="")]] <- renderTable({
          getGenes()[x]
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
  
  ################## Read Results: Applied #####################
  readAppliedResults <- function(){
    res <- input$appl
    samples <- read_xlsx(res$datapath, sheet = "Results", na = "Undetermined")
    df <- data.frame(samples)
    df <- df[-c(1:44),c(2,4,5,7,9)]
    colnames(df) <- c("Well","ID", "Target","Fluor","Cq")
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
  
  output$appliedres <- renderTable({
    readAppliedResults()
  })
  
  ################## Cq Plate Tab: Applied ##################
  cqPlateApp <- function(){
    inp <- readAppliedResults()
    r <- vector()
    c <- vector()
    for (i in inp$Well){
        r <- append(r,str_sub(i,1,1))
        c <- append(c, str_sub(i,-1))
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
  printCqPlateApp <- function(){
    dt <- cqPlateApp()
    defdt <- datatable(dt, rownames = T, options = list(pageLength = 20)) %>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  output$cqplateapp<- renderDataTable(
    printCqPlateApp()
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
    defdt <- datatable(dt, rownames = F,  options = list(pageLength = 20))
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
    inp <- readAppliedResults()
    genes <- readApplied()[[3]]
    
    # Get control rows
    ctrls <- lapply(genes, function(x){
      inp[grepl(x, inp$ID),]
    })
    ctrls <- bind_rows(ctrls)
    
    # Prepare data frame depending on nb of controls and duplicates
    nbctrls <- as.character(unique(ctrls$ID))
    
    ## If user has duplicates (also using the # of positive controls that the users has specified to generate # of lines for table)
    if(input$dupsapplied == TRUE){
      a <- data.frame(matrix(0, nrow = 2+as.numeric(input$posctrlapplied), ncol = 3+(length(genes)*2)+length(genes)*3))
      
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
      copies<- lapply(a$Dilution, function(x){200000*2/x})
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
      
    
      for (i in 1:length(genes)){
        df <- ntc[grep(genes[i], ntc)]
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
      
      ##### LogCopies and Copies
      for (i in 1:nrow(a)){
        for (j in 1:length(avgs)){
          if(is.na(avgs[i,j]) == TRUE){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
          } else if (avgs[i,j] == 0){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- 0
          } else {
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients[2,j])+unlist(coefficients[1,j])
          }
          if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
            a[i, paste(genes[j],"(Copies)", sep="")] <- NA
          } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
            a[i, paste(genes[j],"(Copies)", sep="")] <- 0
          } else{
            a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
          }
        }
      }
      return(a)
      
      ## If user has no duplicates
    } else if (input$dupsapplied == FALSE){
      a <- data.frame(matrix(0, nrow = 2+as.numeric(input$posctrlapplied), ncol = 3+length(genes)*3))
      #Prepare colnames
      colnames(a)[1:3] <- c("Dilution","Copies","logCopies")
      l <- lapply(genes, function(x){
        list(paste(x,"Cq",sep=""),paste(x,"(LogCopies)",sep=""),paste(x,"(Copies)",sep=""))
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
      copies<- lapply(a$Dilution, function(x){200000*2/x})
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
            a[tmp_s[2],paste(unique(tmp$Target),"Cq",sep="")] <- tmp$Cq
          }
        }
      }
      
      ## Negative control
      control <- ctrls_sp[grep(input$endocapplied, ctrls_sp)]
      control <- as.data.frame(control)
      colnames(control) <- c("Well","ID","Target","Fluor","Cq")
      negc <- control[grep("negC", control$ID),]
      negc$Cq <- as.numeric(negc$Cq)
      a["C(-)", paste(unique(negc$Target),"Cq", sep="")] <- negc$Cq
      
      ## NTC
      ntc <- lapply(ctrls_sp, function(x){
        x[grep("NTC", x$ID),]
      })
      
      for (i in 1:length(genes)){
        df <- ntc[grep(genes[i], ntc)]
        df <- as.data.frame(df)
        colnames(df) <- c("Well","ID","Target","Fluor","Cq")
        df$Cq <- as.numeric(df$Cq)
        a["NTC", paste(unique(df$Target),"Cq",sep = "")] <- df$Cq
      }
      
      ## Coefficients
      cp <- a$logCopies
      avgs <- a[,grep("Cq", names(a))]
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
      
      ##### LogCopies and Copies
      for (i in 1:nrow(a)){
        for (j in 1:length(avgs)){
          if(is.na(avgs[i,j]) == TRUE){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- NA
          } else if (avgs[i,j] == 0){
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- 0
          } else {
            a[i, paste(genes[j],"(LogCopies)", sep="")] <- avgs[i,j]*unlist(coefficients[2,j])+unlist(coefficients[1,j])
          }
          if (is.na(a[i, paste(genes[j],"(LogCopies)", sep="")]) == TRUE){
            a[i, paste(genes[j],"(Copies)", sep="")] <- NA
          } else if(a[i, paste(genes[j],"(LogCopies)", sep="")] == 0){
            a[i, paste(genes[j],"(Copies)", sep="")] <- 0
          } else{
            a[i, paste(genes[j],"(Copies)", sep="")] <- 10^a[i, paste(genes[j],"(LogCopies)", sep="")]
          }
        }
      }
      
      return(a)
    }
    
  }
  
  #### Convert to datatable
  stCurveDTApp <- function(){
    a <- stCurveApp()
    ## Transform into data table
    defdt <- datatable(a, rownames = TRUE)
    dt <- defdt %>%
      formatRound(
        columns = c(1:ncol(a)),
        digits = 3
      )
    return(dt)
  }
  
  ##### Render Table
  output$stdcurveapp <- renderDataTable(
    stCurveDTApp()
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
      return(coefficients)
    } else if (input$dupsapplied == FALSE){
      cp <- a$logCopies
      avgs <- a[,grep("Cq", names(a))]
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
      return(coefficients)
    }
    
  }
  
  ################## Standard Curve Tab (Plots): Applied ###########################
  stdPlotsApp <- function(){
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
      avgs <- a[,grep("Cq", names(a))]
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
      if(d[i,"Target"] != "" & is.na.data.frame(d[i,"Cq"]) == TRUE){
        d$Cq[i] <- 0
      } else if (d[i,"Target"] == ""){
        d$Cq[i] <- NA
      } else{
        next
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
            div <- abs(y$Cq[1]/y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else if (div < 0.5 & div > 0){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- mean(c(y$Cq[1], y$Cq[2]))
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Log Copies
      coeffs <- stdCoeffsApp()
      studygenes <- genes[grep(input$endocapplied, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in mean_cq_merged[,paste(studygenes[i],"(MeanCq)",sep = "")]){
          if (j == "Repeat"){
            logcop <- append(logcop, "Repeat")
            cop <- append(cop, "Repeat")
          } else if(is.na(j) == TRUE){
            logcop <- append(logcop, NA)
            cop <- append(cop, NA)
          } else if(j == 0){
            logcop <- append(logcop, 0)
            cop <- append(cop, 0)
          } else {
            lg <- as.numeric(j)*unlist(coeffs[2,paste(studygenes[i],"(Avg)",sep="")])+unlist(coeffs[1,paste(studygenes[i],"(Avg)",sep="")])
            logcop <- append(logcop, lg)
            cop <- append(cop, 10^lg)
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$rangeapplied[1] & as.numeric(paste(x)) < input$rangeapplied[2]){
            "Check copy number"
          } else if (as.numeric(paste(x)) > input$rangeapplied[2]){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          if (any(grepl("Check copy number", ctchecks[i,])) == TRUE){
            ctassig <- append(ctassig, "Check copy number")
          } else if (any(grepl("Check copy number", ctchecks[i,]) == FALSE)){
            pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
            neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
            unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
            if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
              ctassig <- append(ctassig, "Repeat")
            } else if(pos >= input$numposgenapplied){
              ctassig <- append(ctassig, "Positive")
            } else if (pos < input$numposgenapplied){
              ctassig <- append(ctassig, "Negative")
            } else if (neg > pos){
              ctassig <- append(ctassig, "Negative")
            }
          }
        }
      }
      
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) > input$mincnvapplied){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvapplied){
              "Repeat"
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
          if (pos >= input$numposgenapplied){
            "Positive"
          } else if (pos < input$numposgenapplied){
            "Negative"
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
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Cq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Log Copies
      coeffs <- stdCoeffsApp()
      studygenes <- genes[grep(input$endocapplied, genes, invert = TRUE)]
      df <- data.frame(matrix(0,ncol = length(studygenes), nrow = nrow(mean_cq_merged[1])))
      for (i in 1:length(studygenes)){
        logcop <- vector()
        cop <- vector()
        for (j in mean_cq_merged[,paste(studygenes[i],"(Cq)",sep = "")]){
          if (j == "Repeat"){
            logcop <- append(logcop, "Repeat")
            cop <- append(cop, "Repeat")
          } else if(is.na(j) == TRUE){
            logcop <- append(logcop, NA)
            cop <- append(cop, NA)
          } else if(j == 0){
            logcop <- append(logcop, 0)
            cop <- append(cop, 0)
          } else {
            lg <- as.numeric(j)*unlist(coeffs[2,paste(studygenes[i],"Cq",sep="")])+unlist(coeffs[1,paste(studygenes[i],"Cq",sep="")])
            logcop <- append(logcop, lg)
            cop <- append(cop, 10^lg)
          }
        }
        mean_cq_merged[[paste(studygenes[i],"(LogCopies)", sep="")]] <- logcop
        mean_cq_merged[[paste(studygenes[i],"(Copies)", sep="")]] <- cop
      }
      
      mean_cq_merged[mean_cq_merged == "NA"] <- NA
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Cq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) > input$rangeapplied[1] & as.numeric(paste(x)) < input$rangeapplied[2]){
            "Check copy number"
          } else if (as.numeric(paste(x)) > input$rangeapplied[2]){
            "Negative"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
      ## Final assignation according to Ct checks
      ctchecks <- mean_cq_merged[,grep("CtCheck", names(mean_cq_merged))]
      ctassig <- vector()
      for (i in 1:nrow(ctchecks)){
        # NAs
        if (all(is.na(ctchecks[i,])) == TRUE){
          ctassig <- append(ctassig, NA)
        } else if (all(is.na(ctchecks[i,])) == FALSE){
          if (any(grepl("Check copy number", ctchecks[i,])) == TRUE){
            ctassig <- append(ctassig, "Check copy number")
          } else if (any(grepl("Check copy number", ctchecks[i,]) == FALSE)){
            pos <- sum(str_count(ctchecks[i,], "Positive"), na.rm = TRUE)
            neg <- sum(str_count(ctchecks[i,], "Negative"), na.rm = TRUE)
            unc <- sum(str_count(ctchecks[i,], "Repeat"), na.rm = TRUE)
            if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
              ctassig <- append(ctassig, "Repeat")
            } else if(pos >= input$numposgenapplied){
              ctassig <- append(ctassig, "Positive")
            } else if (pos < input$numposgenapplied){
              ctassig <- append(ctassig, "Negative")
            } else if (neg > pos){
              ctassig <- append(ctassig, "Negative")
            }
          }
        }
      }
      
      mean_cq_merged[["FinalCtCheck"]] <- ctassig
      
      ## CopyCheck:gene columns
      ## Copy number (only for duplicate samples) and Ct in range specified (35-40 for Covid)
      p <- as.data.frame(cbind(mean_cq_merged[,grep("FinalCtCheck|(Copies)", names(mean_cq_merged))]))
      
      for (i in studygenes){
        l <- apply(p, MARGIN = 1, function(x){
          if (x[["FinalCtCheck"]] == "Check copy number"){
            if (x[[paste(i,"(Copies)",sep="")]] == "Repeat"){
              "Repeat"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) == 0){
              "Negative"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) > input$mincnvapplied){
              "Positive"
            } else if (as.numeric(x[[paste(i,"(Copies)",sep="")]]) < input$mincnvapplied){
              "Repeat"
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
          if (pos >= input$numposgenapplied){
            "Positive"
          } else if (pos < input$numposgenapplied){
            "Negative"
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
            div <- abs(y$Cq[1]/y$Cq[2])
            if (is.na(div) == TRUE){
              div <- 0
            }
          }
          
          if (is.na(div) == TRUE){
            byidmean <- NA
          } else if (div > 1.5){
            byidemean <- "Repeat"
          } else if (div < 0.5 & div > 0){
            byidemean <- "Repeat"
          } else if (div == 0){
            byidmean <- 0
          } else {
            byidemean <- mean(c(y$Cq[1], y$Cq[2]))
          }
        })
        n <- names(byidmean)
        byidmean <- cbind(n, as.data.frame(melt(as.character(byidmean))))
      })
      
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- as.data.frame(mean_cq[[i]])
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Remove controls and modify colnames
      for (i in 1:length(mean_cq)){
        mean_cq[[i]] <- mean_cq[[i]][grep("10-|NTC|negC|PosCtrl", mean_cq[[i]]$ID, invert = TRUE),]
        colnames(mean_cq[[i]]) <- c("ID", paste(genes[i],"(MeanCq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(mean_cq, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(MeanCq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctapplied){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
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
          if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
            ctassig <- append(ctassig, "Repeat")
          } else if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (pos < input$numposgenapplied){
            ctassig <- append(ctassig, "Negative")
          } else if (neg > pos){
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
        colnames(d_g_clean[[i]]) <- c("ID", paste(genes[i],"(Cq)", sep=""))
      }
      
      # Combine dataframes by ID
      mean_cq_merged <- join_all(d_g_clean, by="ID")
      mean_cq_merged <- mean_cq_merged[order(as.numeric(as.character(mean_cq_merged$ID))),]
      
      ## Positive samples (< than minctbiorad) only for study genes
      ## If mean Repeat -> Repeat
      ## if mean na -> na
      ## if value in range specified -> look copies 
      for (i in studygenes){
        l <- sapply(mean_cq_merged[,paste(i,"(Cq)",sep="")], function(x){
          if(as.character(x) == "Repeat" & is.na(as.character(x)) == FALSE){
            "Repeat"
          } else if (is.na(as.numeric(paste(x))) == TRUE){
            NA
          } else if (as.numeric(paste(x)) < input$minctapplied & as.numeric(paste(x)) != 0){
            "Positive"
          } else if (as.numeric(paste(x)) == 0){
            "Negative"
          } else if (as.numeric(paste(x)) > input$minctapplied){
            "Negative"
          }
        })
        mean_cq_merged[[paste("CtCheck:",i,sep="")]] <- l
      }
      
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
          if (as.numeric(unc) > as.numeric(pos) & as.numeric(unc) > as.numeric(neg)){
            ctassig <- append(ctassig, "Repeat")
          } else if(pos >= input$numposgenapplied){
            ctassig <- append(ctassig, "Positive")
          } else if (pos < input$numposgenapplied){
            ctassig <- append(ctassig, "Negative")
          } else if (neg > pos){
            ctassig <- append(ctassig, "Negative")
          }
        }
      }
      
      mean_cq_merged[["Assignment"]] <- ctassig
      return(mean_cq_merged)
    }
    
  }
  
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
      ) %>%
      formatRound(
        columns = grep("Copies", names(df)),
        digits = 3
      )
    
    return(t)
  }
  
  output$analysisapp <- renderDataTable({
    AnalysisSamplesDTApp()
  })
  
  ################## Interpretation Tab: Applied ####################
  'interpretationApp <- function(){
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
                                     c("seagreen", "tomato", "lightblue", "lemonchiffon"))
      )
    return(f)
  }
  
  ### Render in App ###
  output$interpretapp <- renderDataTable(
    interpretationAppDT()
  )'
  
  ################## ID Well Tab: Applied ###############
  IDWELLtabApp <- function(){
    inp <- readAppliedResults()
    wellid <- cbind(inp$Well, as.character(inp$Target), inp$ID)
    wellid <- as.data.frame(wellid)
    print(wellid)
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
    #int <- interpretationApp()
    realid <- unique(as.numeric(as.character(an$ID)))
    final <- as.data.frame(cbind(realid, an$Assignment))
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
      
      if (length(merge_r_persample) == 0){
        mes <- "No indetermined plots to show"
        return(mes)
      } else {
      ## Plots
      pltList <- list()
      for (z in 1:length(merge_r_persample)){
        pltName <- paste(genes[i],"_",names(merge_r_persample)[[z]], sep = "")
        pltList[[pltName]] <- ggplot(data = merge_pn, aes(x=Cycles, y=RFU, by=Well2, color=Interpretation)) +
          geom_line()+
          geom_line(data = merge_r_persample[[z]], aes(x=Cycles, y=RFU, by=Well2, color=Interpretation))+
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
  
  'output$indetgen1 <- renderUI({
    ls <- indetPlots()
    print(str(ls))
    if (is.character(ls) == TRUE){
      textOutput("indetgen1")
      print(ls)
    } else {
      plot_grid(plotlist = ls[[1]], ncol = 2)
    }
  })
  
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
  output$indetgen2 <- renderUI({
    ls <- indetPlots()
    
    if (is.vector(ls) == TRUE){
      textOutput("indetgen2")
      print(ls)
    } else {
      plot_grid(plotlist = ls[[2]], ncol = 2)
    }
  })
  
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
  output$indetgen3 <- renderUI({
    ls <- indetPlots()
    if (is.vector(ls) == TRUE){
      textOutput("indetgen3")
      print(ls)
    } else {
      plot_grid(plotlist = ls[[3]], ncol = 2)
    }
  })
  
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
  })'

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
        colnames(tmp) <- c("Well", "Fluor", "Target", "Content", "ID", "Cq")
        tbl_list[[1]] <- tmp
        
      } else if (grepl("xlsx", res$datapath[[1]][1]) == TRUE){
        ## Read Data
        dat <- read_xlsx(res$datapath, sheet = "Data")
        df <- as.data.frame(dat)
        df <- cbind(df$Well, df$Fluor, df$Target, df$Content ,df$Sample, df$Cq)
        df <- as.data.frame(df)
        colnames(df) <- c("Well", "Fluor", "Target", "Content","ID", "Cq")
        
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
  printCqPlateSYBR <- function(){
    dt <- cqPlateSYBR()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))%>%
      formatRound(columns = c(1:ncol(dt)), digits = 3)
    return(defdt)
  }
  
  #### Render Table in App ### 
  output$cqplatesybr<- renderDataTable(
    printCqPlateSYBR()
  )
  
  ################## Sample Plate Tab: SYBR ############################
  samplePlateSYBR <- function(){
    inp <- readSYBR()[[1]]
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
  printSamplePlateSYBR <- function(){
    dt <- samplePlateSYBR()
    defdt <- datatable(dt, rownames = F, options = list(pageLength = 20))
    return(defdt)
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
      tbl_list <- lapply(dat$datapath, read.csv, header=TRUE, sep=";", dec=",")
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
    f <- function(df){
      out <- as.character(as.numeric(df))
    }
    data_gene <- as.data.frame(data_gene)
    
    if (input$isderiv == TRUE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=as.numeric(input$cutarea),Tm.border=c(as.numeric(input$lowertmborder),as.numeric(input$uppertmborder)),is.deriv=T)
    } else if (input$isderiv == FALSE){
      res_gene <- meltcurve(data_gene,temps=c(c:d),fluos=c(1:b),cut.Area=input$cutarea,Tm.border=c(input$lowertmborder,input$uppertmborder))
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
  
  'output$sybrN <- renderPlot({
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
    })'
  
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
