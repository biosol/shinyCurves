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
ui <- fluidPage(
  titlePanel(strong("qPCR Analysis")),
  sidebarLayout(
    sidebarPanel(
      h5(strong("Taqman")),
      checkboxInput("taqexcel", "Analysis - BioRad"),
      conditionalPanel(
        condition = "input.taqexcel == true",
        fileInput("biorad", 
                  label = "Upload BioRad CFX results file",
                  multiple = TRUE),
        textInput("endocbiorad", "Endogenous control", value = "RNAseP"),
        numericInput("maxendocbiorad", "Maximum cycle number for endogenous control", value = 35),
        numericInput("posctrlbiorad", "How many SERIAL DILUTIONS are you using for the standard curve?", value = 4),
        numericInput("minctbiorad", "Sample is considered POSITIVE with Ct below:", value = 35),
        radioButtons("dupsbiorad", "Do you use duplicates?", choices = c(Yes = TRUE, No = FALSE), selected = TRUE),
        numericInput("numposgenbiorad", "Number of positive genes to consider a sample POSITIVE", value = 1),
        radioButtons("copiesforassig", "Do you want to use the estimated copy number as a result assignation criteria?", choices = c(Yes=TRUE, No=FALSE), selected = TRUE),
        conditionalPanel(condition = "input.copiesforassig == 'TRUE'",
          numericRangeInput("rangebiorad", "Enter range:", value = c(35,40)),
          numericInput("mincnvbiorad", "Sample is considered POSITIVE with estimated copies value above:", value = 4)
          )
        ),
      checkboxInput("app", "Analysis - Applied Quant Studio"),
      conditionalPanel(
        condition = "input.app == true",
        fileInput("appl",
                  label = "Upload Applied Quant Studio xls/xlsx here",
                  multiple = TRUE),
        textInput("endocapplied", "Endogenous control", value = "RNAseP"),
        numericInput("maxendocapplied", "Maximum cycle number for endogenous control", value = 35),
        numericInput("posctrlapplied", "How many SERIAL DILUTIONS are you using for the standard curve?", value = 4),
        numericInput("minctapplied", "Sample is considered POSITIVE with Ct below:", value = 35),
        radioButtons("dupsapplied", "Do you use duplicates?", choices = c(Yes = TRUE, No = FALSE), selected = TRUE),
        numericInput("numposgenapplied", "Number of positive genes to consider a sample POSITIVE", value = 1),
        radioButtons("copiesforassigapplied", "Do you want to use the estimated copy number as a result assignation criteria?", choices = c(Yes=TRUE, No=FALSE), selected = TRUE),
        conditionalPanel(condition = "input.copiesforassigapplied == 'TRUE'",
                         numericRangeInput("rangeapplied", "Enter range:", value = c(35,40)),
                         numericInput("mincnvapplied", "Sample is considered POSITIVE with estimated copies value above:", value = 4)
        )
      ),
      checkboxInput("taqman", "Ct Curves"),
      conditionalPanel(
        condition = "input.taqman == true",
        fileInput("taqmancsv",
                  label = "Upload CSVs here",
                  multiple = TRUE),
        numericInput("ct", "Enter cycle number", value = 40),
        textInput("endoC", "Enter endogenous control", value = "RNAseP"),
        fileInput("taqwellid", label = "Upload ID_well.csv here"),
        fileInput("taqidres", label = "Upload ID_result.csv here")
      ),
      h5(strong("SYBR Green")),
      checkboxInput("inpsybr", "Analysis - Biorad"),
      conditionalPanel(
        condition = "input.inpsybr == true",
        fileInput("sybr",
                  label = "Upload Excel or CSVs here",
                  multiple = TRUE),
        textInput("endocsybr", "Endogenous control", value = "RNAseP"),
        numericInput("maxendocsybr", "Maximum cycle number for endogenous control", value = 35),
        numericInput("posctrlsybr", "How many SERIAL DILUTIONS are you using for the standard curve?", value = 4),
        numericInput("minctsybr", "Sample is considered POSITIVE with Ct below:", value = 35),
        radioButtons("dupssybr", "Do you use duplicates?", choices = c(Yes = TRUE, No = FALSE), selected = TRUE),
        numericInput("numposgensybr", "Number of positive genes to consider a sample POSITIVE", value = 1),
        radioButtons("copiesforassigsybr", "Do you want to use the estimated copy number as a result assignation criteria?", choices = c(Yes=TRUE, No=FALSE), selected = TRUE),
        conditionalPanel(condition = "input.copiesforassigsybr == 'TRUE'",
                         numericRangeInput("rangesybr", "Enter range:", value = c(35,40)),
                         numericInput("mincnvsybr", "Sample is considered POSITIVE with estimated copies value above:", value = 4)
        )
      ),
      checkboxInput("inpsybrapp", "Analysis - Applied Quant Studio"),
      conditionalPanel(
        condition = "input.inpsybrapp == true",
        fileInput("sybrapp",
                  label = "Upload Applied Quant Studio xls/xlsx here",
                  multiple = TRUE),
        textInput("endocsybrapp", "Endogenous control", value = "RNAseP"),
        numericInput("maxendocsybrapp", "Maximum cycle number for endogenous control", value = 35),
        numericInput("posctrlsybrapp", "How many SERIAL DILUTIONS are you using for the standard curve?", value = 4),
        numericInput("minctsybrapp", "Sample is considered POSITIVE with Ct below:", value = 35),
        radioButtons("dupssybrapp", "Do you use duplicates?", choices = c(Yes = TRUE, No = FALSE), selected = TRUE),
        numericInput("numposgensybrapp", "Number of positive genes to consider a sample POSITIVE", value = 1),
        radioButtons("copiesforassigsybrapp", "Do you want to use the estimated copy number as a result assignation criteria?", choices = c(Yes=TRUE, No=FALSE), selected = TRUE),
        conditionalPanel(condition = "input.copiesforassigsybrapp == 'TRUE'",
                         numericRangeInput("rangesybrapp", "Enter range:", value = c(35,40)),
                         numericInput("mincnvsybrapp", "Sample is considered POSITIVE with estimated copies value above:", value = 4)
        )
      ),
      checkboxInput("meltsybr", "Melting Curves"),
      conditionalPanel(
        condition = "input.meltsybr == true",
        fileInput("sybrcsv",
                  label="Upload CSV here",
                  multiple = TRUE),
        numericInput("cutarea", "Enter cut.Area value", value = 10),
        textInput("lowertmborder", "Enter lower Tm border", value = 0.5),
        textInput("uppertmborder", "Enter upper Tm border", value = 0.5),
        radioButtons("isderiv", "Are your data already in first derivative transformed format?",
                     choices = c(Yes = TRUE,
                                 No = FALSE),
                     selected = TRUE),
        fileInput("wellid",
                  label="Upload ID_well.csv here")
      ),
      downloadButton("downtoy", "Download Toy Dataset", icon("Documentos/Covid19/shiny-app/qPCR/www/virus-solid.svg")),
      downloadButton("downmanual", "Download Manual")
    ),
      mainPanel(
      ###### 1a) Main Panel for Biorad Input ########
      conditionalPanel(
        condition = "input.taqexcel == true",
        tabsetPanel(
          id = "taqex",
          type = "tabs",
          tabPanel("Raw Data", tableOutput("inputdf")),
          tabPanel("Run Information", tableOutput("BRruninfo")),
          tabPanel("Ct Plate", dataTableOutput("ctplate")),
          tabPanel("Sample Plate", dataTableOutput("sampleplate")),
          tabPanel("Standard Curve", fluidRow(tableOutput("stdcurve"), plotOutput("std"))),
          tabPanel("Analysis", downloadButton("downanbiorad", "Download CSV"),dataTableOutput("analysis")),
          tabPanel("ID_Well", downloadButton("downIDWELL", "Download CSV"), tableOutput("IDWELL")),
          tabPanel("ID_Result",downloadButton("downidres", "Download CSV"), tableOutput("IDRESULT"))
        )
      ),
      ###### 1b) Main Panel for Applied Input ########
      conditionalPanel(
        condition = "input.app == true",
        tabsetPanel(
          id = "blabla",
          type = "tabs",
          tabPanel("Raw Data", tableOutput("readapp")),
          tabPanel("Conversion", uiOutput("conversionapp")),
          tabPanel("Run Info", tableOutput("appliedruninfo")),
          tabPanel("Applied Results", tableOutput("appliedres")),
          tabPanel("Ct Plate", dataTableOutput("ctplateapp")),
          tabPanel("Sample Plate", dataTableOutput("sampleplateapp")),
          tabPanel("Standard Curve", tableOutput("stdcurveapp"), plotOutput("stdapp")),
          tabPanel("Analysis", downloadButton("downanapplied", "Download CSV"),dataTableOutput("analysisapp")),
          tabPanel("ID Well", downloadButton("downidwellapp", "Download CSV"),tableOutput("IDWELLApp")),
          tabPanel("ID Result", downloadButton("downidresapp", "Download CSV"),tableOutput("IDRESULTApp"))
        )
      ),
      ###### 2) Main Panel for Taqman ########
      conditionalPanel(
        condition = "input.taqman == true",
        tabsetPanel(
          id = "taq",
          type = "tabs",
          tabPanel("ID_Well", tableOutput("taqmanwell")),
          tabPanel("ID_Result", tableOutput("taqmanidres")),
          tabPanel("General Plots", downloadButton("downgen", "Download PDF"), plotOutput("genplots", height = "1000px")),
          tabPanel("Indet Plots", uiOutput("indetplots"))
        )
      ),
      ###### 3) Main Panel for BIORAD - SYBR ########
      conditionalPanel(
        condition = "input.inpsybr == true",
        tabsetPanel(
          id = "inpsyb",
          type = "tabs",
          tabPanel("Run Information", tableOutput("sybrruninfo")),
          tabPanel("Raw Data", tableOutput("sybrdata")),
          tabPanel("Ct Plate", dataTableOutput("ctplatesybr")),
          tabPanel("Sample Plate", dataTableOutput("sampleplatesybr")),
          tabPanel("Standard Curve", tableOutput("stdcurvesybr"), plotOutput("stdsybr")),
          tabPanel("Analysis", downloadButton("downansybr", "Download CSV"), dataTableOutput("analysissybr")),
          tabPanel("ID_Well", downloadButton("downIDWELLsybr", "Download CSV"), tableOutput("IDWELLsybr")),
          tabPanel("ID_Result", downloadButton("downIDRESsybr", "Download CSV"), tableOutput("IDRESsybr"))
        )
      ),
      ############ SYBR - APPLIED #########
      conditionalPanel(
        condition = "input.inpsybrapp == true",
        tabsetPanel(
          id = "inpsybapp",
          type = "tabs",
          tabPanel("Raw Data", tableOutput("readappsybr")),
          tabPanel("Run Information", tableOutput("sybrruninfoapp")),
          tabPanel("Raw Data", tableOutput("sybrdataapp")),
          tabPanel("Ct Plate", dataTableOutput("ctplatesybrapp")),
          tabPanel("Sample Plate", dataTableOutput("sampleplatesybrapp")),
          tabPanel("Standard Curve", tableOutput("stdcurvesybrapp"), plotOutput("stdsybrapp")),
          tabPanel("Analysis", downloadButton("downansybrapp", "Download CSV"), dataTableOutput("analysissybrapp")),
          tabPanel("ID_Well", downloadButton("downIDWELLsybrapp", "Download CSV"), tableOutput("IDWELLsybrapp")),
          tabPanel("ID_Result", downloadButton("downIDRESsybrapp", "Download CSV"), tableOutput("IDRESsybrapp"))
        )
      ),
      ########### MELTING CURVES ###########
      conditionalPanel(
        condition = "input.meltsybr == true",
        tabsetPanel(
          id = "syb",
          type = "tabs",
          tabPanel("ID_well", tableOutput("well")),
          tabPanel("TM Plots", uiOutput("tmplots")),
          tabPanel("TM Table", downloadButton("downloadtable", "Download CSV"),
                   tableOutput("tmtable")
          )
        )
      )
    )
  )
)
