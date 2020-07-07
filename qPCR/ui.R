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
      checkboxInput("taqexcel", "1a) Input: Taqman - BioRad"),
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
      checkboxInput("app", "1b) Input: Taqman - Applied Quant Studio"),
      conditionalPanel(
        condition = "input.app == true",
        fileInput("appl",
                  label = "Upload Applied Quant Studio xlsx",
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
      checkboxInput("taqman", "2) Analysis: Taqman Ct Curves"),
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
      checkboxInput("inpsybr", "3a) Input: SYBR"),
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
      checkboxInput("meltsybr", "3b) Analysis: SYBR Melting Curves"),
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
      ###### 3) Main Panel for SYBR ########
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
