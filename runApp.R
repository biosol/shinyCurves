################################# SHINY APP (COMPLETE) FOR QPCR ANALYSIS ############################################
########################### Sonia Olaechea-Lázaro (UPV/EHU, May 2020) ############################################  

########### runApp.R #################
## 0) Move to your working directory
setwd("/home/sonia/Documentos/Covid19/shiny-app/")

## 1) First specify the packages of interest
packages = c("shiny", "data.table", "qpcR","reshape2","cowplot","rlist","DT", 
             "readxl","plyr", "ggplot2","dplyr", "gridExtra","matrixStats",
             "ggpmisc","outviz", "tidyverse", "tidyr", "janitor")

## 2) Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## 3) Load the ui and server functions (here you have to change the path for your computer)
source("/home/sonia/Documentos/Covid19/shiny-app/ui.R")
source("/home/sonia/Documentos/Covid19/shiny-app/server.R")

## 4) Call app (store app in variable)
app <- shinyApp(ui = ui, server = server)

## 5) Run app
runApp()

# Example of shiny app paper  
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229330
# SYBR
'tabsetPanel(
  tabPanel("N", fluidRow(downloadButton("downsybrN", "Download PDF"), plotOutput("sybrN", height = "1000px"))),
  tabPanel("Rdrp",fluidRow(downloadButton("downsybrRdrp", "Download PDF"), plotOutput("sybrRdrp", height = "1000px"))),
  tabPanel("Rpp30", fluidRow(downloadButton("downsybrRpp30", "Download PDF"), plotOutput("sybrRpp30", height = "1000px"))),
  tabPanel("S", fluidRow(downloadButton("downsybrS", "Download PDF"), plotOutput("sybrS", height = "1000px")))
),'

#Taqman
'tabsetPanel(
  tabPanel(title=uiOutput("gen1"), fluidRow(uiOutput("downlgen1"), uiOutput("indetgen1"))), #downloadButton("downlgen1", "Download PDF"), plotOutput("indetgen1", height = "500px") 
  tabPanel(title=uiOutput("gen2"), fluidRow(uiOutput("downlgen2"), uiOutput("indetgen2"))), #downloadButton("downlgen2", "Download PDF"), plotOutput("indetgen2", height = "500px")
  tabPanel(title=uiOutput("gen3"), fluidRow(uiOutput("downlgen3"), uiOutput("indetgen3")))  #downloadButton("downlgen3", "Download PDF"), plotOutput("indetgen3", height = "500px")
)'