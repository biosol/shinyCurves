library(shiny)
library(data.table)
library(qpcR)
library(ggplot2)
library(googledrive)
library(reshape2)
library(dplyr)
library(cowplot)
library(plyr)
library(matrixStats)
library(rlist)
library(DT)
library(readxl)
library(ggpmisc)
library(gridExtra)
library(outviz)
library(tidyverse)
library(janitor)
source("/home/sonia/Documentos/Covid19/shiny-app/ui.R")
source("/home/sonia/Documentos/Covid19/shiny-app/server.R")

setwd("/home/sonia/Documentos/Covid19/shiny-app/qPCR-shiny-1653ffa849b74fd4312faca5729450bf71602a9b/")
#Call app
app <- shinyApp(ui = ui, server = server)
#Run app
runApp()

# Example of shiny app paper  
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229330
