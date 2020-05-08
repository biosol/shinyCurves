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
source("/home/sonia/Documentos/Covid19/shiny-app/ui.R")
source("/home/sonia/Documentos/Covid19/shiny-app/server.R")

setwd("/home/sonia/Documentos/Covid19/shiny-app/")
#Call app
app <- shinyApp(ui = ui, server = server)
#Run app
runApp()
