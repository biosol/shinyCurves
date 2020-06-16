################################# SHINY APP FOR QPCR ANALYSIS ############################################
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
                  multiple = TRUE)
      ),
      checkboxInput("app", "1b) Input: Taqman - Applied Quant Studio"),
      conditionalPanel(
        condition = "input.app == true",
        fileInput("appl",
                  label = "Upload Applied Quant Studio xlsx",
                  multiple = TRUE)
      ),
      checkboxInput("taqman", "2) Analysis: Taqman Cq Curves"),
      conditionalPanel(
        condition = "input.taqman == true",
        fileInput("taqmancsv",
                  label = "Upload CSVs here",
                  multiple = TRUE),
        numericInput("ct", "Enter Ct value", value = 40),
        textInput("endoC", "Enter endogenous control", value = "RNAseP"),
        #fileInput("dataendoC", "Enter qPCR results for endoC"),
        fileInput("taqwellid", label = "Upload ID_well.csv here"),
        fileInput("taqidres", label = "Upload ID_results.csv here")
      ),
      checkboxInput("inpsybr", "3a) Input: SYBR"),
      conditionalPanel(
        condition = "input.inpsybr == true",
        fileInput("sybr",
                  label = "Upload Excel or CSVs here",
                  multiple = TRUE),
      ),
      checkboxInput("meltsybr", "3b) Analysis: SYBR Melting Curves"),
      conditionalPanel(
        condition = "input.meltsybr == true",
        fileInput("sybrcsv",
                  label="Upload CSV here",
                  multiple = TRUE),
        radioButtons("isderiv", "Isderiv?",
                     choices = c(True = TRUE,
                                 False = FALSE),
                     selected = TRUE),
        fileInput("wellid",
                  label="Upload ID_well.csv here")
      )
    ),
      mainPanel(
      ###### 1a) Main Panel for Biorad Input ########
      conditionalPanel(
        condition = "input.taqexcel == true",
        tabsetPanel(
          id = "taqex",
          type = "tabs",
          tabPanel("Empty Plate", dataTableOutput("plate")),
          tabPanel("Raw Data", tableOutput("inputdf")),
          tabPanel("Run Information", tableOutput("BRruninfo")),
          #tabPanel("New Data", tableOutput("newdatabiorad")),
          tabPanel("Cq Plate", dataTableOutput("cqplate")),
          tabPanel("Sample Plate", dataTableOutput("sampleplate")),
          tabPanel("Plate Setup MultiChanel", dataTableOutput("setupmultic")),
          tabPanel("Sample Order Check", dataTableOutput("samplecheck")),
          tabPanel("Standard Curve", fluidRow(dataTableOutput("stdcurve"), plotOutput("std"))),
          tabPanel("Analysis Samples", dataTableOutput("analysis")),
          tabPanel("Interpretation", dataTableOutput("interpret")),
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
          tabPanel("Conversion",
                   tabsetPanel(
                     tabPanel("Tranposed", tableOutput("trans")),
                     tabPanel(title=uiOutput("gen1app"), downloadButton("downtransgen1", "Download CSV"), tableOutput("transgen1")),
                     tabPanel(title=uiOutput("gen2app"), downloadButton("downtransgen2", "Download CSV"), tableOutput("transgen2")),
                     tabPanel(title=uiOutput("gen3app"), downloadButton("downtransgen3", "Download CSV"), tableOutput("transgen3"))
                   )
          ),
          tabPanel("Run Info", tableOutput("appliedruninfo")),
          tabPanel("Applied Results", tableOutput("appliedres")),
          tabPanel("Cq Plate", dataTableOutput("cqplateapp")),
          tabPanel("Sample Plate", dataTableOutput("sampleplateapp")),
          tabPanel("Check Samples", dataTableOutput("samplecheckapp")),
          tabPanel("Plate Setup Multichanel", dataTableOutput("setupmulticapp")),
          tabPanel("Std Curves", dataTableOutput("stdcurveapp"), plotOutput("stdapp")),
          tabPanel("Analysis Samples", dataTableOutput("analysisapp")),
          tabPanel("Interpretation", dataTableOutput("interpretapp")),
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
          tabPanel("C", tableOutput("c")),
          tabPanel("Info", tableOutput("info")),
          tabPanel("General Plots", downloadButton("downgen", "Download PDF"), plotOutput("genplots", height = "1000px")),
          tabPanel("Indet Plots", 
                   tabsetPanel(
                     tabPanel(title=uiOutput("gen1"),fluidRow(downloadButton("downlgen1", "Download PDF"), plotOutput("indetgen1", height = "500px"))),
                     tabPanel(title=uiOutput("gen2"), fluidRow(downloadButton("downlgen2", "Download PDF"), plotOutput("indetgen2", height = "500px"))),
                     tabPanel(title=uiOutput("gen3"), fluidRow(downloadButton("downlgen3", "Download PDF"), plotOutput("indetgen3", height = "500px")))
                   )
          )
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
          tabPanel("New Data", tableOutput("newdatasybr")),
          tabPanel("Cq Plate", dataTableOutput("cqplatesybr")),
          tabPanel("Sample Plate", dataTableOutput("sampleplatesybr")),
          tabPanel("Check Sample Order", dataTableOutput("samplechecksybr")),
          tabPanel("Plate Setup Multichanel", dataTableOutput("setupmulticsybr")),
          tabPanel("Standard Curve", dataTableOutput("stdcurvesybr"),plotOutput("stdsybr")),
          tabPanel("Analysis Samples", dataTableOutput("analysissybr")),
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
          tabPanel("Fluos_Gene", tableOutput("fluosgene")),
          tabPanel("Fluo_Temp", tableOutput("fluotemp")),
          tabPanel("TM Plots", 
                   tabsetPanel(
                       tabPanel("N", fluidRow(downloadButton("downsybrN", "Download PDF"), plotOutput("sybrN", height = "1000px"))),
                       tabPanel("Rdrp",fluidRow(downloadButton("downsybrRdrp", "Download PDF"), plotOutput("sybrRdrp", height = "1000px"))),
                       tabPanel("Rpp30", fluidRow(downloadButton("downsybrRpp30", "Download PDF"), plotOutput("sybrRpp30", height = "1000px"))),
                       tabPanel("S", fluidRow(downloadButton("downsybrS", "Download PDF"), plotOutput("sybrS", height = "1000px")))
                   ),
          ),
          tabPanel("TM Table", downloadButton("downloadtable", "Download CSV"), actionButton("send2drive", "Send to Google Drive", icon = icon("google-drive", lib="font-awesome")),
                   tableOutput("tmtable")
          )
        )
      )
    )
  )
)
