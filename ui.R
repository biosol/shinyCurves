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
                  label = "Upload BioRad CFX results file")
        ),
       checkboxInput("app",  "1b) Input: Taqman - Applied Quant Studio"),
       conditionalPanel(
        condition = "input.app == true",
        fileInput("appl",
                    label = "Upload Applied Quant Studio xlsx")
        ),
       checkboxInput("taqman", "2) Analysis: Taqman"),
       conditionalPanel(
        condition = "input.taqman == true",
        fileInput("taqmancsv",
                  label = "Upload CSVs here",
                  multiple = TRUE),
        numericInput("ct", "Enter Ct value", value = 40),
        textInput("endoC", "Enter endogenous control", value = "RNAseP"),
        fileInput("dataendoC", "Enter qPCR results for endoC"),
        fileInput("taqwellid", label = "Upload ID_well.csv here"),
        fileInput("taqidres", label = "Upload ID_results.csv here")
      ),
      ###### Sidebar panel appears only of SYBR is checked ########
      checkboxInput("sybr", "3) Input + Analysis: SYBR Melting Curves"),
      conditionalPanel(
        condition = "input.sybr == true",
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
      ###### Main panel appears only of Taqman Excel is checked ########
      conditionalPanel(
        condition = "input.taqexcel == true",
        tabsetPanel(
          id = "taqex",
          type = "tabs",
          tabPanel("Empty Plate", dataTableOutput("plate")),
          tabPanel("Input DF", tableOutput("inputdf")),
          tabPanel("Run Information", tableOutput("BRruninfo")),
          tabPanel("Cq Plate", dataTableOutput("cqplate")),
          tabPanel("Sample Plate", dataTableOutput("sampleplate")),
          tabPanel("Plate Setup MultiChanel", dataTableOutput("setupmultic")),
          tabPanel("Sample Order Check", dataTableOutput("samplecheck")),
          tabPanel("Standard Curve", 
                   fluidRow(
                     DT::dataTableOutput("stdcurve"),
                     plotOutput("std"))
          ),
          tabPanel("Analysis Samples", dataTableOutput("analysis")),
          tabPanel("ID_Well", tableOutput("IDWELL"))
        )
      ),
      ########### Taqman Amplification Curves ###############
      conditionalPanel(
        condition = "input.taqman == true",
        tabsetPanel(
          id = "taq",
          type = "tabs",
          tabPanel("ID_Well", tableOutput("taqmanwell")),
          tabPanel("ID_Result", tableOutput("taqmanidres")),
          tabPanel("C", tableOutput("c")),
          tabPanel("Info", tableOutput("info")),
          tabPanel("General Plots",
                   downloadButton("downgen", "Download PDF"),
                   plotOutput("genplots", height = "1000px")
                   ),
          tabPanel("Indet Plots", 
                   tabsetPanel(
                     tabPanel("N1",
                              fluidRow(downloadButton("downln1", "Download PDF"),
                                       plotOutput("n1", height = "500px"))
                     ),
                     tabPanel("N2",
                              fluidRow(downloadButton("downln2", "Download PDF"),
                                       plotOutput("n2", height = "500px"))
                     ),
                     tabPanel("RNAseP",
                              fluidRow(downloadButton("downrnasep", "Download PDF"),
                                       plotOutput("rnasep", height = "500px"))
                     )
                   )
          )
        )
      ),
      conditionalPanel(
        condition = "input.app == true",
        tabsetPanel(
          id = "blabla",
          type = "tabs",
          tabPanel("Raw Data", tableOutput("prueba")),
          tabPanel("Transposed", tableOutput("trans")),
          tabPanel("N1",
                   downloadButton("downtransN1", "Download CSV"),
                   tableOutput("transN1"),),
          tabPanel("N2",
                   downloadButton("downtransN2", "Download CSV"),
                   tableOutput("transN2")),
          tabPanel("RNAseP",
                   downloadButton("downtransRNAseP", "Download CSV"),
                   tableOutput("transRNAsep"))
        )
      ),
      ###### Main panel appears only of SYBR is checked ########
      conditionalPanel(
        condition = "input.sybr == true",
        tabsetPanel(
          id = "syb",
          type = "tabs",
          tabPanel("ID_well", tableOutput("well")),
          tabPanel("Fluos_Gene", tableOutput("fluosgene")),
          tabPanel("Fluo_Temp", tableOutput("fluotemp")),
          tabPanel("TM Plots", 
                   tabsetPanel(
                       tabPanel("N",
                              fluidRow(downloadButton("downsybrN", "Download PDF"),
                                       plotOutput("sybrN", height = "1000px"))
                       ),
                       tabPanel("Rdrp",
                              fluidRow(downloadButton("downsybrRdrp", "Download PDF"),
                                       plotOutput("sybrRdrp", height = "1000px"))
                      ),
                      tabPanel("Rpp30",
                              fluidRow(downloadButton("downsybrRpp30", "Download PDF"),
                                       plotOutput("sybrRpp30", height = "1000px"))
                      ),
                      tabPanel("S",
                             fluidRow(downloadButton("downsybrS", "Download PDF"),
                                      plotOutput("sybrS", height = "1000px"))
                      )
                  )), 
          tabPanel("TM Table",
                   downloadButton("downloadtable", "Download"),
                   actionButton("send2drive", "Send to Google Drive", icon = icon("google-drive", lib="font-awesome")),
                   tableOutput("tmtable")
          )
        ))
    )
  )
)
