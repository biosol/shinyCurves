################################# SHINY APP FOR QPCR ANALYSIS ############################################
########################### Sonia Olaechea-Lázaro (UPV/EHU, May 2020) ############################################
ui <- fluidPage(
  titlePanel(strong("qPCR Analysis")),
  sidebarLayout(
    sidebarPanel(
      checkboxInput("taqman", "Taqman"),
      ###### Sidebar panel appears only of Taqman is checked ########
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
      checkboxInput("sybr", "SYBR"),
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
      ###### Main panel appears only of Taqman is checked ########
      conditionalPanel(
        condition = "input.taqman == true",
        tabsetPanel(
          id = "taq",
          type = "tabs",
          tabPanel("ID_Well", tableOutput("taqmanwell")),
          tabPanel("ID_Result", tableOutput("taqmanidres")),
          tabPanel("C", tableOutput("c")),
          tabPanel("Info", tableOutput("info")),
          tabPanel("General Plots", plotOutput("genplots", height = "1000px"),
                   downloadButton("downgen", "Download PDF")),
          tabPanel("Indet Plots", 
                   tabsetPanel(
                     tabPanel("N1",
                              fluidRow(plotOutput("n1", height = "1000px"),
                                       downloadButton("downln1", "Download PDF"))
                     ),
                     tabPanel("N2",
                              fluidRow(plotOutput("n2", height = "1000px"),
                                       downloadButton("downln2", "Download PDF"))
                     ),
                     tabPanel("RNAseP",
                              fluidRow(plotOutput("rnasep", height = "1000px"),
                                       downloadButton("downrnasep", "Download PDF"))
                     )
                   )
           )
        )
      ),
      ###### Main panel appears only of SYBR is checked ########
      conditionalPanel(
        condition = "input.sybr == true",
        tabsetPanel(
          id = "syb",
          type = "tabs",
          #tabPanel("Gene", tableOutput("tables")),
          tabPanel("ID_well", tableOutput("well")),
          tabPanel("Fluos_Gene", tableOutput("fluosgene")),
          tabPanel("Fluo_Temp", tableOutput("fluotemp")),
          tabPanel("TM Plots", 
                   tabsetPanel(
                     tabPanel("N",
                              fluidRow(plotOutput("sybrN", height = "1000px"),
                                       downloadButton("downsybrN", "Download PDF"))
                       
                     ),
                     tabPanel("Rdrp",
                              fluidRow(plotOutput("sybrRdrp", height = "1000px"),
                                       downloadButton("downsybrRdrp", "Download PDF"))
                     ),
                     tabPanel("Rpp30",
                              fluidRow(plotOutput("sybrRpp30", height = "1000px"),
                                       downloadButton("downsybrRpp30", "Download PDF"))
                    ),
                    tabPanel("S",
                             fluidRow(plotOutput("sybrS", height = "1000px"),
                                      downloadButton("downsybrS", "Download PDF"))
                    )
                  )), 
          tabPanel("TM Table", tableOutput("tmtable"),
                   downloadButton("downloadtable", "Download"),
                   actionButton("send2drive", "Send to Google Drive", icon = icon("google-drive", lib="font-awesome"))
          )
        ))
    )
  )
)
