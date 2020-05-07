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
                  label = "Upload CSV here",
                  multiple = FALSE),
        numericInput("ct", "Enter PCR number of cycles", value = 40),
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
                  multiple = FALSE),
        radioButtons("isderiv", "Isderiv?",
                     choices = c(True = TRUE,
                                 False = FALSE),
                     selected = TRUE),
        fileInput("wellid",
                  label="Upload ID_well.csv here"),
      )
    ),
    mainPanel(
      ###### Main panel appears only of Taqman is checked ########
      conditionalPanel(
        condition = "input.taqman == true",
        tabsetPanel(
          id = "taq",
          type = "tabs",
          tabPanel("Data", tableOutput("incycles")),
          tabPanel("ID_Well", tableOutput("taqmanwell")),
          tabPanel("ID_Result", tableOutput("taqmanidres")),
          tabPanel("Info", tableOutput("info")),
          tabPanel("Indet Plots", plotOutput("indetplots"),
                   downloadButton("downloadIndet", "Download PNG"))
        )
      ),
      ###### Main panel appears only of SYBR is checked ########
      conditionalPanel(
        condition = "input.sybr == true",
        tabsetPanel(
          id = "syb",
          type = "tabs",
          tabPanel("Gene", tableOutput("tables")),
          tabPanel("Plate Map", tableOutput("well")),
          tabPanel("TM Plots", plotOutput("TM"),
                   downloadButton("downloadplot", "Download PNG")),
          tabPanel("TM Table", tableOutput("tmtable"),
                   downloadButton("downloadtable", "Download"),
                   actionButton("send2drive", "Send to Google Drive", icon = icon("google-drive", lib="font-awesome"))
          )
        )),
    )
  )
)
