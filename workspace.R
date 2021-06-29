library(shiny)
library(tercen)
library(dplyr)
library(tidyr)
library(shinyjs)
library(parallel)
source("upstream_analysis.R")

############################################
#### This part should not be included in ui.R and server.R scripts
getCtx <- function(session) {
  # Set appropriate options
  #options("tercen.serviceUri"="http://tercen:5400/api/v1/")
  #options("tercen.workflowId"= "4133245f38c1411c543ef25ea3020c41")
  #options("tercen.stepId"= "2b6d9fbf-25e4-4302-94eb-b9562a066aa5")
  #options("tercen.username"= "admin")
  #options("tercen.password"= "admin")
  
  options("tercen.serviceUri"="https://tercen.com/")
  options("tercen.username"= "ginberg")
  options("tercen.password"= 'w:~9u203-@,uL[zi5q{!N_$uN+_"R:y6FSmZ&6`mujgTE/]=')
  options("tercen.workflowId"= "4629c134b09c53160ea461e75fe42d12")
  options("tercen.stepId"= "2d17de25-3d52-4964-8ac4-1f514af62a93")
  ctx <- tercenCtx()
  return(ctx)
}
####
############################################

ui <- shinyUI(fluidPage(
  shinyjs::useShinyjs(),
  tags$script(HTML('setInterval(function(){ $("#hiddenButton").click(); }, 1000*30);')),
  tags$footer(shinyjs::hidden(actionButton(inputId = "hiddenButton", label = "hidden"))),
  uiOutput("body"),
  title = "Kinase Analysis"
))

server <- shinyServer(function(input, output, session) {
  
  dataInput <- reactive({
    getValues(session)
  })
  
  propertiesInput <- reactive({
    getProperties(session)
  })
  
  output$body <- renderUI({
    sidebarLayout(
      sidebarPanel(
        tags$div(HTML("<strong><font color = #6895d1>Upstream Kinase Analysis</font></strong>")),
        tags$hr(),
        actionButton("start", "Start"),
        tags$hr(),
        textOutput("status")
        
      ),
      mainPanel(tabsetPanel(
                  tabPanel("Basic Settings",
                           tags$hr(),
                           selectInput("kinasefamily", label = "Kinase family", choices = c("PTK", "STK")),
                           tags$hr(),
                           checkboxInput("usePairing", label = "Use a Pairing Factor"),
                           conditionalPanel("input.usePairing == true",
                                            selectInput("pairingFactor",label = "Select the Pairing Factor", choices = "" )
                           )
                  ),
                  
                  tabPanel("Advanced Settings",
                           sliderInput("scan", "Scan Rank From-To", min = 1, max = 20, value = c(4,12),round = TRUE),
                           sliderInput("nperms", "Number of permutations", min = 100, max = 1000, value = 500, step = 100, round = TRUE),
                           tags$hr(),
                           helpText("Set weights for database types:"),
                           sliderInput("wIviv", "In Vitro / In Vivo", min = 0, max = 1, value = 1),
                           sliderInput("wPhosphoNET", "In Silico (PhosphoNET)", min = 0, max = 1, value = 1),
                           tags$hr(),
                           helpText("Set minimal sequence homology required to link a peptide to a phosphosite"),
                           sliderInput("seqHom", "Minimal sequence homology", min =0, max = 1, value =0.9),
                           tags$hr(),
                           numericInput("minPScore", "Minimal PhosphoNET prediction score", value = 300),
                           tags$hr(),
                           checkboxInput("seed", "Set seed for random permutations"))
                )
      )
    )
  })
  
  getGroupingLabel = function(values) {
    if (length(values$colorLabels) > 1) {
      stop("Need at most one Color to define the grouping")
    }
    else if(length(values$colorLabels) == 1) {
      grouping_label <- values$colorLabels
      group_factor   <- as.factor(values$data$color)
      if (length(levels(group_factor)) != 2){
        stop(paste("Wrong number of groups, found:", levels(group_factor)))
      }
      
      df           <- subset(values$data, .ri = 0)
      group_factor <- as.factor(df$color)
      if (min(summary(group_factor)) < 2){
        stop("Need at least 2 observations per group.")
      }
      grouping_label
    }
    else {
      if (max(values$data$.ci) > 0) {
        stop(paste("Can not run this data without a grouping factor.", values$colorLabels))
      } else return(NULL)
    }
  }
  
  # read database file
  DB  <- readRDS("db2.rds")
  nid <- showNotification("Press Start to start the analysis.", duration = NULL, type = "message", closeButton = FALSE)
  updateSliderInput(session, "seqHom", min = min(DB$PepProtein_SeqHomology))
  
  observe({
    ctx  <- getCtx(session)
    data_input <- dataInput()
    properties <- propertiesInput()
    
    if (is.null(properties)) return()
    updateSelectInput(session, "kinasefamily", selected = properties$Kinase_family)
    
    if (is.null(data_input)) return()
    
    if (length(ctx$cnames) > 0) {
      # poptions = data$colorColumnNames
      # names(poptions) = data$arrayLabels
      poptions <- unlist(ctx$cnames)
      updateSelectInput(session, "pairingFactor", choices = poptions)
    }
    output$status = renderText({
      grp = getGroupingLabel(data_input)
      if (input$start == 0){
        if(properties$Lock_kinase_family == "Yes"){
          shinyjs::disable("kinasefamily")
        }
        if (!is.null(grp)){
          return(paste("Grouping factor:", grp))
        } else {
          return("Grouping factor: none")
        }
      }
      
      shinyjs::disable("kinasefamily")
      shinyjs::disable("scan")
      shinyjs::disable("nperms")
      shinyjs::disable("wIviv")
      shinyjs::disable("wPhosphoNET")
      shinyjs::disable("seqHom")
      shinyjs::disable("minPScore")
      shinyjs::disable("usePairing")
      shinyjs::disable("pairingFactor")
      
      # TODO checks
      # if(!data$hasUniqueDataMapping) stop("Mapping error: not all input data is unique")
      # if(data$hasMissingCells) stop("Missing values are not allowed.")
      # if(data$getMaxNPerCell() > 1) stop("More than one value per cell in the cross-tab view is not allowed.")
      
      df <- data_input$data
      
      if (!is.null(data_input$hasZeroScaleRows) && data_input$hasZeroScaleRows) {
        zIdx <- data_input$getZeroScaleRows()
        df   <- df %>% filter (!(rowSeq %in% zIdx))
        msg  <- paste("Warning:", length(zIdx), "peptides with zero scale have been removed from the data")
        showNotification(ui = msg, duration = NULL, type = "warning")
      } 
      
      if (input$kinasefamily == "PTK") {
        DB <- DB %>% filter(PepProtein_Residue == "Y")
      } else if(input$kinasefamily == "STK") {
        DB <- DB %>% filter(PepProtein_Residue == "S" |  PepProtein_Residue == "T")
      } else {
        stop("Unknown value for kinase family")
      }
      
      DB <- DB %>% 
        filter(PepProtein_SeqHomology >= input$seqHom)  %>%
        filter(Kinase_PKinase_PredictorVersion2Score >= input$minPScore | Database == "iviv") %>%
        filter(Kinase_Rank <= input$scan[2])
      
      nCores <- detectCores()
      msg    <- paste("Please wait ... running analysis. Using", nCores, "cores.")
      showNotification(ui = msg, id = nid, type = "message", closeButton = FALSE, duration = NULL)
      
      if (input$seed) {
        set.seed(42)
      }
      
      if (!is.null(grp)) {
        df$grp <- as.factor(df$color)
        result <- pgScanAnalysis2g(df, dbFrame = DB,
                                   scanRank = input$scan[1]:input$scan[2],
                                   nPermutations = input$nperms,
                                   dbWeights = c(iviv = input$wIviv,
                                                 PhosphoNET = input$wPhosphoNET
                                   ))
      } else {
        result <- pgScanAnalysis0(df, dbFrame = DB,
                                  scanRank = input$scan[1]:input$scan[2],
                                  nPermutations = input$nperms,
                                  dbWeights = c(iviv = input$wIviv,
                                                PhosphoNET = input$wPhosphoNET
                                  ))
      }
      showNotification(ui = "Done", id = nid, type = "message", closeButton = FALSE)
      full_result <- ldply(result, .fun = function(.) return(data.frame(.$result, mxRank = .$mxRank)))
      
      settings = data.frame(setting = c("Kinase family", "ScanRank Min", "ScanRank Max", "Number of Permutations", "In Vitro In Vitro weight", "PhosphoNET weight", "Min PhosphoNet score", "Min Sequence Homology"),
                            value   = c(input$kinasefamily, input$scan[1] , input$scan[2], input$nperms, input$wIviv, input$wPhosphoNET, input$minPScore, input$seqHom) )
      
      # TODO save objects in tercen context!
      
      # spath = file.path(getFolder(), "runData.RData")
      # save(file = spath, df, result, full_result, settings)
      # dpath = file.path(getFolder(), "runDb.RData")
      # save(file = dpath, DB)
      # out = data.frame(rowSeq = 1, colSeq = 1, dummy = NaN)
      # meta = data.frame(labelDescription = c("rowSeq", "colSeq", "dummy"), groupingType = c("rowSeq", "colSeq", "QuantitationType"))
      # result = AnnotatedData$new(data = out, metadata = meta)
      # context$setResult(result)
      return("Done")
    })
  })
})

getValues <- function(session){
  ctx <- getCtx(session)
  values <- list()
  
  # TODO check if ID in rnames
  data     <- ctx %>% 
    select(.y, .ri, .ci) %>%
    mutate(color = ctx$select(ctx$colors) %>% pull)
  row_data <- ctx %>% 
    rselect("ID") %>% 
    mutate(ID = as.factor(ID)) %>%
    mutate(.ri = seq(0, length(unique(data$.ri))-1))
  
  values$data <- data %>%
    left_join(., row_data)
  
  values$colorLabels <- colnames(ctx$select(ctx$colors))
  
  return(values)
}

getProperties <- function(session) {
  ctx <- getCtx(session)
  
  list(Kinase_family      = ifelse(is.null(ctx$op.value("Kinase_family")), "PTK", ctx$op.value("Kinase_family")), 
       Lock_kinase_family = ifelse(is.null(ctx$op.value("Lock_kinase_family")), "Yes", ctx$op.value("Lock_kinase_family")))
}

runApp(shinyApp(ui, server))  