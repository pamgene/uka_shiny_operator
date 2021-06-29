library(shiny)
library(tercen)
library(dplyr)
library(tidyr)

############################################
#### This part should not be modified
getCtx <- function(session) {
  # retreive url query parameters provided by tercen
  query <- parseQueryString(session$clientData$url_search)
  token <- query[["token"]]
  taskId <- query[["taskId"]]
  
  # create a Tercen context object using the token
  ctx <- tercenCtx(taskId = taskId, authToken = token)
  return(ctx)
}
####
############################################

shinyServer(function(input, output, session) {
  
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
  
  # TODO database connection
  #DB = UpstreamApp::UpstreamDatabase
  #nid = showNotification("Press Start to start the analysis.", duration = NULL, type = "message", closeButton = FALSE)
  #updateSliderInput(session, "seqHom", min = min(DB$PepProtein_SeqHomology))
  
  observe({
    ctx  <- getCtx(session)
    data <- dataInput()
    prop <- propertiesInput()
    
    if (is.null(getProperties)) return()
    updateSelectInput(session, "kinasefamily", selected = prop$Kinase_family)
    
    if (is.null(data)) return()
    
    if (length(ctx$cnames) > 0) {
      # poptions = data$colorColumnNames
      # names(poptions) = data$arrayLabels
      poptions <- unlist(ctx$cnames)
      updateSelectInput(session, "pairingFactor", choices = poptions)
    }
    output$status = renderText({
      grp = getGroupingLabel(data)
      if (input$start == 0){
        if(prop$Lock_kinase_family == "Yes"){
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
      
      if(!data$hasUniqueDataMapping) stop("Mapping error: not all input data is unique")
      
      if(data$hasMissingCells) stop("Missing values are not allowed.")
      
      if(data$getMaxNPerCell() > 1) stop("More than one value per cell in the cross-tab view is not allowed.")
      
      
      
      df = data$getData()
      
      if(data$hasZeroScaleRows){
        zIdx = data$getZeroScaleRows()
        df = df %>% filter (!(rowSeq %in% zIdx))
        msg = paste("Warning:", length(zIdx), "peptides with zero scale have been removed from the data")
        showNotification(ui = msg, duration = NULL, type = "warning")
      } 
      
      if(input$kinasefamily == "PTK"){
        DB = DB %>% filter(PepProtein_Residue == "Y")
      } else if(input$kinasefamily == "STK") {
        DB = DB %>% filter(PepProtein_Residue == "S" |  PepProtein_Residue == "T")
      } else {
        stop("Unknown value for kinase family")
      }
      
      DB = DB %>% filter(PepProtein_SeqHomology >= input$seqHom) 
      DB = DB %>% filter(Kinase_PKinase_PredictorVersion2Score >= input$minPScore | Database == "iviv")
      DB = DB %>% filter(Kinase_Rank <= input$scan[2])
      
      nCores = detectCores()
      msg = paste("Please wait ... running analysis. Using", nCores, "cores.")
      showNotification(ui = msg, id = nid, type = "message", closeButton = FALSE, duration = NULL)
      
      if(input$seed){
        set.seed(42)
      }
      
      if(!is.null(grp)){
        df$grp = as.factor(df[[grp]])
        aResult = pgScanAnalysis2g(df,  dbFrame = DB,
                                   scanRank = input$scan[1]:input$scan[2],
                                   nPermutations = input$nperms,
                                   dbWeights = c(iviv = input$wIviv,
                                                 PhosphoNET = input$wPhosphoNET
                                   ))
      } else {
        aResult = pgScanAnalysis0(df,  dbFrame = DB,
                                  scanRank = input$scan[1]:input$scan[2],
                                  nPermutations = input$nperms,
                                  dbWeights = c(iviv = input$wIviv,
                                                PhosphoNET = input$wPhosphoNET
                                  ))
      }
      showNotification(ui = "Done", id = nid, type = "message", closeButton = FALSE)
      aFull = ldply(aResult, .fun = function(.)return(data.frame(.$aResult, mxRank = .$mxRank) ))
      
      
      settings = data.frame(setting = c("Kinase family", "ScanRank Min", "ScanRank Max", "Number of Permutations", "In Vitro In Vitro weight", "PhosphoNET weight", "Min PhosphoNet score", "Min Sequence Homology"),
                            value   = c(input$kinasefamily, input$scan[1] , input$scan[2], input$nperms, input$wIviv, input$wPhosphoNET, input$minPScore, input$seqHom) )
      
      spath = file.path(getFolder(), "runData.RData")
      save(file = spath, df, aResult, aFull, settings)
      dpath = file.path(getFolder(), "runDb.RData")
      save(file = dpath, DB)
      out = data.frame(rowSeq = 1, colSeq = 1, dummy = NaN)
      meta = data.frame(labelDescription = c("rowSeq", "colSeq", "dummy"), groupingType = c("rowSeq", "colSeq", "QuantitationType"))
      result = AnnotatedData$new(data = out, metadata = meta)
      context$setResult(result)
      return("Done")
    })
  })
})

getValues <- function(session){
  ctx <- getCtx(session)
  values <- list()
  
  values$data <- ctx %>% select(.y, .ri, .ci) %>%
    mutate(color = ctx$select(ctx$colors) %>% pull) %>%
    group_by(.ci, .ri, color) %>%
    summarise(.y = mean(.y)) # take the mean of multiple values per cell
  
  values$colorLabels <- colnames(ctx$select(ctx$colors))
  
  return(values)
}

getProperties <- function(session) {
  ctx <- getCtx(session)
  
  list(Kinase_family      = ifelse(is.null(ctx$op.value("Kinase_family")), "PTK", ctx$op.value("Kinase_family")), 
       Lock_kinase_family = ifelse(is.null(ctx$op.value("Lock_kinase_family")), "Yes", ctx$op.value("Lock_kinase_family")))
}
