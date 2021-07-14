library(shiny)
library(tercen)
library(plyr)
library(dplyr)
library(tidyr)
library(shinyjs)
library(parallel)
library(reshape2)
library(pgCheckInput)
library(colourpicker)
library(rmarkdown)
library(DT)
source("upstream_analysis.R")
source("upstream_output.R")

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

server <- shinyServer(function(input, output, session) {
  
  results <- reactiveValues()
  viewer  <- reactiveValues(view = "run")
  # read database file
  DB      <- readRDS("data/db.rds")
  
  dataInput <- reactive({
    getValues(session)
  })
  
  propertiesInput <- reactive({
    getProperties(session)
  })
  
  getComputedResults <- reactive({
    getCtxResults(session)
  })
  
  modeInput <- reactive({ 
    getMode(session) 
  })
  
  getView <- reactive({
    result          <- "run"
    ctx             <- getCtx(session)
    computedResults <- getComputedResults()
    switchToRun     <- input$switchToRun
    if (!is.null(computedResults) && inputDataEqual(ctx, computedResults) &&
        (is.null(switchToRun) || !is.null(switchToRun) && switchToRun == 0)) {
      result <- "showResult"
    }
    result
  })
  
  output$body <- renderUI({
    if (isRunView(getView())) {
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
    } else {
      computedResults     <- getComputedResults()
      results$DB          <- computedResults$DB
      results$full_result <- computedResults$full_result
      results$df          <- computedResults$df
      
      results$kinase2uniprot <- results$DB %>%
        group_by(Kinase_UniprotID) %>%
        dplyr::summarise(Kinase_Name = Kinase_Name[1])
      
      sidebarLayout(
        sidebarPanel(
          tags$div(HTML("<strong><font color = #6895d1>Upstream Kinase Analysis</font></strong>")),
          tags$hr(),
          textOutput("grpText"),
          tags$hr(),
          sliderInput("minsetsize", "Include results based on peptide sets with minimal size of:",
                      min = 1,
                      max = 5,
                      step = 1,
                      value = 3),
          tags$hr(),
          selectInput("spsort", "Rank kinases on", choices = list("Median Score", "Max Score", "Statistic")),
          actionButton("switchToRun", "Switch to Run Analysis View"),
          width = 3),
        mainPanel(
          tabsetPanel(
            tabPanel("Upstream Kinase Score",
                     helpText("This plot shows putative upstream kinases ranked by their Final Score (median) or the value of the Kinase Statistic (median)."),
                     helpText("Each point is the result of an individual analysis with a different rank cut-off for adding upstream kinases for peptides."),
                     helpText("The size of the points indicates the size of the peptide set used for a kinase in the corresponding analysis."),
                     helpText("The color of the points indicates the specificity score resulting from the corresponding analysis."),
                     plotOutput("scorePlot", height = "1400px")
            ),
            tabPanel("Kinase Volcano",
                     helpText("This plot shows the median Final Score (y-axis) versus the mean Kinase Statistic (x-axis) of putative upstream kinases."),
                     helpText("The size of the kinase names indicates the size of the peptide set used for a kinase in the corresponding analysis."),
                     helpText("The color of the kinase names indicates the specificity score resulting from the corresponding analysis."),
                     plotOutput("volcanoPlot", height = "800px")
            ),
            tabPanel("Kinase Details",
                     helpText("These are kinase details"),
                     selectInput("showKinase", "Select kinase", choices = ""),
                     tags$hr(),
                     actionLink("uniprot", "..."),
                     tags$hr(),
                     dataTableOutput("kinaseSummary"),
                     tabsetPanel(
                       tabPanel("Details Table",
                                helpText(""),
                                dataTableOutput("kinaseDetails")
                       ),
                       tabPanel("Per peptide plot",
                                helpText(""),
                                plotOutput("perPeptidePlot", height = "800px")
                                
                       )
                     )
            ),
            tabPanel("Report",
                     helpText(""),
                     downloadButton("report", "Generate report"),
                     tags$hr(),
                     helpText("The table below shows the settings that were used for this analysis."),
                     dataTableOutput("InfoSettings"),
                     helpText("The table below shows the summary results of the analysis"),
                     dataTableOutput("SummaryTable")
            ),
            tabPanel("Kinase tree",
                     tags$a(href = "http://kinhub.org/kinmap/","The summary data can be saved as a file that can be used to map the data to a phylogenetic tree using the external websit Kinmap."),
                     helpText(""),
                     downloadButton("saveKinMap", "Save data as KinMap file"),
                     tags$hr(),
                     helpText("The specificty score will be mapped to the symbol size"),
                     textInput("sclow", "Scale score from:", 0),
                     textInput("schigh","Scale score to:", 2),
                     textInput("scmax", "max symbol size", 50),
                     tags$hr(),
                     helpText("The kinase statistic will be mapped to the symbol color"),
                     textInput("stlow", "Scale kinase statistic from:", -1),
                     textInput("stmid", "Scale kinase statistic midpoint:", 0),
                     textInput("sthigh","Scale kinase statistic to:", 1),
                     colourpicker::colourInput("cllow", label = "Low color", value = "green"),
                     colourpicker::colourInput("clmid", label = "Mid color", value = "black"),
                     colourpicker::colourInput("clhigh", label = "High color", value = "red")
            )
          )
        )
      )
    }
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
  
  observeEvent(getView(), {
    if (isRunView(getView())) {
      nid <<- showNotification("Press Start to start the analysis.", duration = NULL, type = "message", closeButton = FALSE)
      updateSliderInput(session, "seqHom", min = min(DB$PepProtein_SeqHomology))
    }
  })
  
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
      
      df <- data_input$data
      
      # checks
      check(NonUniqueDataMapping, df, openUrlOnError = TRUE)
      check(MultipleValuesPerCell, df)
      check(EmptyCells, df)
      
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
      
      showNotification(ui = "Please wait ... running analysis.", id = nid, type = "message", closeButton = FALSE, duration = NULL)
      
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
      
      # save objects in tercen context
      saveData(session, list(df = df, full_result = full_result, settings = settings, DB = DB), "results")
      return("Done")
    })
  })
  
  # RESULT
  grpText = reactive({
    grp = results$df[["grp"]]
    if (!is.null(grp)){
      txt = paste("Grouping factor with levels", levels(grp)[1], "and", levels(grp)[2])
    } else {
      txt = "Grouping factor: none"
    }
  })
  
  
  output$grpText = renderText({
    grpText()
  })
  
  xaxText = function(){
    grp = results$df[["grp"]]
    if (!is.null(grp)){
      txt = paste(">0 indicates higher activity in the",as.character(levels(grp)[2]), "group")
    } else {
      return(NULL)
    }
  }
  
  scorePlot = reactive({
    aSub <- results$full_result %>% filter(nFeatures >= input$minsetsize)
    cs   <- makeScorePlot(aSub, input$spsort)
    xax  <- paste("Normalized kinase statistic (", xaxText(),")", sep = "")
    cs   <- cs + ylab(xax)
  })
  
  perPeptidePlot = reactive({
    makePerPeptidePlot(results$df, results$DB %>% filter(Kinase_Name == input$showKinase))
  })
  
  output$scorePlot = renderPlot({
    print(scorePlot())
  })
  
  aSummary = reactive({
    aSub = results$full_result %>% filter(nFeatures >= input$minsetsize)
    aSum = makeSummary(aSub)
    updateSelectInput(session, "showKinase", choices = aSum$ClassName)
    aSum
  })
  
  volcanoPlot = reactive({
    vp = makeVolcanoPlot(aSummary())
    vp = vp + xlab("Mean Kinase Statistic") + ylab("Median Final Score")
  })
  
  output$volcanoPlot = renderPlot({
    print(volcanoPlot())
  })
  
  output$perPeptidePlot = renderPlot({
    print(perPeptidePlot())
  })
  
  output$kinaseSummary = renderDataTable({
    aSum = aSummary()
    aKin = aSum %>% filter(ClassName == input$showKinase)
    aTable = data.frame(name  = c("Mean Kinase Statistic", "Mean Specificity Score", "Mean Significance Score", "Median Combined Score"),
                        value = c(aKin$meanStat, aKin$meanFeatScore, aKin$meanPhenoScore, aKin$medianScore))
    colnames(aTable) = c(input$showKinase, "value")
    create_datatable(aTable)
  })
  
  observeEvent(input$showKinase, {
    upid = results$kinase2uniprot %>% filter(Kinase_Name == input$showKinase)%>%.$Kinase_UniprotID
    aLabel = paste(input$showKinase," (",upid,") on Uniprot.org", sep = "")
    updateActionButton(session, "uniprot", label = aLabel)
  })
  
  
  detailsTable = reactive({
    makeDetailsTable(results$df, results$DB %>% filter(Kinase_Name == input$showKinase))
  })
  
  output$kinaseDetails = renderDataTable({
    create_datatable(detailsTable())
  })
  
  output$InfoSettings = renderDataTable({
    create_datatable(getSettingsInfo(results$settings))
  })
  
  observeEvent(input$saveDetailsTable, {
    aTable = detailsTable()
    filename = file.path(getFolder(), paste(gsub("/", "-",input$showKinase), format(Sys.time(), "%Y%m%d-%H%M.txt")) )
    write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
    shell.exec(getFolder())
  })
  
  observeEvent(input$uniprot, {
    upid = results$kinase2uniprot %>% filter(Kinase_Name == input$showKinase)%>%.$Kinase_UniprotID
    browseURL(paste("http://www.uniprot.org/uniprot/", upid, sep = ""))
  })
  
  summaryResultTable = reactive({
    df = aSummary()
    df = left_join(results$kinase2uniprot, df, by = c("Kinase_Name" = "ClassName"))
    df = df%>% filter(!is.na(medianScore))
    df = df %>% select(Kinase_UniprotID, Kinase_Name, meanFeatScore, meanPhenoScore, medianScore, maxScore, meanStat, sdStat, meanSetSize)
    #list("Median Score", "Max Score", "Statistic")
    if (input$spsort == "Median Score"){
      df = df %>% arrange(-medianScore)
    } else if(input$spsort == "Max Score"){
      df = df %>% arrange(-maxScore)
    } else if (input$spsort == "Statistic"){
      df = df%>%arrange(meanStat)
    }
    
    df = df %>% dplyr::rename("Kinase Uniprot ID" = Kinase_UniprotID,
                              "Kinase Name" = Kinase_Name,
                              "Mean Specificity Score" = meanFeatScore,
                              "Mean Significance Score" = meanPhenoScore,
                              "Median Final score" = medianScore,
                              "Max Final Score" = maxScore,
                              "Mean Kinase Statistic" = meanStat,
                              "SD Kinase Statitistic" = sdStat,
                              "Mean peptide set size" = meanSetSize)
  })
  
  output$SummaryTable = renderDataTable({
    create_datatable(summaryResultTable())
  })
  
  output$saveKinMap <- downloadHandler(
    filename = function() {
      paste("KinMap file", format(Sys.time(), "%Y%m%d-%H%M.txt"))
    },
    content = function(file) {
      df = summaryResultTable()
      colnames(df) = make.names(colnames(df)) # convert formatted column names to valid variable names
      szFrom = c(as.numeric(input$sclow), as.numeric(input$schigh))
      szTo = c(0, as.numeric(input$scmax))
      szScale = list(from = szFrom, to = szTo)
      clScale = list(low = as.numeric(input$stlow), mid = as.numeric(input$stmid), high = as.numeric(input$sthigh))
      clr = c(input$cllow, input$clmid, input$clhigh)
      treeFile(fname = file, mappings = Kinase.Name ~ Mean.Kinase.Statistic + Mean.Specificity.Score, data = df,szScale = szScale, clScale = clScale, clValues = clr)
    }
  )
  
  output$report <- downloadHandler(
    filename = function() {
      paste('Report', Sys.Date(), '.html', sep='')
    },
    content = function(file) {
      tmp_dir    <- tempdir()
      tempReport <- file.path(tmp_dir, "report.Rmd")
      tempLogo   <- file.path(tmp_dir, "pglogo.png")
      file.copy("rmd/report.Rmd", tempReport, overwrite = TRUE)
      file.copy("rmd/pglogo.png", tempLogo, overwrite = TRUE)
      params <- list(
        scorePlot     = scorePlot(),
        volcanoPlot   = volcanoPlot(),
        summaryTable  = summaryResultTable(),
        settingsTable = getSettingsInfo(results$settings),
        grpText       = grpText()
      )
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
})

getValues <- function(session){
  ctx <- getCtx(session)
  values <- list()
  
  if (!"ID" %in% names(ctx$rnames)) {
    stop("ID field should be available as row name")
  }
  
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

getWorkflowId = function(session){
  workflowId = getOption("tercen.workflowId")
  if (!is.null(workflowId)) return(workflowId)
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["workflowId"]])
}

getStepId = function(session){
  stepId = getOption("tercen.stepId")
  if (!is.null(stepId)) return(stepId)
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["stepId"]])
}

getMode = function(session){
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["mode"]])
}

saveData <- function(session, result_list, name) {
  ctx        <- getCtx(session)
  workflowId <- getWorkflowId(session)
  stepId     <- getStepId(session)
  workflow   <- ctx$client$workflowService$get(workflowId)
  
  fileDoc = FileDocument$new()
  fileDoc$name = name
  fileDoc$projectId = workflow$projectId
  fileDoc$acl$owner = workflow$acl$owner
  fileDoc$metadata$contentType = 'application/octet-stream'
  
  metaWorkflowId = Pair$new()
  metaWorkflowId$key = 'workflow.id'
  metaWorkflowId$value = workflowId
  
  metaStepId = Pair$new()
  metaStepId$key = 'step.id'
  metaStepId$value = stepId
  
  fileDoc$meta <- list(metaWorkflowId, metaStepId)
  con          <- rawConnection(raw(0), "r+")
  # store the query object to be able to check if input data has changed
  result_list  <- append(result_list, list(query = ctx$query))
  saveRDS(result_list, con)
  bytes        <- rawConnectionValue(con)
  fileDoc      <- ctx$client$fileService$upload(fileDoc, bytes)
  fileDoc
}

getResultsFile = function(session, name) {
  result     <- NULL
  ctx        <-  getCtx(session)
  workflowId <- getWorkflowId(session)
  stepId     <- getStepId(session)
  
  files <- ctx$client$fileService$findFileByWorkflowIdAndStepId(
    startKey   = list(workflowId, stepId),
    endKey     = list(workflowId, ''),
    descending = TRUE, limit = 4)
  
  if (length(files) > 0) {
    files <- files[unlist(lapply(files, FUN = function(x) x$name == name))]
    if (!identical(files, list())) {
      result <- files[[1]]
    }
  } 
  result
}

getCtxResults <- function(session) {
  ctx      <-  getCtx(session)
  result   <- NULL
  file     <- getResultsFile(session, "results")
  if (!is.null(file)) {
    bytes    <- ctx$client$fileService$download(file$id)
    raw_con  <- rawConnection(object = bytes, open = "r")
    result   <- readRDS(raw_con)
  }
  result
}

inputDataEqual <- function(ctx, results) {
  result <- FALSE
  if (!is.null(results) && class(results) == "list" && "query" %in% names(results) &&
      ctx$query$qtHash == results$query$qtHash && ctx$query$rowHash == results$query$rowHash && ctx$query$columnHash == results$query$columnHash) {
    result <- TRUE
  }
  result
}

isRunView <- function(view) {
  view == "run"
}
