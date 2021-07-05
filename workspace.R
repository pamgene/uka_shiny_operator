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
  
  modeInput <- reactive({ 
    getMode(session) 
  })
  
  output$body <- renderUI({
    
    mode <- modeInput()
    if (mode == "run") {
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
    } else if (mode == "showResult") {
      DB <- getResults(session, "DB")
      
      kinase2uniprot <- DB %>%
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
          width = 3),
        mainPanel(
          tabsetPanel(
            tabPanel("Upstream Kinase Score",
                     helpText("This plot shows putative upstream kinases ranked by their Final Score (median) or the value of the Kinase Statistic (median)."),
                     helpText("Each point is the result of an individual analysis with a different rank cut-off for adding upstream kinases for peptides."),
                     helpText("The size of the points indicates the size of the peptide set used for a kinase in the corresponding analysis."),
                     helpText("The color of the points indicates the specificity score resulting from the corresponding analysis."),
                     actionLink("saveScorePlot", "Save score plot"),
                     plotOutput("scorePlot", height = "1400px")
            ),
            tabPanel("Kinase Volcano",
                     helpText("This plot shows the median Final Score (y-axis) versus the mean Kinase Statistic (x-axis) of putative upstream kinases."),
                     helpText("The size of the kinase names indicates the size of the peptide set used for a kinase in the corresponding analysis."),
                     helpText("The color of the kinase names indicates the specificity score resulting from the corresponding analysis."),
                     actionLink("saveVolcanoPlot", "Save volcano plot"),
                     plotOutput("volcanoPlot", height = "800px")
            ),
            tabPanel("Kinase Details",
                     helpText("These are kinase details"),
                     selectInput("showKinase", "Select kinase", choices = ""),
                     tags$hr(),
                     actionLink("uniprot", "..."),
                     tags$hr(),
                     tableOutput("kinaseSummary"),
                     tabsetPanel(
                       tabPanel("Details Table",
                                helpText(""),
                                actionLink("saveDetailsTable", "Save details table"),
                                tableOutput("kinaseDetails")
                       ),
                       tabPanel("Per peptide plot",
                                helpText(""),
                                actionLink("savePerpeptidePlot", "Save per peptide plot"),
                                plotOutput("perPeptidePlot", height = "800px")
                                
                       )
                     )
            ),
            tabPanel("Report",
                     helpText(""),
                     downloadButton("report", "Generate report"),
                     tags$hr(),
                     helpText("The table below shows the settings that were used for this analysis."),
                     tableOutput("InfoSettings"),
                     helpText("The table below shows the summary results of the analysis"),
                     actionLink("saveSummaryResults", "Save summary results"),
                     tableOutput("SummaryTable")
            ),
            tabPanel("Kinase tree",
                     tags$a(href = "http://kinhub.org/kinmap/","The summary data can be saved as a file that can be used to map the data to a phylogenetic tree using the external websit Kinmap."),
                     helpText(""),
                     actionButton("saveKinMap", "Save data as KinMap file"),
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
  
  # read database file
  DB  <- readRDS("db.rds")
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
      
      nCores <- detectCores()
      msg    <- paste("Please wait ... running analysis. Using", nCores, "cores.")
      showNotification(ui = msg, id = nid, type = "message", closeButton = FALSE, duration = NULL)
      
      if (input$seed) {
        set.seed(42)
      }
      
      # if (!is.null(grp)) {
      #   df$grp <- as.factor(df$color)
      #   result <- pgScanAnalysis2g(df, dbFrame = DB,
      #                              scanRank = input$scan[1]:input$scan[2],
      #                              nPermutations = input$nperms,
      #                              dbWeights = c(iviv = input$wIviv,
      #                                            PhosphoNET = input$wPhosphoNET
      #                              ))
      # } else {
      #   result <- pgScanAnalysis0(df, dbFrame = DB,
      #                             scanRank = input$scan[1]:input$scan[2],
      #                             nPermutations = input$nperms,
      #                             dbWeights = c(iviv = input$wIviv,
      #                                           PhosphoNET = input$wPhosphoNET
      #                             ))
      # }
      showNotification(ui = "Done", id = nid, type = "message", closeButton = FALSE)
      #full_result <- ldply(result, .fun = function(.) return(data.frame(.$result, mxRank = .$mxRank)))
      
      settings = data.frame(setting = c("Kinase family", "ScanRank Min", "ScanRank Max", "Number of Permutations", "In Vitro In Vitro weight", "PhosphoNET weight", "Min PhosphoNet score", "Min Sequence Homology"),
                            value   = c(input$kinasefamily, input$scan[1] , input$scan[2], input$nperms, input$wIviv, input$wPhosphoNET, input$minPScore, input$seqHom) )
      
      # TODO save objects in tercen context!
      # spath = file.path(getFolder(), "runData.RData")
      saveData(session, df, "data")
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
  
  # RESULT
  grpText = reactive({
    grp = df[["grp"]]
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
    grp = df[["grp"]]
    if (!is.null(grp)){
      txt = paste(">0 indicates higher activity in the",as.character(levels(grp)[2]), "group")
    } else {
      return(NULL)
    }
  }
  
  scorePlot = reactive({
    aSub = aFull %>% filter(nFeatures >= input$minsetsize)
    cs = makeScorePlot(aSub, input$spsort)
    xax = paste("Normalized kinase statistic (", xaxText(),")", sep = "")
    cs = cs + ylab(xax)
  })
  
  perPeptidePlot = reactive({
    
    aPlot = makePerPeptidePlot(df, DB %>% filter(Kinase_Name == input$showKinase))
    return(aPlot)
  })
  
  output$scorePlot = renderPlot({
    print(scorePlot())
  })
  
  aSummary = reactive({
    aSub = aFull %>% filter(nFeatures >= input$minsetsize)
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
  
  output$kinaseSummary = renderTable({
    aSum = aSummary()
    aKin = aSum %>% filter(ClassName == input$showKinase)
    aTable = data.frame(name  = c("Mean Kinase Statistic", "Mean Specificity Score", "Mean Significance Score", "Median Combined Score"),
                        value = c(aKin$meanStat, aKin$meanFeatScore, aKin$meanPhenoScore, aKin$medianScore))
    colnames(aTable) = c(input$showKinase, "value")
    return(aTable)
  })
  
  observeEvent(input$showKinase, {
    upid = kinase2uniprot %>% filter(Kinase_Name == input$showKinase)%>%.$Kinase_UniprotID
    aLabel = paste(input$showKinase," (",upid,") on Uniprot.org", sep = "")
    updateActionButton(session, "uniprot", label = aLabel)
  })
  
  
  detailsTable = reactive({
    aTable = makeDetailsTable(df, DB %>% filter(Kinase_Name == input$showKinase))
  })
  
  output$kinaseDetails = renderTable({
    aTable = detailsTable()
  })
  
  output$InfoSettings = renderTable({
    getSettingsInfo(settings)
  })
  
  observeEvent(input$saveDetailsTable, {
    aTable = detailsTable()
    filename = file.path(getFolder(), paste(gsub("/", "-",input$showKinase), format(Sys.time(), "%Y%m%d-%H%M.txt")) )
    write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
    shell.exec(getFolder())
  })
  
  observeEvent(input$saveScorePlot, {
    filename = file.path(getFolder(), paste("UpstreamScorePlot", format(Sys.time(), "%Y%m%d-%H%M.png")) )
    ggsave(filename, scorePlot(), device = "png" ,units = "cm", height = 40, width = 28)
    shell.exec(getFolder())
  })
  
  observeEvent(input$saveVolcanoPlot, {
    filename = file.path(getFolder(), paste("UpstreamVolcanoPlot", format(Sys.time(), "%Y%m%d-%H%M.png")) )
    ggsave(filename, volcanoPlot(), device = "png" ,units = "cm", height = 30, width = 30)
    shell.exec(getFolder())
  })
  
  observeEvent(input$savePerpeptidePlot, {
    filename = file.path(getFolder(), paste(gsub("/", "-",input$showKinase),"_peptides", format(Sys.time(), "%Y%m%d-%H%M.png")) )
    ggsave(filename, perPeptidePlot(), device = "png" ,units = "cm", height = 40, width = 28)
    shell.exec(getFolder())
  })
  
  observeEvent(input$uniprot, {
    upid = kinase2uniprot %>% filter(Kinase_Name == input$showKinase)%>%.$Kinase_UniprotID
    browseURL(paste("http://www.uniprot.org/uniprot/", upid, sep = ""))
  })
  
  summaryResultTable = reactive({
    df = aSummary()
    df = left_join(kinase2uniprot, df, by = c("Kinase_Name" = "ClassName"))
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
  
  observeEvent(input$saveSummaryResults, {
    aTable = summaryResultTable()
    filename = file.path(getFolder(), paste("Summaryresults", format(Sys.time(), "%Y%m%d-%H%M.txt")) )
    write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
    shell.exec(getFolder())
  })
  
  output$SummaryTable = renderTable({
    summaryResultTable()
  })
  
  observeEvent(input$saveKinMap, {
    df = summaryResultTable()
    colnames(df) = make.names(colnames(df)) # convert formatted column names to valid variable names
    filename = file.path(getFolder(), paste("KinMap file", format(Sys.time(), "%Y%m%d-%H%M.txt")) )
    szFrom = c(as.numeric(input$sclow), as.numeric(input$schigh))
    szTo = c(0, as.numeric(input$scmax))
    szScale = list(from = szFrom, to = szTo)
    clScale = list(low = as.numeric(input$stlow), mid = as.numeric(input$stmid), high = as.numeric(input$sthigh))
    clr = c(input$cllow, input$clmid, input$clhigh)
    print(colnames(df))
    treeFile(fname = filename, mappings = Kinase.Name ~ Mean.Kinase.Statistic + Mean.Specificity.Score, data = df,szScale = szScale, clScale = clScale, clValues = clr)
    shell.exec(getFolder())
  })
  
  output$report <- downloadHandler(
    filename = function() {
      paste('Report', Sys.Date(), '.html', sep='')
    },
    content = function(file) {
      tempReport <- file.path(getFolder(), "report.Rmd")
      templatePath = system.file("rmd","Report.Rmd",package = "UpstreamApp")
      print(tempReport)
      file.copy(templatePath, tempReport, overwrite = TRUE)
      params = list(
        scorePlot = scorePlot(),
        volcanoPlot = volcanoPlot(),
        summaryTable = summaryResultTable(),
        settingsTable = getSettingsInfo(settings),
        grpText = grpText()
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

saveData <- function(session, df, name) {
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
  write.csv(df, file=con, row.names = F)
  bytes        <- rawConnectionValue(con)
  fileDoc      <- ctx$client$fileService$upload(fileDoc, bytes)
  fileDoc
}

getResultsFile = function(session) {
  ctx        <-  getCtx(session)
  workflowId <- getWorkflowId(session)
  stepId     <- getStepId(session)
  
  files <- ctx$client$fileService$findFileByWorkflowIdAndStepId(
    startKey   = list(workflowId, stepId),
    endKey     = list(workflowId, ''),
    descending = TRUE, limit = 1)
  
  if (length(files) > 0) {
    return (files[[1]])
  } 
  
  return (NULL)
}

getResults <- function(session, name) {
  ctx      <-  getCtx(session)
  result   <- NULL
  file     <- getResultsFile()
  if (!is.null(file)) {
    bytes    <- ctx$client$fileService$download(file$id)
    raw_con  <- rawConnection(object = bytes, open = "r")
    bin_data <- readBin(raw_con, what = "raw", n = length(bytes))
    result   <- read.csv(file = textConnection(object = rawToChar(bin_data)))
  }
  result
}

runApp(shinyApp(ui, server))  