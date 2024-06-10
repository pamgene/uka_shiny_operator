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
source("kinase_tree_parser.R")
source("supporting.R")

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
  message <- reactiveValues(text = "")
  DB <- readRDS("data/240516_86312_86402_86412_87102_87202_UpstreamDb.rds") %>% select(-PepProtein_SeqIdentity)

  dataInput <- reactive({
    getValues(session)
  })

  propertiesInput <- reactive({
    getProperties(session)
  })

  mode <- reactive({
    getMode(session)
  })

  summaryResultTable = reactive({
    df = aSummary()
    df = left_join(results$kinase2uniprot, df, by = c("Kinase_Name" = "ClassName"))
    group_names <- DB %>% ungroup() %>% distinct(Kinase_Name, Kinase_group, Kinase_family)
    df <- left_join(df, group_names, by = c("Kinase_Name" = "Kinase_Name"))
    df = df %>% filter(!is.na(medianScore))
    df = df %>% select(Kinase_UniprotID, Kinase_Name, Kinase_group, Kinase_family, meanFeatScore, meanPhenoScore, medianScore, maxScore, meanStat, medianStat, sdStat, meanSetSize)
    #list("Median Score", "Max Score", "Statistic")
    if (input$spsort == "Median Score") {
      df = df %>% arrange(-medianScore)
    } else if (input$spsort == "Max Score") {
      df = df %>% arrange(-maxScore)
    } else if (input$spsort == "Statistic") {
      df = df %>% arrange(meanStat)
    }

    df = df %>% dplyr::rename("Kinase Uniprot ID" = Kinase_UniprotID,
                              "Kinase Name" = Kinase_Name,
                              "Kinase Group" = Kinase_group,
                              "Kinase Family" = Kinase_family,
                              "Mean Specificity Score" = meanFeatScore,
                              "Mean Significance Score" = meanPhenoScore,
                              "Median Final score" = medianScore,
                              "Max Final Score" = maxScore,
                              "Mean Kinase Statistic" = meanStat,
                              "Median Kinase Statistic" = medianStat,
                              "SD Kinase Statitistic" = sdStat,
                              "Mean peptide set size" = meanSetSize)
  })

  scorePlot = reactive({
    aSub <- results$full_result %>% filter(nFeatures >= input$minsetsize)
    group_names <- DB %>% distinct(Kinase_Name, Kinase_group)
    aSub <- left_join(aSub, group_names, by = c("ClassName" = "Kinase_Name"))

    cs <- makeScorePlot(aSub, input$spsort)
    xax <- paste("Normalized kinase statistic (", xaxText(results), ")", sep = "")

    axis_order <- ggplot_build(cs)$layout$panel_params[[1]]$y$get_labels()
    kinLabelGroup <- aSub %>%
          distinct(ClassName, Kinase_group) %>%
          arrange(match(ClassName, axis_order)) %>%
          pull(Kinase_group)
    
    properties <- propertiesInput()

    if (properties$Kinase_family == "PTK") {
      kinLabelColors <- ifelse((kinLabelGroup == "TYR") | (kinLabelGroup == "TKL"), "black", "red")
    } else if (properties$Kinase_family == "STK") {
      kinLabelColors <- ifelse((kinLabelGroup == "TYR") | (kinLabelGroup == "TKL"), "red", "black")
    }
    
    cs + ylab(xax) + theme(axis.text.y = element_text(colour = kinLabelColors))
  })

  perPeptidePlot = reactive({
    makePerPeptidePlot(results$df, results$DB %>% filter(Kinase_Name == input$showKinase))
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

  detailsTable = reactive({
    makeDetailsTable(results$df, results$DB %>% filter(Kinase_Name == input$showKinase))
  })

  grpText = reactive({
    grp = results$df[["grp"]]
    if (!is.null(grp)) {
      txt = paste("Grouping factor with levels", levels(grp)[1], "and", levels(grp)[2])
    } else {
      txt = "Grouping factor: none"
    }
  })

  ## Outputs

  output$body <- renderUI({
    mode <- mode()
    if (isShowView(mode)) {
      createInitialView(disableRun = TRUE)
    } else if (isRunView(mode)) {
      createInitialView()
    } else if (isResultView(mode)) {
      computedResults <- getCtxResults(session)
      if (!is.null(computedResults)) {
        results$DB <- computedResults$DB
        results$full_result <- computedResults$full_result
        results$df <- computedResults$df
        results$settings <- computedResults$settings

        results$kinase2uniprot <- results$DB %>%
          group_by(Kinase_UniprotID) %>%
          dplyr::summarise(Kinase_Name = Kinase_Name[1])

        sidebarLayout(
          sidebarPanel(
            create_UI_title(),
            tags$hr(),
            textOutput("grpText"),
            tags$hr(),
            sliderInput("minsetsize", "Include results based on peptide sets with minimal size of:",
                        min = 1,
                        max = 5,
                        step = 1,
                        value = 3),
            width = 3),
          mainPanel(
            tabsetPanel(
              tabPanel("Upstream Kinase Score",
                       helpText("This plot shows putative upstream kinases ranked by their Final Score (median) or the value of the Kinase Statistic (median)."),
                       helpText("Each point is the result of an individual analysis with a different rank cut-off for adding upstream kinases for peptides."),
                       helpText("The size of the points indicates the size of the peptide set used for a kinase in the corresponding analysis."),
                       helpText("The color of the points indicates the specificity score resulting from the corresponding analysis."),
                       helpText("Red colored kinase labels on the y-axis mean the kinases could be dual-specificity (both PTK and STK). Caution is advised, check literature to see whether the result is valid."),
                       selectInput("spsort", "Sort score plot on", choices = list("Median Score", "Max Score", "Statistic")),
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
                       textOutput("kinGroup"),
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
                       helpText("The table below shows the settings that were used for this analysis."),
                       dataTableOutput("InfoSettings"),
                       helpText("The table below shows the summary results of the analysis"),
                       dataTableOutput("SummaryTable")
              )
        # ,
        # tabPanel("Kinase tree",
        #          tags$a(href = "http://kinhub.org/kinmap/","The summary data can be saved as a file that can be used to map the data to a phylogenetic tree using the external websit Kinmap."),
        #          helpText(""),
        #          downloadButton("saveKinMap", "Save data as KinMap file"),
        #          tags$hr(),
        #          helpText("The specificty score will be mapped to the symbol size"),
        #          textInput("sclow", "Scale score from:", 0),
        #          textInput("schigh","Scale score to:", 2),
        #          textInput("scmax", "max symbol size", 50),
        #          tags$hr(),
        #          helpText("The kinase statistic will be mapped to the symbol color"),
        #          textInput("stlow", "Scale kinase statistic from:", -1),
        #          textInput("stmid", "Scale kinase statistic midpoint:", 0),
        #          textInput("sthigh","Scale kinase statistic to:", 1),
        #          colourpicker::colourInput("cllow", label = "Low color", value = "green"),
        #          colourpicker::colourInput("clmid", label = "Mid color", value = "black"),
        #          colourpicker::colourInput("clhigh", label = "High color", value = "red")
        # )
            )
          )
        )
      }
    }
  })

  output$grpText = renderText({
    grpText()
  })

  output$scorePlot = renderPlot({
    print(scorePlot())
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
    aTable = data.frame(name = c("Mean Kinase Statistic", "Mean Specificity Score", "Mean Significance Score", "Median Combined Score"),
                        value = c(aKin$meanStat, aKin$meanFeatScore, aKin$meanPhenoScore, aKin$medianScore))
    colnames(aTable) = c(input$showKinase, "value")
    create_datatable(aTable)
  })

  output$kinaseDetails = renderDataTable({
    create_datatable(detailsTable())
  })

  output$InfoSettings = renderDataTable({
    create_datatable(getSettingsInfo(results$settings))
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
      treeFile(fname = file, mappings = Kinase.Name ~ Mean.Kinase.Statistic + Mean.Specificity.Score, data = df, szScale = szScale, clScale = clScale, clValues = clr)
    }
  )

  ## Observe (event)
  observe({
    ctx <- getCtx(session)
    data_input <- dataInput()
    properties <- propertiesInput()

    if (is.null(properties)) return()
    updateSelectInput(session, "kinasefamily", selected = properties$Kinase_family)

    if (is.null(data_input)) return()

    if (length(ctx$cnames) > 0) {
      poptions <- unlist(ctx$cnames)
      updateSelectInput(session, "pairingFactor", choices = poptions)
    }

    if (!is.null(input$start) && input$start == 0) {
      grp = getGroupingLabel(data_input)
      if (properties$Lock_kinase_family == "Yes") {
        shinyjs::disable("kinasefamily")
      }
      if (!is.null(grp)) {
        message$text <- paste("Grouping factor:", grp)
      } else {
        message$text <- "Grouping factor: none"
      }
    }
    output$status = renderText({ message$text })
  })

  observeEvent(input$showKinase, {
    upid = results$kinase2uniprot %>% filter(Kinase_Name == input$showKinase) %>% .$Kinase_UniprotID
    aLabel = paste(input$showKinase, " (", upid, ") on Uniprot.org", sep = "")
    updateActionButton(session, "uniprot", label = aLabel)
  })

  output$kinGroup <- renderText({
    upid <- results$kinase2uniprot %>%
        filter(Kinase_Name == input$showKinase) %>%
        .$Kinase_UniprotID

    kinGroup <- DB %>%
        filter(Kinase_UniprotID == upid) %>%
        ungroup() %>%
        distinct(Kinase_group) %>%
        pull()
    kinGroupStr <- paste("Kinase", input$showKinase, "belongs to the", kinGroup, "kinase group.")
  })

  observeEvent(mode(), {
    mode <- mode()
    if (isRunView(mode)) {
      nid <<- showNotification("Press Start to start the analysis.", duration = NULL, type = "message", closeButton = FALSE)
      updateSliderInput(session, "seqHom", min = min(DB$PepProtein_SeqSimilarity))
    } else if (isResultView(mode)) {
      if (exists("nid")) {
        removeNotification(nid)
      }
    }
  })

  observeEvent(input$saveDetailsTable, {
    aTable = detailsTable()
    filename = file.path(getFolder(), paste(gsub("/", "-", input$showKinase), format(Sys.time(), "%Y%m%d-%H%M.txt")))
    write.table(aTable, filename, sep = "\t", quote = FALSE, row.names = FALSE)
    shell.exec(getFolder())
  })

  observeEvent(input$uniprot, {
    upid = results$kinase2uniprot %>% filter(Kinase_Name == input$showKinase) %>% .$Kinase_UniprotID
    browseURL(paste("http://www.uniprot.org/uniprot/", upid, sep = ""))
  })

  # calculation
  observeEvent(input$start, {
    lapply(c("start", "kinasefamily", "scan", "nperms", "wIviv", "wPhosphoNET", "seqHom", "minPScore", "usePairing", "pairingFactor"), shinyjs::disable)

    ctx <- getCtx(session)
    properties <- propertiesInput()
    data_input <- dataInput()
    grp <- getGroupingLabel(data_input)
    df <- data_input$data

    # remove peptides with missing data
    zIdx <- df %>% group_by(.ri) %>% tally() %>% filter(n != max(n)) %>% pull(.ri)
    if (!is.null(zIdx) && length(zIdx) > 0) {
      df <- df %>% filter(!(.ri %in% zIdx))
      msg <- paste("Warning:", length(zIdx), "peptides with missing data have been removed from the data")
      showNotification(ui = msg, duration = NULL, type = "warning")
    }

    # checks
    check(NonUniqueDataMapping, df, openUrlOnError = TRUE)
    check(MultipleValuesPerCell, df)
    check(EmptyCells, df)

    if (input$kinasefamily == "PTK") {
      DB <- DB %>% filter(PepProtein_Residue == "Y")
    } else if (input$kinasefamily == "STK") {
      DB <- DB %>% filter(PepProtein_Residue == "S" | PepProtein_Residue == "T")
    } else {
      stop("Unknown value for kinase family")
    }

    DB <- DB %>%
      filter(PepProtein_SeqSimilarity >= input$seqHom) %>%
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
                          value = c(input$kinasefamily, input$scan[1], input$scan[2], input$nperms, input$wIviv, input$wPhosphoNET, input$minPScore, input$seqHom))

    # save objects in tercen context
    saveData(session, list(df = df, full_result = full_result, settings = settings, DB = DB), "results")

    lapply(c("start", "scan", "nperms", "wIviv", "wPhosphoNET", "seqHom", "minPScore", "usePairing", "pairingFactor"), shinyjs::enable)
    if (properties$Lock_kinase_family == "No") { shinyjs::enable("kinasefamily") }
    message$text = "Done"
  })

})

getValues <- function(session) {
  ctx <- getCtx(session)
  values <- list()


  if (length(names(ctx$rnames)) != 1) {
    stop("A single row name should be defined.")
  }

  row_name <- ctx %>%
    rselect() %>%
    colnames()
  data <- ctx %>%
    select(.y, .ri, .ci, .sids, .tlbIdx) %>%
    mutate(color = ctx$select(ctx$colors) %>% pull)
  row_data <- ctx %>%
    rselect() %>%
    rename(ID = !!row_name) %>%
    mutate(.ri = seq(0, length(unique(data$.ri)) - 1))

  values$data <- data %>%
    left_join(., row_data)

  colors <- ctx$colors
  if (identical(colors, list())) {
    values$colorLabels <- c()
  } else {
    values$colorLabels <- colnames(ctx$select(colors))
  }

  return(values)
}

getProperties <- function(session) {
  ctx <- getCtx(session)

  list(Kinase_family = ifelse(is.null(ctx$op.value("Kinase_family")), "PTK", ctx$op.value("Kinase_family")),
       Lock_kinase_family = ifelse(is.null(ctx$op.value("Lock_kinase_family")), "Yes", ctx$op.value("Lock_kinase_family")))
}

getWorkflowId = function(session) {
  workflowId = getOption("tercen.workflowId")
  if (!is.null(workflowId)) return(workflowId)
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["workflowId"]])
}

getStepId = function(session) {
  stepId = getOption("tercen.stepId")
  if (!is.null(stepId)) return(stepId)
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["stepId"]])
}

getMode = function(session) {
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["mode"]])
}

saveData <- function(session, result_list, name) {
  ctx <- getCtx(session)
  workflowId <- getWorkflowId(session)
  stepId <- getStepId(session)
  workflow <- ctx$client$workflowService$get(workflowId)

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
  con <- rawConnection(raw(0), "r+")
  # store the query object to be able to check if input data has changed
  result_list <- append(result_list, list(query = ctx$query))
  saveRDS(result_list, con)
  bytes <- rawConnectionValue(con)
  removeCurrentFiles(session, workflowId, stepId)
  fileDoc <- ctx$client$fileService$upload(fileDoc, bytes)

  # call save to indicate step has a result
  ctx %>%
    select(.y, .ci, .ri) %>%
    group_by(.ci, .ri) %>%
    summarise(mean = mean(.y)) %>%
    ctx$addNamespace() %>%
    ctx$save()
}

getFilesByWorkflowAndStep <- function(ctx, workflowId, stepId) {
  ctx$client$fileService$findFileByWorkflowIdAndStepId(
    startKey = list(workflowId, stepId),
    endKey = list(workflowId, ''),
    descending = TRUE, limit = 10)
}

getResultsFile = function(session, name) {
  result <- NULL
  ctx <- getCtx(session)
  workflowId <- getWorkflowId(session)
  stepId <- getStepId(session)

  files <- getFilesByWorkflowAndStep(ctx, workflowId, stepId)

  if (length(files) > 0) {
    files <- files[unlist(lapply(files, FUN = function(x) x$name == name))]
    if (!identical(files, list())) {
      result <- files[[1]]
    }
  }
  result
}

removeCurrentFiles <- function(session, workflowId, stepId) {
  ctx <- getCtx(session)
  files <- getFilesByWorkflowAndStep(ctx, workflowId, stepId)
  for (file in files) {
    ctx$client$fileService$delete(file$id, file$rev)
  }
}

getCtxResults <- function(session) {
  ctx <- getCtx(session)
  result <- NULL
  file <- getResultsFile(session, "results")
  if (!is.null(file)) {
    bytes <- ctx$client$fileService$download(file$id)
    raw_con <- rawConnection(object = bytes, open = "r")
    result <- readRDS(raw_con)
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

isShowView <- function(mode) {
  !is.null(mode) && mode == "show"
}
isRunView <- function(mode) {
  !is.null(mode) && mode == "run"
}
isResultView <- function(mode) {
  !is.null(mode) && mode == "showResult"
}
