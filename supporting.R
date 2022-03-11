# Supporting functions

create_UI_title <- function() {
  tags$div(HTML("<strong><font color = #6895d1>Upstream Kinase Analysis - January 2022</font></strong>"))
}

# initial view for show and run mode
createInitialView <- function(disableRun = FALSE) {
  button <- actionButton("start", "Start")
  if (disableRun) {
    button <- disabled(button)
  }
  sidebarLayout(
    sidebarPanel(
      create_UI_title(),
      tags$hr(),
      button,
      tags$hr(),
      textOutput("status")
    ),
    mainPanel(tabsetPanel(
      tabPanel("Basic Settings",
               tags$hr(),
               selectInput("kinasefamily", label = "Kinase family", choices = c("PTK", "STK")),
               sliderInput("scan", "Scan Rank From-To", min = 1, max = 20, value = c(4,12),round = TRUE),
               sliderInput("nperms", "Number of permutations", min = 100, max = 1000, value = 500, step = 100, round = TRUE)
      ),
      
      tabPanel("Advanced Settings",
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
}

xaxText <- function(results) {
  grp = results$df[["grp"]]
  if (!is.null(grp)){
    txt = paste(">0 indicates higher activity in the", as.character(levels(grp)[2]), "group")
  } else {
    return(NULL)
  }
}

getGroupingLabel <- function(values) {
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