# Supporting functions

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