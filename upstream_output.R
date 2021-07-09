# Upstream output functions
library(ggplot2)
library(scales)

scorePlot = function(aResult, plotorder = "score"){
  aFull = ldply(aResult, .fun = function(.)return(data.frame(.$aResult, mxRank = .$mxRank) ))
  sp = makeScorePlot(aFull, plotorder)
}

makeScorePlot = function (aFrame, plotorder = "score") {
  if (plotorder == "Median Score") {
    aFrame = aFrame %>% mutate(rankedClassName = reorder(ClassName, combinedScore, median))
  } else if(plotorder == "Max Score") {
    aFrame = aFrame %>% mutate(rankedClassName = reorder(ClassName, combinedScore, max))
  } else if (plotorder == "Statistic") {
    aFrame = aFrame %>% mutate(rankedClassName = reorder(ClassName, -NormalizedSetStat, median))
  } else stop("Unknown Ranking Method")
  
  prt2 = ggplot(aFrame, aes(x= rankedClassName, y = NormalizedSetStat)) + geom_boxplot()
  prt2 = prt2 + geom_point(aes(colour = rescale(pFeatureScore, from = c(0,2), clip = TRUE),
                               size   = rescale(nFeatures, from = c(1,30), to = c(1,30), clip = TRUE)))
  prt2 = prt2 + geom_abline(slope = 0)
  prt2 = prt2 + coord_flip() + theme_bw()
  prt2 = prt2 + scale_colour_gradientn(name = "Specificity Score",space = "rgb",
                                       colours = c("black", "white", "red"),
                                       values = c(0,1.3, 2)/2,
                                       breaks = c(0,1.3, 2)/2,
                                       labels = c(0, 1.3, 2),
                                       limits = c(0,1))
  prt2 = prt2 + scale_size(breaks = c(10,20,31), labels = c("10", "20", ">30"))
  prt2 = prt2 + ylab("Normalized kinase statistic") + xlab("Kinase name")
  yRange = layer_scales(prt2)$y$range$range
  prt2 = prt2 + scale_y_continuous(position = "top", limits = c( min(0, yRange[1]), max(0, yRange[2])) )
  prt2 = prt2 + guides( size = guide_legend("Peptide set size"), color = guide_colorbar("Specificity score"))
  return(prt2)
}

makeVolcanoPlot = function(aSummary) {
  vp = ggplot(aSummary, aes(
    x = meanStat,
    y = medianScore,
    colour = rescale(meanFeatScore, from = c(0,2), clip = TRUE),
    label = ClassName,
    size = meanSetSize,
    alpha =  1-rescale(sdStat, from = c(0,1), clip = TRUE )
  )
  )
  
  vp = vp + geom_text() + scale_colour_gradientn(name = "Mean Specificity Score",space = "rgb",
                                                 colours = c("black", "white", "red"),
                                                 values = c(0,1.3, 2)/2,
                                                 breaks = c(0,1.3, 2)/2,
                                                 labels = c(0, 1.3, 2),
                                                 limits = c(0,1))
  
  vp = vp + scale_alpha("Consistency", limits = c(0,1))
  vp = vp + scale_size("Number of peptides")
  vp = vp + theme_bw() + theme(panel.background = element_rect(fill = 'lightyellow', colour = 'black'))
  return(vp)
}

lowestRank = function(dbf) {
  lrf = dbf %>% summarize(lowestRank = min(Kinase_Rank))
  return(lrf$lowestRank)
}

makePerPeptidePlot = function(df, dbFrame, scanRank = NULL, minPScore = NULL) {
  if (!is.null(minPScore)) {
    dbFrame <- subset(dbFrame, Database != "phosphoNET" |  Kinase_PKinase_PredictorVersion2Score > minPScore)
  }
  if (!is.null(scanRank)) {
    dbFrame <- subset(dbFrame, Kinase_Rank <= scanRank)
  }
  ixList  <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df      <- ixList[[2]]
  
  if (!is.null(df[["grp"]])) {
    perPepStats = peptideAnalysis(df)
    perPepStats = perPepStats %>% group_by(ID) %>% mutate(lowestRank = )
    #ppp = ggplot(perPepStats, aes(x = reorder(ID, -lowRank), colour = as.factor(lowRank), y = pes, ymin = cil, ymax = ciu)) + geom_point() + geom_errorbar()
    ppp = ggplot(perPepStats, aes(x = ID, y = pes, ymin = cil, ymax = ciu)) + geom_point() + geom_errorbar()
    ppp = ppp + coord_flip()
    ppp = ppp + xlab("Peptide ID") + ylab("Group difference")
    ppp = ppp + geom_abline(slope = 0)
    ppp = ppp + guides( color = guide_legend("Lowest Kinase Rank"))
    
  } else {
    perPepStats = df %>% group_by(ID) %>% do({
      thisPep = subset(dbFrame, ID == .$ID[1])
      aResult = data.frame(value = .$value,
                           lowRank = thisPep$lowestRank
      )
    })
    ppp = ggplot(perPepStats, aes(x = reorder(ID, -lowRank), colour = as.factor(lowRank), y = value)) + geom_point()
    ppp = ppp + coord_flip()
    ppp = ppp + xlab("Peptide ID") + ylab("Value")
    ppp = ppp + geom_abline(slope = 0)
    ppp = ppp + guides( color = guide_legend("Lowest Kinase Rank"))
  }
  yRange = layer_scales(ppp)$y$range$range
  
  ppp + scale_y_continuous(position = "top", limits = c(-max(abs(yRange)), max(abs(yRange))))
}

makeDetailsTable = function(df, dbFrame, scanRank = NULL, minPScore = NULL) {
  if (!is.null(minPScore)) {
    dbFrame= subset(dbFrame, Database != "phosphoNET" |  Kinase_PKinase_PredictorVersion2Score > minPScore)
  }
  if (!is.null(scanRank)) {
    dbFrame = subset(dbFrame, Kinase_Rank <= scanRank)
  }
  ixList = intersectById(dbFrame, df)
  dbFrame = ixList[[1]]
  outFrame = dbFrame %>% group_by(ID, PepProtein_PhosLink, PepProtein_UniprotName) %>%
    dplyr::summarise(Database = paste(Database, collapse = " / "), Kinase_Rank = min(Kinase_Rank) )
  
  outFrame %>% arrange(Kinase_Rank, -as.integer(ID))
}

makeSummary = function(df) {
  aSum = df %>% group_by(ClassName) %>% dplyr::summarise(meanFeatScore = mean(pFeatureScore, na.rm = TRUE),
                                                         maxFeatScore = max(pFeatureScore, na.rm = TRUE),
                                                         meanPhenoScore = mean(pPhenoScore, na.rm = TRUE),
                                                         maxPhenoScore = max(pPhenoScore, na.rm = TRUE),
                                                         meanScore = mean(combinedScore, na.rm = TRUE),
                                                         medianScore = median(combinedScore, na.rm = TRUE),
                                                         maxScore = max(combinedScore, na.rm = TRUE),
                                                         meanStat = mean(NormalizedSetStat, na.rm = TRUE),
                                                         medianStat = median(NormalizedSetStat, na.rm = TRUE),
                                                         sdStat   = sd(NormalizedSetStat, na.rm = TRUE),
                                                         meanSetSize =round(mean(nFeatures, na.rm = TRUE)))
  aSum %>% arrange(-medianScore)
}

getSettingsInfo = function(settings) {
  return(settings)
}

create_datatable <- function(data) {
  datatable(data, extensions = 'Buttons', options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
}
