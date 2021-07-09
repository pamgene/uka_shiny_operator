# Upstream Analysis functionality
library(combinat)
library(data.table)

setStats = function(M, classMat, grp, statfun){
  # statistic per set
  setStat = t(classMat) %*% statfun(M, grp)
  return(setStat)
}

stat.snr = function(M, grp, pair = NULL){
  # snr statistic per peptide
  bGrp1 = grp == levels(grp)[1]
  bGrp2 = grp == levels(grp)[2]
  m1 = apply(M[bGrp1,], 2, mean)
  s1 = apply(M[bGrp1,], 2, var)
  m2 = apply(M[bGrp2,], 2, mean)
  s2 = apply(M[bGrp2,], 2, var)
  aStat = (m2 - m1) / sqrt(s1+s2)
  # mask any constant columns by 0
  aStat[is.nan(aStat)] = 0
  return(aStat)
}

pScore = function(dataMatrix, classMatrix, phenoGrp, upStat, permFun, statFun = stat.snr, nPerms = DFT.PERMS){
  pSetStat = permFun(dataMatrix, classMatrix, grp = phenoGrp, statFun = statFun, nPerms = nPerms)
  nPerms = dim(pSetStat)[2]
  myScore = vector(length = length(upStat))
  for (i in 1:length(upStat)){
    myScore[i] = max(sum(abs(pSetStat[i,])>=abs(upStat[i]))/nPerms, 1/nPerms)
  }
  return(myScore)
}

mkPerms = function(grp, maxPerms = DFT.PERMS){
  if ( length(grp) < 7){
    pIdx = permn(length(grp))
    p = matrix(nrow = length(grp), ncol = length(pIdx), data = unlist(pIdx))
    if (length(pIdx) > maxPerms){
      p = p[, sample(1:maxPerms)]
    }
  } else {
    mIdx = matrix(nrow = length(grp), ncol = maxPerms, data = 1:length(grp))
    p = apply(mIdx, 2, sample)
  }
  return(p)
}

fcsPhenoPerms = function(dataMatrix, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  grpPerms = mkPerms(grp, nPerms)
  nPerms = dim(grpPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pGrp = grp[grpPerms[,i]]
    pStat[,i] = setStats(dataMatrix, classMatrix, pGrp, statFun)
  }
  return(pStat)
}

fcsSpotPerms = function(X, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  if(is.matrix(X)){
    dataMatrix = X
  } else {
    dataMatrix = matrix(nrow = 1, ncol = length(X), data = X)
  }
  
  spotPerms = mkPerms(1:dim(dataMatrix)[2], nPerms)
  nPerms = dim(spotPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pMatrix  = dataMatrix[, spotPerms[,i]]
    pStat[,i] = setStats(pMatrix, classMatrix, grp, statFun)
  }
  return(pStat)
}

fcs <- function(dataMatrix, classMatrix, phenoGrp, statFun = stat.snr, phenoPerms = TRUE, featurePerms = TRUE, nPerms = DFT.PERMS){
  upStats = setStats(dataMatrix, classMatrix, phenoGrp, statFun)
  nClass = apply(classMatrix > 0 ,2, sum)
  aFcsResult = data.frame(ClassName = colnames(classMatrix), nFeatures = nClass, SetStat = upStats, NormalizedSetStat = upStats/nClass)
  if(phenoPerms){
    aPhenoScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = fcsPhenoPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aPhenoScore = NaN
  }
  if(featurePerms){
    aFeatureScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = fcsSpotPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aFeatureScore = NaN
  }
  aFcsResult = data.frame(aFcsResult, phenoScore = aPhenoScore,
                          featureScore = aFeatureScore,
                          pPhenoScore = -log10(aPhenoScore),
                          pFeatureScore = -log10(aFeatureScore) )
  combinedScore = apply(as.matrix(aFcsResult[, c("pPhenoScore", "pFeatureScore")]),1,sum, na.rm = TRUE)
  aFcsResult = data.frame(aFcsResult, combinedScore = combinedScore)
  return(aFcsResult)
}

pgScanAnalysis2g = function(df, dbFrame,
                            dbWeights = c(iviv = 1,PhosphoNET = 1),
                            scanRank,
                            nPermutations = 500) {
  #run two group.
  #add dbWeight
  dbFrame <- dbFrame %>% 
    group_by(Database) %>% 
    do({data.frame(., dbWeight = dbWeights[[.$Database[1]]])})
  
  ixList  <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df      <- ixList[[2]]
  X       <- acast(df, .ci~ID, fun.aggregate = mean, value.var = ".y")
  grp     <- acast(df,  .ci~ID, value.var = "grp")[,1]
  grp     <- as.factor(grp)
  
  result <- lapply(scanRank, FUN = function(i) {
    db_top <- subset(dbFrame, Kinase_Rank <= i)
    M      <- acast(db_top, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] = 0
    inx2   <- intersect(colnames(X), rownames(M))
    Xi     <- X[, colnames(X) %in% inx2]
    M      <- M[rownames(M) %in% inx2,]
    M      <- M[order(rownames(M)),]
    Xi     <- Xi[, order(colnames(Xi))]
    if (!all(rownames(M)== colnames(Xi))) stop("Mismatch between data matrix and upstream kinase matrix")
    
    result <- fcs(Xi, M, grp, nPerms = nPermutations)
    result <- result[order(result$combinedScore, decreasing = TRUE),]
    data.table(mxRank = i, result = list(result), X = list(Xi), M = list(M))
  })
  result
}

pgScanAnalysis0 = function(df, dbFrame,
                           dbWeights = c(HPRD = 1,PhosphoNET = 1,Phosphosite = 1,Reactome = 1),
                           scanRank,
                           nPermutations = 500){
  
  # run a sinlge column, without grouping
  #add dbWeight
  dbFrame <- dbFrame %>% 
    group_by(Database) %>% do({ data.frame(., dbWeight = dbWeights[[.$Database[1]]]) })
  
  ixList  <- intersectById(dbFrame, df)
  dbFrame <- ixList[[1]]
  df      <- ixList[[2]]
  
  result <- lapply(scanRank, FUN = function(i) {
    aTop = subset(dbFrame, Kinase_Rank <= i)
    M = acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
    M[!is.finite(M)] = 0
    inx2 = intersect(df$ID, rownames(M))
    dfx =  df%>%filter(ID %in% inx2)
    X = matrix(ncol = dim(dfx)[1], nrow = 1, data = dfx$value)
    colnames(X) = dfx$ID
    M = M[rownames(M) %in% inx2,]
    M = M[order(rownames(M)),]
    X = X[, order(colnames(X))]
    ##
    if (!all(rownames(M)== colnames(X))) stop("Mismatch between data matrix and upstream kinase matrix")
    result <- fcs(X, M, statFun = stat.identity, phenoGrp = NULL, phenoPerms = FALSE, nPerms = nPermutations)
    result <- result[order(result$combinedScore, decreasing = TRUE),]
    data.table(mxRank = i, result = list(result), X = list(X), M = list(M))
  })
  result
}


intersectById = function(df1, df2){
  df1$ID = droplevels(as.factor(df1$ID))
  df2$ID = droplevels(as.factor(df2$ID))
  isct = intersect(df1$ID, df2$ID)
  if (length(isct) == 0){
    stop("No matching peptide IDs!")
  }
  df1 = subset(df1, ID %in% isct)
  df2 = subset(df2, ID %in% isct)
  df1$ID = droplevels(df1$ID)
  df2$ID = droplevels(df2$ID)
  return(list(df1, df2))
}

peptideAnalysis = function(df, type = "unpaired"){
  if( type == "unpaired"){
    result = unpaired(df)
  } else {
    stop(paste("unknown type for peptide analysis:", type))
  }
}

unpaired = function(df) {
  df = df %>% group_by(ID) %>% do({
    aTest = t.test(.y ~ grp, data = ., var.equal = TRUE)
    result = data.frame(pes = diff(aTest$estimate),
                        ciu = -aTest$conf.int[1],
                        cil = -aTest$conf.int[2],
                        pval = aTest$p.value)
  })
}
