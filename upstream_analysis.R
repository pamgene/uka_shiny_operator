# Upstream Analysis functionality

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
  nCores  <- detectCores()
  cl      <- makeCluster(nCores)
  registerDoParallel(cl)
  #print(paste("Please wait, processing on ", nCores, " CPU cores ..."))
  aPackagesList = c("pgFCS", "reshape2" )
  
  aScanResult =  foreach(i = scanRank, .packages = aPackagesList) %dopar%
    {
      aTop = subset(dbFrame, Kinase_Rank <= i)
      
      M = acast(aTop, ID ~ Kinase_Name, value.var = "dbWeight", fun.aggregate = max)
      M[!is.finite(M)] = 0
      inx2 = intersect(colnames(X), rownames(M))
      Xi = X[, colnames(X) %in% inx2]
      M = M[rownames(M) %in% inx2,]
      # cannot depend on dcast for correct order
      M = M[order(rownames(M)),]
      Xi = Xi[, order(colnames(Xi))]
      if(!all(rownames(M)== colnames(Xi))) stop("Mismatch between data matrix and upstream kinase matrix")
      aResult = fcs(Xi, M, grp, nPerms = nPermutations)
      aResult = aResult[order(aResult$combinedScore, decreasing = TRUE),]
      return(data.table(mxRank = i, aResult = list(aResult), X = list(Xi), M = list(M)))
    }
  stopCluster(cl)
  return(aScanResult)
}

#'@export
pgScanAnalysis0 = function(df     ,dbFrame,
                           dbWeights = c(HPRD = 1,PhosphoNET = 1,Phosphosite = 1,Reactome = 1),
                           scanRank,
                           nPermutations = 500){
  
  # run a sinlge column, without grouping
  #add dbWeight
  dbFrame <- dbFrame %>% 
    group_by(Database) %>% do({ data.frame(., dbWeight = dbWeights[[.$Database[1]]]) })
  
  browser()
  ixList = intersectById(dbFrame, df)
  dbFrame = ixList[[1]]
  df = ixList[[2]]
  
  nCores = detectCores()
  cl = makeCluster(nCores)
  registerDoParallel(cl)
  #print(paste("Please wait, processing on ", nCores, " CPU cores ..."))
  aPackagesList = c("pgFCS", "reshape2" )
  
  aScanResult =  foreach(i = scanRank, .packages = aPackagesList) %dopar%
    #aScanResult = for(i in scanRank)#debug
    {
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
      if(!all(rownames(M)== colnames(X))) stop("Mismatch between data matrix and upstream kinase matrix")
      aResult = fcs(X, M, statFun = stat.identity, phenoGrp = NULL, phenoPerms = FALSE, nPerms = nPermutations)
      aResult = aResult[order(aResult$combinedScore, decreasing = TRUE),]
      return(data.table(mxRank = i, aResult = list(aResult), X = list(X), M = list(M)))
    }
  stopCluster(cl)
  return(aScanResult)
}


intersectById = function(df1, df2){
  df1$ID = droplevels(df1$ID)
  df2$ID = droplevels(df2$ID)
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

#'@export
peptideAnalysis = function(df, type = "unpaired"){
  if( type == "unpaired"){
    result = unpaired(df)
  } else {
    stop(paste("unknown type for peptide analysis:", type))
  }
}

unpaired = function(df){
  df = df %>% group_by(ID) %>% do({
    aTest = t.test(value ~ grp, data = ., var.equal = TRUE)
    result = data.frame(pes = diff(aTest$estimate),
                        ciu = -aTest$conf.int[1],
                        cil = -aTest$conf.int[2],
                        pval = aTest$p.value)
  })
}
