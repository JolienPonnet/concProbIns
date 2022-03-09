concProbKMeansExpPredSep <- function(inputDT, lowCat, highCat, nClusMaxLow, nClusMaxHigh, gamma, fine = FALSE, letsTime = TRUE){
  checkDT(inputDT, c('exposure', 'observed', 'predicted'))
  checkLogicVec(list(letsTime, fine))
  checkNumOrIntVec(list(lowCat, highCat, nClusMaxLow, nClusMaxHigh, gamma))
  lowCat <- round(lowCat, digits = 0)
  highCat <- round(highCat, digits = 0)
  nClusMaxLow <- round(nClusMaxLow, digits = 0)
  nClusMaxHigh <- round(nClusMaxHigh, digits = 0)
  checkRanges(list(lowCat[1], highCat[1], nClusMaxLow, nClusMaxHigh), list(c('>=', 0), c('>=', 0), c('>=', 1), c('>=', 1)))
  checkLength(list(nClusMaxLow), c(1))
  checkLength(list(nClusMaxHigh), c(1))
  checkLength(list(gamma), c(1))
  if(!length(lowCat) %in% c(1, 2)) stop('The lowCat argument can only be of length 1 or 2.')
  if(!length(highCat) %in% c(1, 2)) stop('The highCat argument can only be of length 1 or 2.')
  if(length(lowCat) == 2) checkRanges(list(lowCat[2]), list(c('>=', 0)))
  if(length(highCat) == 2) checkRanges(list(highCat[2]), list(c('>=', 0)))

  observed <- predicted <- exposure <- NULL

  lowCat <- sort(unique(lowCat))
  highCat <- sort(unique(highCat))

  if(!(max(lowCat) < min(highCat))) stop('Please make sure that the intervals of the lowCat and the highCat intervals do not intersect.')

  if (letsTime) ptm <- proc.time()

  dataTemp <- inputDT

  if(length(lowCat)==1){
    if(length(highCat)==1){
      inputZero <- dataTemp[observed == lowCat, ]
      inputOne <- dataTemp[observed == highCat, ]
    }else{
      inputZero <- dataTemp[observed == lowCat, ]
      inputOne <- dataTemp[observed >= highCat[1] & observed <= highCat[2], ]
    }
  }else{
    if(length(highCat)==1){
      inputZero <- dataTemp[observed >= lowCat[1] & observed <= lowCat[2], ]
      inputOne <- dataTemp[observed == highCat, ]
    }else{
      inputZero <- dataTemp[observed >= lowCat[1] & observed <= lowCat[2], ]
      inputOne <- dataTemp[observed >= highCat[1] & observed <= highCat[2], ]
    }
  }
  inputZero[,observed:=NULL]
  inputOne[,observed:=NULL]

  nClusZeroPred <- max(min(nClusMaxLow, (dim(inputZero)[1]-1)),1)
  nClusOnePred <- max(min(nClusMaxHigh, (dim(inputOne)[1]-1)),1)
  nClusZeroExp <- max(min(nClusMaxLow, (dim(unique(inputZero[,'exposure']))[[1]]-1)),1)
  nClusOneExp <- max(min(nClusMaxHigh, (dim(unique(inputOne[,'exposure']))[[1]]-1)),1)

  zeroClusPred <- suppressWarnings(kmeans(inputZero[,'predicted'], nClusZeroPred))
  oneClusPred <- suppressWarnings(kmeans(inputOne[,'predicted'], nClusOnePred))
  zeroClusExp <- suppressWarnings(kmeans(inputZero[,'exposure'], nClusZeroExp))
  oneClusExp <- suppressWarnings(kmeans(inputOne[,'exposure'], nClusOneExp))

  clusZeroPred <- zeroClusPred$cluster # The cluster to which each prediction belongs
  clusOnePred <- oneClusPred$cluster # The cluster to which each prediction belongs
  clusZeroExp <- zeroClusExp$cluster # The cluster to which each exposure belongs
  clusOneExp <- oneClusExp$cluster # The cluster to which each exposure belongs

  predZero <- zeroClusPred$centers # The center of each cluster
  predOne <- oneClusPred$centers # The center of each cluster
  expZero <- zeroClusExp$centers # The center of each cluster
  expOne <- oneClusExp$centers # The center of each cluster

  if((length(expZero)>=length(expOne) & fine == TRUE) | (length(expZero)<length(expOne) & fine == FALSE)){
    compPairs = rep(0, length(expZero))
    concPairs = rep(0, length(expZero))
    concProbs = rep(0, length(expZero))
    for(expInd in 1:length(expZero)){
      selExp0 = which(clusZeroExp == expInd)
      selPredClusZero = clusZeroPred[selExp0]
      selNPredZero <- table(selPredClusZero) # Number of predictions that belong to each cluster
      selPredZero <- predZero[as.numeric(rownames(selNPredZero))] # The prediction center of each selected cluster

      selExp1 = which((expOne>= round(expZero[expInd] - gamma,4)) & (expOne<= round(expZero[expInd] + gamma,4)))
      if(length(selExp1)==0) next()

      selExp1 = which(clusOneExp %in% selExp1)
      selPredClusOne = clusOnePred[selExp1]
      selNPredOne <- table(selPredClusOne) # Number of predictions that belong to each cluster
      selPredOne <- predOne[as.numeric(rownames(selNPredOne))] # The prediction center of each selected cluster

      for (iZero in 1:length(selPredZero)) {
        for (iOne in 1:length(selPredOne)) {
          if (selPredZero[iZero] < selPredOne[iOne]) {
            concPairs[expInd]  <- concPairs[expInd] + as.numeric(selNPredZero[iZero]) * as.numeric(selNPredOne[iOne])
          } else if (selPredZero[iZero] > selPredOne[iOne]) {
            compPairs[expInd]  <- compPairs[expInd] + as.numeric(selNPredZero[iZero]) * as.numeric(selNPredOne[iOne])
          }
        }
      }
      compPairs[expInd] = concPairs[expInd] + compPairs[expInd]
      concProbs[expInd] = concPairs[expInd] / compPairs[expInd]
    }
    uniqExp = expZero
  }else{
    compPairs = rep(0, length(expOne))
    concPairs = rep(0, length(expOne))
    concProbs = rep(0, length(expOne))
    for(expInd in 1:length(expOne)){
      selExp1 = which(clusOneExp == expInd)
      selPredClusOne = clusOnePred[selExp1]
      selNPredOne <- table(selPredClusOne) # Number of predictions that belong to each cluster
      selPredOne <- predOne[as.numeric(rownames(selNPredOne))] # The prediction center of each selected cluster

      selExp0 = which((expZero>= round(expOne[expInd] - gamma,4)) & (expZero<= round(expOne[expInd] + gamma,4)))
      if(length(selExp0)==0) next()

      selExp0 = which(clusZeroExp %in% selExp0)
      selPredClusZero = clusZeroPred[selExp0]
      selNPredZero <- table(selPredClusZero) # Number of predictions that belong to each cluster
      selPredZero <- predZero[as.numeric(rownames(selNPredZero))] # The prediction center of each selected cluster

      for (iZero in 1:length(selPredZero)) {
        for (iOne in 1:length(selPredOne)) {
          if (selPredZero[iZero] < selPredOne[iOne]) {
            concPairs[expInd]  <- concPairs[expInd] + as.numeric(selNPredZero[iZero]) * as.numeric(selNPredOne[iOne])
          } else if (selPredZero[iZero] > selPredOne[iOne]) {
            compPairs[expInd]  <- compPairs[expInd] + as.numeric(selNPredZero[iZero]) * as.numeric(selNPredOne[iOne])
          }
        }
      }
      compPairs[expInd] = concPairs[expInd] + compPairs[expInd]
      concProbs[expInd] = concPairs[expInd] / compPairs[expInd]
    }
    uniqExp = expOne
  }
  globalConcProb = sum(concPairs)/sum(compPairs)
  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = globalConcProb, concProbs = concProbs, compPairs = compPairs, uniqExp = uniqExp, time = time))
  } else {
    return(list(concProbGlobal = globalConcProb, concProbs = concProbs, compPairs = compPairs, uniqExp = uniqExp))
  }
}
