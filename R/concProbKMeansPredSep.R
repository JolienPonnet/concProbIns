concProbKMeansPredSep <- function(inputDT, lowCat, highCat, nClusMaxLow, nClusMaxHigh, gamma, fine = FALSE, letsTime = TRUE){
  checkDT(inputDT, c('exposure', 'observed', 'predicted'))
  checkLogicVec(list(letsTime))
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
  
  expZero = round(inputZero$exposure,4)
  expOne = round(inputOne$exposure,4)
  expZeroUniq = unique(expZero)
  expOneUniq = unique(expOne)
  
  if((length(expZeroUniq)>=length(expOneUniq) & fine == TRUE) | (length(expZeroUniq)<length(expOneUniq) & fine == FALSE)){
    compPairs = rep(0, length(expZeroUniq))
    concPairs = rep(0, length(expZeroUniq))
    concProbs = rep(0, length(expZeroUniq))
    for(expInd in 1:length(expZeroUniq)){
      selExp0 = which(expZero==expZeroUniq[expInd])
      selPred0 = inputZero$predicted[selExp0]
      
      selExp1 = which((expOne>= round(expZeroUniq[expInd] - gamma,4)) & (expOne<= round(expZeroUniq[expInd] + gamma,4)))
      if(length(selExp1)==0) next()
      
      selPred1 <- inputOne$predicted[selExp1]
      
      nClusZero <- max(min(nClusMaxLow, (length(selPred0)-1)),1)
      nClusOne <- max(min(nClusMaxHigh, (length(selPred1)-1)),1)
      zeroClus <- suppressWarnings(kmeans(selPred0, nClusZero))
      oneClus <- suppressWarnings(kmeans(selPred1, nClusOne))
      
      selPredZero <- zeroClus$centers
      selPredOne <- oneClus$centers
      
      predClusZero <- zeroClus$cluster # The cluster to which each prediction belongs
      predClusOne <- oneClus$cluster # The cluster to which each prediction belongs
      selNPredZero <- table(predClusZero) # Number of predictions that belong to each cluster
      selNPredOne <- table(predClusOne) # Number of predictions that belong to each cluster
      
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
    uniqExp = expZeroUniq
  }else{
    compPairs = rep(0, length(expOneUniq))
    concPairs = rep(0, length(expOneUniq))
    concProbs = rep(0, length(expOneUniq))
    for(expInd in 1:length(expOneUniq)){
      selExp1 = which(expOne==expOneUniq[expInd])
      selPred1 = inputOne$predicted[selExp1]
      
      selExp0 = which((expZero>= round(expOneUniq[expInd] - gamma,4)) & (expZero<= round(expOneUniq[expInd] + gamma,4)))
      if(length(selExp0)==0) next()
      
      selPred0 <- inputZero$predicted[selExp0]
      
      nClusZero <- max(min(nClusMaxLow, (length(selPred0)-1)),1)
      nClusOne <- max(min(nClusMaxHigh, (length(selPred1)-1)),1)
      zeroClus <- suppressWarnings(kmeans(selPred0, nClusZero))
      oneClus <- suppressWarnings(kmeans(selPred1, nClusOne))
      
      selPredZero <- zeroClus$centers
      selPredOne <- oneClus$centers
      
      predClusZero <- zeroClus$cluster # The cluster to which each prediction belongs
      predClusOne <- oneClus$cluster # The cluster to which each prediction belongs
      selNPredZero <- table(predClusZero) # Number of predictions that belong to each cluster
      selNPredOne <- table(predClusOne) # Number of predictions that belong to each cluster
      
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
    uniqExp = expOneUniq = unique(expOne)
  }
  globalConcProb = sum(concPairs)/sum(compPairs) 
  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = globalConcProb, concProbs = concProbs, compPairs = compPairs, uniqExp = uniqExp, time = time))
  } else {
    return(list(concProbGlobal = globalConcProb, concProbs = concProbs, compPairs = compPairs, uniqExp = uniqExp))
  }
}