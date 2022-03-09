concProbGridLowHigh <- function(inputDT, lowCat, highCat, quantSplitsLow, quantSplitsHigh, gamma, fine = FALSE, letsTime = TRUE){
  checkDT(inputDT, c('exposure', 'observed', 'predicted'))
  checkLogicVec(list(letsTime))
  checkNumOrIntVec(list(quantSplitsLow, quantSplitsHigh))
  checkNumOrIntVec(list(lowCat, highCat, gamma))
  lowCat <- round(lowCat, digits = 0)
  highCat <- round(highCat, digits = 0)
  checkRanges(list(lowCat[1], highCat[1]), list(c('>=', 0), c('>=', 0)))
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
    for(iExp in 1:length(expZeroUniq)){
      selExp0 = which(expZero==expZeroUniq[iExp])
      selPred0 = inputZero$predicted[selExp0]

      selExp1 = which((expOne>= round(expZeroUniq[iExp] - gamma,4)) & (expOne<= round(expZeroUniq[iExp] + gamma,4)))
      if(length(selExp1)==0) next()

      selPred1 <- inputOne$predicted[selExp1]

      splitsLow <- sort(unique(quantile(selPred0, quantSplitsLow))) # Define the split points
      splitsHigh <- sort(unique(quantile(selPred1, quantSplitsHigh))) # Define the split points
      groups <- groupConstructLowHigh(selPred0, selPred1, splitsLow, splitsHigh) # Count the number of observations in each interval
      counts <- countConcLowHigh(groups$zeroGroup, groups$oneGroup, splitsLow, splitsHigh) # Determine the elements of the sum in (3.3), without division by the total nb of elements in each group
      concPairs[iExp] <- sum(counts$nHigher)
      compPairs[iExp] <- sum(counts$nLower) + concPairs[iExp]
      concProbs[iExp] <- concPairs[iExp]/compPairs[iExp]
    }
    uniqExp = expZeroUniq
  }else{
    compPairs = rep(0, length(expOneUniq))
    concPairs = rep(0, length(expOneUniq))
    concProbs = rep(0, length(expOneUniq))
    for(iExp in 1:length(expOneUniq)){
      selExp1 = which(expOne==expOneUniq[iExp])
      selPred1 = inputOne$predicted[selExp1]

      selExp0 = which((expZero>= round(expOneUniq[iExp] - gamma,4)) & (expZero<= round(expOneUniq[iExp] + gamma,4)))
      if(length(selExp0)==0) next()

      selPred0 <- inputZero$predicted[selExp0]

      splitsLow <- sort(unique(quantile(selPred0, quantSplitsLow))) # Define the split points
      splitsHigh <- sort(unique(quantile(selPred1, quantSplitsHigh))) # Define the split points
      groups <- groupConstructLowHigh(selPred0, selPred1, splitsLow, splitsHigh) # Count the number of observations in each interval
      counts <- countConcLowHigh(groups$zeroGroup, groups$oneGroup, splitsLow, splitsHigh) # Determine the elements of the sum in (3.3), without division by the total nb of elements in each group
      concPairs[iExp] <- sum(counts$nHigher)
      compPairs[iExp] <- sum(counts$nLower) + concPairs[iExp]
      concProbs[iExp] <- concPairs[iExp]/compPairs[iExp]
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

groupConstructLowHigh <- function(selPred0, selPred1, splitsLow, splitsHigh){
  vecX0 <- sort(selPred0)
  vecX1 <- sort(selPred1)

  zeroGroup = rep(NA,(length(splitsLow) + 1))
  oneGroup = rep(NA,(length(splitsHigh) + 1))

  stepSize <- 100

  for(iSplit in 1:min(length(splitsLow),length(splitsHigh))){
    seed <- max(round(iSplit/length(splitsLow)*length(vecX0)),1)
    ind0 <- findFirstPos(vecX0, splitsLow[iSplit], seed, stepSize)
    seed <- max(round(iSplit/length(splitsHigh)*length(vecX1)),1)
    ind1 <- findFirstPos(vecX1, splitsHigh[iSplit], seed, stepSize)

    zeroGroup[iSplit] <- ind0 - 1
    oneGroup[iSplit] <- ind1 - 1
  }

  if(length(splitsLow)<length(splitsHigh)){
    for(iSplit in (length(splitsLow)+1):length(splitsHigh)){
      seed <- max(round(iSplit/length(splitsHigh)*length(vecX1)),1)
      ind1 <- findFirstPos(vecX1, splitsHigh[iSplit], seed, stepSize)

      oneGroup[iSplit] <- ind1 - 1
    }
  }else if(length(splitsLow)>length(splitsHigh)){
    for(iSplit in (length(splitsHigh)+1):length(splitsLow)){
      seed <- max(round(iSplit/length(splitsLow)*length(vecX0)),1)
      ind0 <- findFirstPos(vecX0, splitsLow[iSplit], seed, stepSize)

      zeroGroup[iSplit] <- ind0 - 1
    }
  }
  zeroGroup <- c(zeroGroup[1], diff(zeroGroup[-length(zeroGroup)]), length(vecX0) - zeroGroup[length(splitsLow)])
  oneGroup <- c(oneGroup[1], diff(oneGroup[-length(oneGroup)]), length(vecX1) - oneGroup[length(splitsHigh)])
  list(zeroGroup = zeroGroup, oneGroup = oneGroup)
}

findFirstPos <- function(vec, targetVal, seed, stepSize){
  if(stepSize == 1) stepSize = 2
  vec = as.vector(vec)
  indTemp <- seed
  if(vec[indTemp] == targetVal) return(indTemp)
  if(vec[indTemp] < targetVal){
    while((vec[indTemp] <= targetVal) & (indTemp < length(vec)) ){
      indTemp <- min(indTemp + stepSize, length(vec))
    }

    indTemp <- max(indTemp - stepSize + 1, 1)
    while((vec[indTemp] < targetVal) & (indTemp < length(vec)) ){
      indTemp <- indTemp + 1
    }
  } else {
    while((vec[indTemp] > targetVal) & (indTemp > 1)){
      indTemp <- max(indTemp - stepSize, 1)
    }

    indTemp <- min(indTemp + stepSize - 1, length(vec))
    while((vec[indTemp] >= targetVal) & indTemp >= 1){
      indTemp <- indTemp - 1
      if(indTemp == 0) break
    }
    indTemp = indTemp + 1
  }
  return(indTemp)
}

countConcLowHigh <- function(zeroGroup, oneGroup, splitsLow, splitsHigh){
  if(length(splitsLow)<=length(splitsHigh)){
    nHigher <- rep(0,length(splitsLow)+1)
    nLower <- rep(0,length(splitsLow)+1)
    cumsumOne <- cumsum(oneGroup)
    cumsumRevOne <- rev(cumsum(rev(oneGroup)))
    stepSize = max(round(length(splitsHigh)/20),2)
    for (k in 1:length(splitsLow)){
      ind <- findFirstPos(splitsHigh, splitsLow[k], min(k, length(splitsHigh)), stepSize)#Geeft index van eerste getal >= targetVal
      if(k != length(zeroGroup)){
        if(splitsHigh[ind]>= splitsLow[k]){
          nHigher[k] = zeroGroup[k] * cumsumRevOne[ind+1]
        }
      }

      if((ind==1) & (splitsLow[k] < splitsHigh[ind])) next()
      if(splitsLow[k] >= splitsHigh[ind]){ #findfirstpos faalde, kon geen hoger of gelijk aan vinden / gelijk aan
        nLower[k+1] = zeroGroup[k+1] * cumsumOne[ind]
      }else{
        if(ind > 1) nLower[k+1] = zeroGroup[k+1] * cumsumOne[ind-1]
      }
    }
  }else{
    nHigher <- rep(0,length(splitsHigh)+1)
    nLower <- rep(0,length(splitsHigh)+1)
    cumsumZero <- cumsum(zeroGroup)
    cumsumRevZero <- rev(cumsum(rev(zeroGroup)))
    stepSize = max(round(length(splitsHigh)/20),2)
    for (k in 1:length(splitsHigh)){
      ind <- findFirstPos(splitsLow, splitsHigh[k], min(k, length(splitsLow)), stepSize)#Geeft index van eerste getal >= targetVal
      if((ind==1) & (splitsLow[ind] > splitsHigh[k])) next()
      if(splitsLow[ind] <= splitsHigh[k]){
        nHigher[k+1] = oneGroup[k+1] * cumsumZero[ind]
      }else{
        if(ind > 1) nHigher[k+1] = oneGroup[k+1] * cumsumZero[ind-1]
      }

      if(k != length(oneGroup)){
        if(splitsLow[ind]>= splitsHigh[k]){
          nLower[k] = oneGroup[k] * cumsumRevZero[ind+1]
        }
      }
    }
  }
  return(list(nHigher = nHigher, nLower = nLower))
}
