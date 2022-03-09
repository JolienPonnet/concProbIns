concProbGrid <- function(inputDT, lowCat, highCat, quantSplits, gamma, fine = FALSE, letsTime = TRUE){
  checkDT(inputDT, c('exposure', 'observed', 'predicted'))
  checkLogicVec(list(letsTime))
  checkNumOrIntVec(list(quantSplits))
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

      splits <- sort(unique(quantile(c(selPred0, selPred1), quantSplits))) # Define the split points
      groups <- groupConstruct(selPred0, selPred1, splits) # Count the number of observations in each interval
      counts <- countConc(groups$zeroGroup, groups$oneGroup) # Determine the elements of the sum in (3.3), without division by the total nb of elements in each group
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

      splits <- sort(unique(quantile(c(selPred0, selPred1), quantSplits))) # Define the split points
      groups <- groupConstruct(selPred0, selPred1, splits) # Count the number of observations in each interval
      counts <- countConc(groups$zeroGroup, groups$oneGroup) # Determine the elements of the sum in (3.3), without division by the total nb of elements in each group
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

groupConstruct <- function(selPred0, selPred1, splits){
  vecX0 <- sort(selPred0)
  vecX1 <- sort(selPred1)

  zeroGroup = rep(NA,(length(splits) + 1))
  oneGroup = rep(NA,(length(splits) + 1))

  stepSize <- 100

  for(iSplit in 1:(length(splits))){
    seed <- max(round(iSplit/length(splits)*length(vecX0)),1)
    ind0 <- findFirstPos(vecX0, splits[iSplit], seed, stepSize)
    seed <- max(round(iSplit/length(splits)*length(vecX1)),1)
    ind1 <- findFirstPos(vecX1, splits[iSplit], seed, stepSize)

    zeroGroup[iSplit] <- ind0 - 1
    oneGroup[iSplit] <- ind1 - 1
  }
  zeroGroup <- c(zeroGroup[1], diff(zeroGroup[-length(zeroGroup)]), length(vecX0) - zeroGroup[length(splits)])
  oneGroup <- c(oneGroup[1], diff(oneGroup[-length(oneGroup)]), length(vecX1) - oneGroup[length(splits)])
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

countConc <- function(zeroGroup, oneGroup){
  nSplits <- length(zeroGroup)
  nLower <- zeroGroup * c(0,cumsum(oneGroup)[-nSplits])
  nHigher <- zeroGroup * c(rev(cumsum(rev(oneGroup)))[-1],0)
  list(nHigher = nHigher, nLower = nLower)
}
