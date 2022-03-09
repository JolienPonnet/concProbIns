concProbContKMeans <- function(inputDT, nu = 0, nClus, letsTime = TRUE, log = FALSE, samplePart = FALSE, sampleSize = 10000, allData = TRUE){

  checkDT(inputDT, c('observed', 'predicted'))
  checkLogicVec(list(letsTime, log))
  checkNumOrIntVec(list(nu, nClus))
  nClus <- round(nClus, digits = 0)
  checkRanges(list(nu, nClus), list(c('>=', 0), c('>', 0)))

  observed <- predicted <- "." <- NULL

  if(log == TRUE){
    inputDT$observed <- log(inputDT$observed)
    inputDT$predicted <- log(inputDT$predicted)
  }

  if (letsTime) ptm <- proc.time()

  dataTemp <- inputDT
  clusFit <- suppressWarnings(kmeans(dataTemp[, .(predicted, observed)], nClus, iter.max = 20))

  clusCents <- clusFit$centers
  clusSize <- clusFit$size

  obsCol <- which(colnames(clusCents) == 'observed')
  predCol <- which(colnames(clusCents) == 'predicted')

  ch <- rep(0, nClus)
  dh <- rep(0, nClus)

  for(iClus in 1:(nClus-1)){
    for(jClus in (iClus+1):nClus){
      if(abs(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) > nu){
        if( (sign(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) == sign(clusCents[iClus, predCol] - clusCents[jClus, predCol])) & ((sign(clusCents[iClus, obsCol] - clusCents[jClus, obsCol]) != 0)) ){
          ch[iClus] <- ch[iClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
          ch[jClus] <- ch[jClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
        } else if ((clusCents[iClus, obsCol] != clusCents[jClus, obsCol]) &( clusCents[iClus, predCol] != clusCents[jClus, predCol])){ #remove ties
          dh[iClus] <- dh[iClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
          dh[jClus] <- dh[jClus] + as.numeric(clusSize[iClus])*as.numeric(clusSize[jClus])
        }
      }
    }
  }

  concProb <- sum(ch)/(sum(ch) + sum(dh))

  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = concProb, time = time))
  } else {
    return(list(concProbGlobal = concProb))
  }
}
