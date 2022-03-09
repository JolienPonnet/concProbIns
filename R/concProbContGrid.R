concProbContGrid <- function (inputDT, nu = 0, quantSplits, direction = 'left', letsTime = TRUE, log = FALSE){
  checkDT(inputDT, c('observed', 'predicted'))
  checkLogicVec(list(letsTime, log))
  checkNumOrIntVec(list(nu, quantSplits))
  checkRanges(list(nu), list(c('>=', 0)))
  checkCharVec(list(direction))
  checkValues(list(direction), list(c('up', 'down', 'left', 'right', 'all')))
  
  observed <- predicted <- NULL
  
  if(log == TRUE){
    inputDT$observed <- log(inputDT$observed)
    inputDT$predicted <- log(inputDT$predicted)
  } 
  
  if (letsTime)
    ptm <- proc.time()
  
  resMat = list()
  splits = list()
  
  dataTemp <- inputDT
  
  splits <- sort(unique(quantile(dataTemp[, observed], quantSplits)))
  
  numbMat <- groupConstructCont(dataTemp, splits)
  resMat <- countConcCont(numbMat, direction, nu, splits)
  
  nConc <- as.numeric(resMat$concMat)
  nDisc <- as.numeric(resMat$discMat)
  
  concProb <- sum(nConc)/(sum(nDisc) + sum(nConc))
  
  if (letsTime) {
    time <- proc.time() - ptm
    return(list(concProbGlobal = concProb, time = time, concProbMat = resMat, splits = splits)) 
  } else {
    return(list(concProbGlobal = concProb, concProbMat = resMat, splits = splits))
  }
}

countConcCont <- function(numbMat, direction = 'left', nu, splits){
  checkCharVec(list(direction))
  checkValues(list(direction), list(c('up', 'down', 'left', 'right', 'all')))
  
  dimMat <- nrow(numbMat)
  concMat <- matrix(0, nrow = dimMat, ncol = dimMat)
  discMat <- matrix(0, nrow = dimMat, ncol = dimMat)
  
  for(iRow in 1:(dimMat)){
    for(iCol in 1:dimMat){
      if(numbMat[iRow, iCol] != 0){
        if(nu>0){
          if(iRow == 1){
            ind <- which(splits >= (splits[1] + nu))[1]
            if(is.na(ind)) ind <- dimMat
            rows2BeRemoved <- 1:pmin(ind, dimMat)
          }else if(iRow == dimMat){
            ind <- which(splits >= (splits[dimMat-1] - nu) )[1]
            rows2BeRemoved <- pmax(1,ind):dimMat
          }else{
            indMin <- which((splits) >= (splits[iRow] + nu))[1]
            if(is.na(indMin)) indMin <- dimMat
            indMax <- which((splits) >= (splits[iRow-1] - nu))[1]
            rows2BeRemoved <- pmax(1,indMax):pmin(indMin, dimMat)
          }
        }else{
          rows2BeRemoved <- iRow
        }
        
        if(direction == 'left'){
          nConc <- sum(numbMat[-(min(rows2BeRemoved):dimMat),-(iCol:dimMat)])*as.numeric(numbMat[iRow, iCol])
          nDisc <- sum(numbMat[-(min(rows2BeRemoved):dimMat),-(1:iCol)])*as.numeric(numbMat[iRow, iCol])
        } else if (direction == 'right'){
          nConc <- sum(numbMat[-(1:max(rows2BeRemoved)),-(1:iCol)])*as.numeric(numbMat[iRow, iCol])
          nDisc <- sum(numbMat[-(1:max(rows2BeRemoved)),-(iCol:dimMat)])*as.numeric(numbMat[iRow, iCol])
        } else if (direction == 'down'){
          nConc <- sum(numbMat[-(min(rows2BeRemoved):dimMat),-(iCol:dimMat)])*as.numeric(numbMat[iRow, iCol])
          nDisc <- sum(numbMat[-(1:max(rows2BeRemoved)),-(iCol:dimMat)])*as.numeric(numbMat[iRow, iCol])
        } else if (direction == 'up'){
          nConc <- sum(numbMat[-(1:max(rows2BeRemoved)),-(1:iCol)])*as.numeric(numbMat[iRow, iCol])
          nDisc <- sum(numbMat[-(min(rows2BeRemoved):dimMat),-(1:iCol)])*as.numeric(numbMat[iRow, iCol])
        } else if(direction == 'all'){
          nConc <-  (sum(numbMat[-(min(rows2BeRemoved):dimMat),-(iCol:dimMat)]) + sum(numbMat[-(1:max(rows2BeRemoved)),-(1:iCol)]))*as.numeric(numbMat[iRow, iCol])
          nDisc <- (sum(numbMat[-(1:max(rows2BeRemoved)),-(iCol:dimMat)]) + sum(numbMat[-(min(rows2BeRemoved):dimMat),-(1:iCol)]))*as.numeric(numbMat[iRow, iCol])
        }
        concMat[iRow, iCol] <- nConc
        discMat[iRow, iCol] <- nDisc
      }
    }
  }
  return(list(concMat = concMat, discMat = discMat))
}

groupConstructCont <- function(inputDT, splits){
  nSplits <- length(splits)
  groupMat <- matrix(0, nrow = nSplits + 1, ncol = nSplits + 1)
  orderedPred <- setorder(inputDT, "predicted")
  ind0 <- 0
  stepSize <- 100
  stopBool <- FALSE
  for(iSplit in 1:(nSplits+1)){
    seed <- max(round(iSplit/length(splits)*nrow(orderedPred)),1)
    if(iSplit == (nSplits + 1)){
      ind1 <- nrow(orderedPred)  
    }else{
      testIndPred <- findFirstPos(unlist(orderedPred[,"predicted"]), splits[iSplit], seed, stepSize)
      if(testIndPred == 1) next
      if(testIndPred != nrow(orderedPred)){
        ind1 <- testIndPred -1
      }else{
        if(orderedPred$predicted[testIndPred] >= splits[iSplit]){
          ind1 <- testIndPred -1
        }else{
          ind1 <- testIndPred
          stopBool <- TRUE
        }
      }
    }
    dataTemp = orderedPred[(ind0+1):ind1,]
    indObs0 <- 0
    orderedObs <- sort(dataTemp$observed)
    for(jSplit in 1:(nSplits+1)){
      seedObs <- max(round(jSplit/length(splits)*length(orderedObs)),1)
      if(jSplit == (nSplits + 1)){
        indObs <- length(orderedObs)
      }else{
        testInd <- findFirstPos(orderedObs, splits[jSplit], seedObs, stepSize)
        if(testInd == 1) next
        if(testInd != length(orderedObs)){
          indObs <- testInd-1 
        }else{
          if(orderedObs[testInd]>= splits[jSplit]){
            indObs <- testInd-1 
          }else{
            indObs <- testInd
            groupMat[jSplit, iSplit] <- indObs - indObs0
            break
          }
        }
      }
      groupMat[jSplit, iSplit] <- indObs - indObs0
      indObs0 <- indObs
    }
    ind0 <- ind1
    if(stopBool) break
  }
  return(groupMat)
}
