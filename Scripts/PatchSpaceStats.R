##BirdSpaceSignal

birdPatchSpaceStats <- function(cleanedDataSet, patch = NULL)
{
  flag <- F
  patchList <- data.frame(stringsAsFactors = F)
  if(!is.null(patch))
  {
    thisSpaceCol <- na.omit(cleanedDataSet$RawDat[which(cleanedDataSet$RawDat$Patch == patch), ])
    thisSpacecolNoDupes <- thisSpaceCol[!duplicated(thisSpaceCol$mappedVals),]
    
    rownames(thisSpacecolNoDupes) <- thisSpacecolNoDupes$mappedVals
    datSum <- thisSpacecolNoDupes
    return(datSum)
    
  } else {for(i in unique(cleanedDataSet$RawDat$mappedVals))
  {
  matchLocList <- grep(pattern = i, ignore.case = T, value = F, x = cleanedData$RawDat$mappedVals)
  thisSpaceCol <- na.omit(cleanedDataSet$ColSpace[matchLocList, c(20:23)])
  tempSpace <- colspace(qcatch = "Qi", vismodeldata = thisSpaceCol, space = "tcs")
  datSum <- try(summary((tempSpace)))
  patchList <- rbind(patchList, datSum)
  }
  pcaDat<- prcomp(patchList, retx = T, scale. = T)
  rownames(patchList) <- (unique(cleanedDataSet$RawDat$mappedVals))
  patchListFinal <- cbind(pcaDat$x[ , 1], patchList)
  return(patchList)
  }
  }





