

relativeFitTable <- function(cleanDataSet = cleanedData, pcVal = 1, models = c("OU", "BM", "delta", "white"), plot = T)
{
  aicList <- data.frame(stringsAsFactors = F)
  
  tipLabeledData <- as.data.frame(cleanDataSet$RawDat, stringsAsFactors = T)
  pcaData <- prcomp(cleanDataSet$ColSpace[, c(20:23)], center = T, scale. = T)
  
  firstFlag <- T
  for (i in unique(cleanDataSet$RawDat$Patch))
  {
    aicListTMP <- data.frame(stringsAsFactors = F)
    firstFlag = F
    if(stri_length(i) > 2)
    
    {
      
      pcaDataX <- pcaData$x[which(tipLabeledData$Patch == i) , pcVal]
      pcaScaled <- as.numeric(scale.default(pcaDataX, center = min(pcaDataX), scale = (max(pcaDataX) - min(pcaDataX))))
      pcaList <- pcaScaled
      nameList <- cleanDataSet$RawDat[which(tipLabeledData$Patch == i), 1]   
      names(pcaList) <- nameList
      matchData <- match.phylo.data(data = pcaList, cleanDataSet$CleanTree)
      modelChoice <- models
      
      
      
      for(k in modelChoice)
      {
        thisFit <- try(fitContinuous(dat = matchData$data, phy = force.ultrametric(matchData$phy), model = k,  bounds = c(0, 100)))
        if(any(class(thisFit) %in% 'try-error'))
        {
          cat('ERROR: fitting model again \n')
          pc1Sig <- try(fitContinuous(phy= thisPhy, dat=matchData$data, model=k, bounds = c(0, 100)))
          print("Error, trying again")
        }
        if(k == modelChoice[1]){aicListTMP <- thisFit$opt$aicc} else {aicListTMP <- cbind(aicListTMP, thisFit$opt$aicc)}
      }
      if(firstFlag){aicList <- aicListTMP} else {aicList <- rbind(aicList, aicListTMP)}
      
    }
  }
  ##rownames(aicList) <- patchesTested
  ##colnames(aicList) <- models
  weightTable <- NULL
  weightTable <- data.frame()
  print(aicList)
  for(j in c(1:length(aicList[ , 1])))
  {
    print(j)
    thisWeight <- aicw(as.numeric(unlist(aicList[j, ])))
    print(thisWeight)
    if(length(weightTable) != 0){weightTable <- cbind(weightTable, thisWeight$w)
    } else {weightTable <- thisWeight$w}
  
  }
  colnames(weightTable) <- unique(cleanDataSet$RawDat$Patch)
  print(weightTable)
  rownames(weightTable) <- models
  print("2")
  
  if(plot == T)
  {
  par(mar = c(9, 3, 3, 1))
  barplot(weightTable, col = c("coral4", "darkgoldenrod3", "skyblue3", "white"), #legend.text = c(rownames(weightTable)), 
          main = "Relative AICc Weights by Patch", pch = .4, cex.names = .8, las = 2, space = .5)
  } else {return(weightTable)}
}