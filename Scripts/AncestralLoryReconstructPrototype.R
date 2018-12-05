#Find best-fit models for each patch

longColListAllPatch <- matrix()
pcList <- prcomp(cleanedData$RawDat[ , c(21:24)], center = T, scale = T)

ancestralCol <- function(componentnum = NULL)
{
flag = F
for(i in unique(cleanedData$RawDat$mappedVals))
{
  #grab35patchesfrom1bird
  fullPatchPCs <- list()
  matchList <- cleanedData$RawDat[which(cleanedData$RawDat$mappedVals == i), c(20:24)]
  thisBird <- t(matchList)
  thisBird2 <- c(thisBird[componentnum, c(1:35)])
  thisBirdScaled <- scale.default(thisBird2, center = min(thisBird2), scale = max(thisBird2)-min(thisBird2))
  
  #bindbirdsintodataframe
  if(!flag)
  {
    longColListAllPatch <- (as.numeric(thisBirdScaled)) ; flag = T;
  } else {
    longColListAllPatch <- rbind(longColListAllPatch, (as.numeric(thisBirdScaled)))
  }
  
}
colnames(longColListAllPatch) <- (unique(cleanedData$RawDat$Patch))
rownames(longColListAllPatch) <- unique(cleanedData$RawDat$mappedVals)
return(longColListAllPatch)
}

collateMeasure <- function(matchedSet = NULL)
{
  ancList <- NULL
  for(i in c(1:35))
  {
    ancRes <- fastAnc(tree = matchedSet$phy, x = matchedSet$data[, i])
    ancList <- c(ancList, ancRes["100"])
    
  }
  return(ancList)
}
  

reconCol <- function(ancResult = NULL)
{
  matchedAllPatch <- match.phylo.data(data = ancResult, phy = force.ultrametric(cleanedData$CleanTree))
  totalAnc <- collateMeasure(matchedSet = matchedAllPatch)
  return(totalAnc)
}



fullPipeLinefunction <- function(thisNum = NULL)
{
  thisOne <- ancestralCol(componentnum = thisNum)
  thisTwo <- reconCol(ancResult = thisOne)
  return(thisTwo)
}

ancResults <- cbind(fullPipeLinefunction(2), fullPipeLinefunction(3), fullPipeLinefunction(4))
colList <- rgb(blue = scales::rescale(to = c(0, 1), ancResults[ ,1]), green = ancResults[,2], red = ancResults[,3])
matchList <- cbind(unique(cleanedData$RawDat$Patch), (colList))
names(matchList) <- NULL
colnames(matchList) <- c("i", "coList")
matchList <- as.data.frame(matchList)
joined <- merge(x = patchMapTotal, matchList, by.x = "PatchName", by.y =  "i")
plot(backgroundMap, border = "black", col = "white")
#coList <- rgb(maxColorValue = 1, red = rescale(to = c(.1, .6), (log(1+as.numeric(paste(joined@data$varPC1))))), blue = rescale(to = c(.1, .6), (log(1+as.numeric(paste(joined@data$varPC1))))), green =  .1)
plot(add = T, joined, col = joined@data$coList, border = joined@data$coList)
plot(wingMap, add = T, col = alpha("gray", .3))



#Using best fit model, transform tree

#perform fast ANC, find ancestral state of patch

#save to patch data table

#map all patch ancestral states to a bird
