##testing model clustering
longColListAllPatch <- matrix()
pcList <- prcomp(cleanedData$RawDat[ , c(21:24)], center = T, scale = T)

flag = F
for(i in unique(cleanedData$RawDat$mappedVals))
{
  #grab35patchesfrom1bird
  fullPatchPCs <- list()
  matchList <- pcList$x[which(cleanedData$RawDat$mappedVals == i), ]
  thisBird <- t(matchList)
  thisBird2 <- c(thisBird[2, c(1:35)])
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



force.ultrametric(cleanedData$CleanTree)
#match that dataframe with the ultrametric tree
matchedAllPatch <- match.phylo.data(data = longColListAllPatch, phy = force.ultrametric(cleanedData$CleanTree))
##matchedAllPatchPC <- match.phylo.data(data = testVertPC$x, phy = force.ultrametric(cleanedData$CleanTree))
#run ratematrix on that data
heatmap(x = matchedAllPatch$data)
(ratematrix(matchedAllPatch$phy, matchedAllPatch$data))
heatmap(ratematrix(matchedAllPatch$phy, matchedAllPatch$data))
contMap(x = matchedAllPatchPC$data[, 1], tree = matchedAllPatchPC$phy)

##test all under OU 
listVect <- split(matchedAllPatch$data[, c(1:5)], c(col(matchedAllPatch$data[ , c(1:5)])))
for(i in c(1:5)){names(listVect[[i]]) <- rownames(matchedAllPatch$data)}
thisTree <- phytools::ratebytree(trees = as.multiPhylo(matchedAllPatch$phy), x = listVect, type = "continuous", model = "BM")


##test all under BM
sumOpts <- function(optObs = NULL)
{
  thisSum <- 0
  for(i in c(1:35))
  {
    thisSum <- (thisSum + bmOpt[[i]]$opt$lnL)
  }
  return(thisSum)
}


##TestBestMods
bmOpt <- fitContinuous(control=list(niter=50), SE = NA, phy = geiger::rescale(x = force.ultrametric(matchedAllPatch$phy), model = "delta", delta = 3), dat = matchedAllPatch$data, model = "BM")
print(sumOpts(bmOpt))

library(phylocurve)

patchOU <- paste(signalTabPC1$Patch[which(signalTabPC1$bestMod == "OU")])
evoOU <- evo.model(model = "OU", tree = force.ultrametric(matchedAllPatch$phy), Y=  matchedAllPatch$data[, patchOU])

patchBM <- paste(signalTabPC1$Patch[which(signalTabPC1$bestMod == "BM")])
evoBM <- evo.model(model = "BM", tree = force.ultrametric(matchedAllPatch$phy), Y = matchedAllPatch$data[, patchBM])

patchDelta <- paste(signalTabPC1$Patch[which(signalTabPC1$bestMod == "delta")])
evoDelta <- evo.model(model = "Delta", tree = force.ultrametric(matchedAllPatch$phy), Y =  matchedAllPatch$data[, patchDelta])
sum(evoOU$logL, evoBM$logL, evoDelta$logL)


evoAll <- evo.model(multirate = T, ret.level = 3, model = "BM", tree = force.ultrametric(matchedAllPatch$phy), Y = matchedAllPatch$data)
evoAll$logL
corrplot(is.corr = F, evoAll$phylocov)

evoAll$model.par

corrplot(cor(matchedAllPatch$data))


