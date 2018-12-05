#RateMatrixTest


?ratematrix





#combine PC1 accross all patches for each bird, running scaling on each patch to control for variance
  #for each bird, grab all patches
  longColListAllPatch <- matrix()
  pcList <- prcomp(cleanedData$RawDat[ , c(21:24)], center = T, scale = T)
  
  flag = F
  for(i in unique(cleanedData$RawDat$mappedVals))
  {
      #grab35patchesfrom1bird
      fullPatchPCs <- list()
      matchList <- pcList$x[which(cleanedData$RawDat$mappedVals == i), ]
      thisBird <- t(matchList)
      thisBird2 <- thisBird[1, c(1:35)]
      thisBirdScaled <- scale.default(thisBird2, center = min(thisBird2), scale = max(thisBird2)-min(thisBird2))
      
      #bindbirdsintodataframe
      if(!flag)
      {
        longColListAllPatch <- (as.numeric(thisBirdScaled)) ; flag = T;
      } else {
       longColListAllPatch <- rbind(longColListAllPatch, (as.numeric(thisBirdScaled)))
      }
      
  }
  colnames(longColListAllPatch) <- unique(cleanedData$RawDat$Patch)
  rownames(longColListAllPatch) <- unique(cleanedData$RawDat$mappedVals)



  force.ultrametric(cleanedData$CleanTree)
  
  rownames(testVertPC$x) <- unique(cleanedData$RawDat$mappedVals)
#match that dataframe with the ultrametric tree
  matchedAllPatch <- match.phylo.data(data = longColListAllPatch, phy = force.ultrametric(cleanedData$CleanTree))
  matchedAllPatchPC <- match.phylo.data(data = testVertPC$x, phy = force.ultrametric(cleanedData$CleanTree))
#run ratematrix on that data
  heatmap(x = matchedAllPatch$data)
  (ratematrix(matchedAllPatch$phy, matchedAllPatch$data))
  heatmap(ratematrix(matchedAllPatch$phy, matchedAllPatch$data))
  contMap(x = matchedAllPatchPC$data[, 1], tree = matchedAllPatchPC$phy)
  
  evo.model(tree = matchedAllPatch$phy, Y = matchedAllPatch$data, plot.LL.surface = T)
  
  corrplot(is.corr = F, (matchedAllPatchPC$data))
  thisCor <- cor.mtest(matchedAllPatch$data)
  pdf(file = here::here("CorrPlotSigLevel.pdf"), width = 7, height = 7)
  corrplot(rect.col = "black", bg = "white", insig = "blank", sig.level = 0.001428571, p.mat = thisCor$p, cor(matchedAllPatch$data), method = "color", title = "Between-Patch Correlation", mar = c(1, 1, 1, 1), tl.col = "black", order = "hclust", addrect = 3)
  dev.off()
  
  
  
  
  compare.models(evoFitOUMULT, evoFitOUStat, plot = T)
  
  
  
  evoFitOUMULT <- evo.model(multirate = T, ret.level = 3, model = "OU", tree = matchedAllPatch$phy, Y = matchedAllPatch$data, plot.LL.surface = F)
  evoFitOUStat <- evo.model(multirate = F, ret.level = 3, model = "OU", tree = matchedAllPatch$phy, Y = matchedAllPatch$data, plot.LL.surface = F)
  
  evoFitBM$logL
  corrplot(evoFitBM$phylocov, method = "color", is.corr = F)
  
  ancEst <- ultraFastAnc(phy = matchedAllPatch$phy, x = matchedAllPatch$data, CI = T)
  paintBranches(tree = matchedAllPatchancEst)
  heatmap(evoFitBM$phylocov)
  evoFitBM$logL
  K.mult(evoFitBM)
  evoFitOU <- evo.model(model = "delta", tree = matchedAllPatch$phy, Y = matchedAllPatch$data[, c("Lores", "AuricularFeathers", "SideBreast", "Collar", "UpperBreast", "LowerBreast", "MalarRegion", "ThroatSides")])
  evoFitOU$logL
  K.mult(evoFitOU)
  evoFitDelta <- evo.model(model = "delta", tree = matchedAllPatch$phy, Y = matchedAllPatch$data[ ,  c(10:14)], plot.LL.surface = F)
  evoFitDelta$logL
  K.mult(evoFitDelta)
  evoFitDelta <- evo.model(model = "lambda", tree = matchedAllPatch$phy, Y = matchedAllPatch$data[ ,  pList], plot.LL.surface = F)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  