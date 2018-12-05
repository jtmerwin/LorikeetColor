##BirdSpaceModelCompare
fullBirdSpaceStats <- function(cleanedDataSet)
{
  flag <- F
  patchList <- data.frame(stringsAsFactors = F)
  thisSpaceCol <- cleanedDataSet$ColSpace[, c(20:23)]
  tempSpace <- colspace(qcatch = "Qi", vismodeldata = thisSpaceCol, space = "tcs")
  datSum <- summary(tempSpace, by = cleanedDataSet$RawDat$mappedVals)
  pcComp <- prcomp(datSum[ , c(1:4)], center = F, scale. = F)
  plot(colspace(datSum[ , c(1:4)], space = "tcs", qcatch = "Qi"), col = rgb(red = datSum[,4], blue = datSum[,2], green = datSum[,3]))
  colValList <- rgb(red = datSum[,4], blue = datSum[,2], green = datSum[,3])
  matchup <- na.omit(match.phylo.data(cleanedDataSet$CleanTree, datSum))
  
  #png(filename = "~/Output4Contmaps", width = 5, height = 5, units = "in", bg = "lightgray", res = 300)
  #par(mfrow = c(2, 2))
  ##for(i in c(1:4))
  #{
  #contMap(res = 1000, matchup$phy, log(.6+matchup$data[ , i]),ftype = "reg", fsize = .25, lwd = 2, outline = F)
  #OUphenogram
  #}
  
  #dev.off()
  
  
  
  matchup <- na.omit(match.phylo.data(cleanedDataSet$CleanTree, datSum))
  nameData <- matchup$data[, 6]
  names(nameData) <- rownames(matchup$data)
  #nameDat <- (prcomp(na.omit(matchup$data)))
  matchup <- match.phylo.data(cleanedDataSet$CleanTree, datSum)
  modelChoice <- c("BM", "OU", "white", "EB")
  fitMax <- 0
  bestMod <- NULL
  bestScore <- 0
  resultsBig <- data.frame(stringsAsFactors = F)
  
  
  
  for(k in modelChoice)
  {
    scoreNum <- 0
    pc1Sig = try(fitContinuous(phy = matchup$phy, dat = matchup$data, model = k))
    if(any(class(pc1Sig$PC1) %in% 'try-error'))
    {
      cat('ERROR: fitting model again \n')
      pc1Sig <- try(fitContinuous(phy= matchup$phy, dat=matchup$data, model=k))
      print("Error, trying again")
    }
    
    arbResult = try(arbutus(pc1Sig$PC1, nsim = 100))
    if(any(class(arbResult) %in% 'try-error')){ #if error, just try it again
      cat('	ERROR: fitting arbutus again\n')
      pc1Sig = try(fitContinuous(phy= matchup$phy, dat=matchup$data, model=k))		
      arbResult = try(arbutus(pc1Sig, nsim=100))
      print("check")}
    
    thisPic <- try(compare_pic_stat(arbResult$obs, arbResult$sim))
    thisModFit <- try(mahalanobis_arbutus(thisPic))
    
    print(thisModFit)
    if((is.null(bestMod)) && (typeof(thisModFit) != "character"))
    {
      fitMax <- thisModFit
      bestMod <- k
      for(l in arbResult$p.values)
      {if (!is.null(l) && l > .05){scoreNum <- scoreNum + 1}}
      bestScore <- scoreNum
    }
    
    if(((!is.null(bestMod)) && (typeof(thisModFit) != "character")) && (thisModFit < fitMax))
    {
      fitMax <- thisModFit
      bestMod <- k
      for(l in arbResult$p.values)
      {if (!is.null(l) && l > .05){scoreNum <- scoreNum + 1}}
      bestScore <- scoreNum
      
    }
    
    if(typeof(thisModFit) != "character"){print(paste(bestMod, " is the current best model with fit = ", thisModFit))}
    
    ##for(j in c(1:6)){if(arbResult$p.values[j] < .05 || is.na(arbResult$p.values[j])){fitScore <- fitScore - 1}}
    
  }
  
  if(is.null(bigData))
  {
    bigData <- cbind(i, fitMax, bestMod)
    colnames(bigData) <- c("patchName", "fitMaxMahal", "bestMod")
    print(bigData)
  }
  if(bestScore < 3){fitMax <- 0}
  measure <- cbind( fitMax, bestMod)
  return(measure) 
}


  
  
  
  #joinedSig <- sp::merge(x = patchMapTotal, patchListFinal, by.x = "PatchName", by.y = "rownames(x)", duplicateGeoms = TRUE)
  
  #plot(backgroundMap, col = "light gray", border = "transparent")
  #plot(border = "transparent", add = T, joinedSig, col = rgb(red = joinedSig@data$centroid.l, green = joinedSig@data$centroid.m, blue = joinedSig@data$centroid.s, alpha = joinedSig$rel.c.vol*4))
  #return(bigData)










#{
#matchup <- na.omit(match.phylo.data(cleanedDataSet$CleanTree, datSum))
#nameData <- matchup$data[, 8]
#names(nameData) <- rownames(matchup$data)
#nameDat <- (prcomp(na.omit(matchup$data)))
#matchup <- match.phylo.data(cleanedDataSet$CleanTree, datSum)

#fitMax <- 0
#bestMod <- NULL
#bestScore <- 0
#resultsBig <- data.frame(stringsAsFactors = F)

#allMods <- list()
#phyMatch <- force.ultrametric(matchup$phy, "extend")
#for(k in c("BM", "OU", "trend"))
#{
 # mod <- fitContinuous(dat = log(nameData), phy = phyMatch, model = k)
#  print(mod$opt$aicc)
#}
#}

#mod <- fitContinuous(dat = log(nameData), phy = phyMatch, model = "OU")
#arb <- arbutus(x =mod, nsim = 1000)
#plot(arb)





