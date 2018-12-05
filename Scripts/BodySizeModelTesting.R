##BodySizeModelTesting

fullBirdSpaceStats <- function(cleanedDataSet)
{
  flag <- F
  bodySize <- as.numeric(cleanedDataSet$MetaDat$Body.Size)
  names(bodySize) <- (cleanedData$MetaDat$TreeLabel)
  matchup <- match.phylo.data(phy = cleanedData$CleanTree, data = bodySize)
  contMap(res = 1000, matchup$phy, abs(matchup$data)-mean(matchup$data),ftype = "reg", fsize = .25, lwd = 2, outline = F)
  modelChoice <- c("BM", "OU", "kappa", "lambda")
  
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
}
  