#library(rgdal)
#library(stringi)
#matchedTree <- cleanedData$CleanTree
#patchMapTotal <- readOGR(dsn = path.expand("C://Users/jmerwin/Documents/LoryProject/PatchMaps"), layer = "PatchNamedPatchMapFeb7")
#backgroundMap <- readOGR(dsn = path.expand("C://Users/jmerwin/Documents/LoryProject/PatchMaps"), layer = "patches1")
#wingMap <- readOGR(dsn = path.expand("C://Users/jmerwin/Documents/LoryProject"), layer = "wingLine")
#shp <- patchMapTotal

patchSignalMap <- function(cleanDataSet, pcVal = 1, cutoff = 4, retTab = F, retAICTab = F)
{
   library(geiger)
   library(arbutus)
   library(phytools)
    bigData <- NULL
    bigDataReturn <- NULL
    tipLabeledData <- as.data.frame(cleanDataSet$RawDat, stringsAsFactors = T)
    pcaData <- prcomp(cleanDataSet$ColSpace[ , c(20:23)], center = T, scale. = T)
    ##bigMin <- min(pcaData$x[ , pcVal])
    ##bigMax <- max(pcaData$x[ , pcVal])
    fullDataTable <- data.frame()
                    for (i in unique(cleanDataSet$RawDat$Patch))
                    {
                      aicListTMP <- data.frame(stringsAsFactors = F)
                      if(stri_length(i) > 2)
                      {
                        pcaScaled <- pcaData$x[which(tipLabeledData$Patch == i), pcVal]
                        bigMin <- min(pcaScaled)
                        bigMax <- max(pcaScaled)
                        pcaList <- as.numeric(scale.default(x = pcaScaled, center = bigMin, scale = (bigMax - bigMin)))
                        
                        ##pcData <- cleanDataSet$ColSpace[which(tipLabeledData$Patch == i), c(20:23)]
                        ##pcaComp<- prcomp(pcData, center = T, scale. = T)
                        
                        ##pcaList <- pcaComp$x[ , pcVal]
                        nameList <- cleanDataSet$RawDat[which(tipLabeledData$Patch == i), 1]   
                        names(pcaList) <- nameList
                        matchData <- match.phylo.data(data = pcaList, cleanDataSet$CleanTree)
                        #modelChoice <- c("white", "delta", "OU", "BM")
                        modelChoice <- c("OU", "white")
                        fitMax <- 0
                        bestMod <- NULL
                        bestScore <- 0
                        resultsBig <- data.frame(stringsAsFactors = F)
                        thisPhy <- force.ultrametric(matchData$phy, method = "extend")
                        for(k in modelChoice)
                        {
                          scoreNum <- 0
                          
                          pc1Sig = try(fitContinuous(phy = force.ultrametric(thisPhy), dat = matchData$data, model = k))
                          if(any(class(pc1Sig) %in% 'try-error'))
                          {
                            cat('ERROR: fitting model again \n')
                            pc1Sig <- try(fitContinuous(phy= force.ultrametric(thisPhy), dat=matchData$data, model=k))
                            print("Error, trying again")
                          }
                          
                          arbResult = try(arbutus(pc1Sig, nsim = 1000))
                          if(any(class(arbResult) %in% 'try-error')){ #if error, just try it again
                            cat('	ERROR: fitting arbutus again\n')
                            pc1Sig = try(fitContinuous(phy= thisPhy, dat=matchData$data, model=k))		
                            arbResult = try(arbutus(pc1Sig, nsim=1000))
                            print("check")}
                          
                          thisPic <- try(compare_pic_stat(arbResult$obs, arbResult$sim))
                          thisMahal <- (try(mahalanobis_arbutus(thisPic)))
                          thisModFit <- pc1Sig$opt$aic
                          
                          print(thisModFit)
                          
                          print(pc1Sig$opt)
                          
                          
                          if((is.null(bestMod)) && (typeof(thisModFit) != "character"))
                          {
                            fitMax <- thisModFit
                            bestMod <- k
                            try(for(l in arbResult$p.values)
                              {if (!is.null(l) && (l > .05)){scoreNum <- scoreNum + 1}})
                            bestScore <- scoreNum
                            bestArb <- arbResult
                          
                          }
                          
                          bestPic <- try(compare_pic_stat(bestArb$obs, bestArb$sim))
                          bestMahal <- (try(mahalanobis_arbutus(bestPic)))
                          
                          
                          if(((!is.null(bestMod)) && (typeof(thisModFit) != "character")) && (thisModFit < fitMax))
                          {
                            if(((thisModFit) >= (fitMax-2)) & (thisModFit<=(fitMax+2)) && (k == "OU" || k == "delta") && bestMahal < thisMahal){ print(paste("best:", bestMahal, " ", k, thisMahal))
                            
                              
                            } else {fitMax <- thisModFit; bestMod <- k;
                            try(for(l in arbResult$p.values)
                              {if(!is.null(l)){if(l > .05){scoreNum <- scoreNum + 1}}})
                            bestScore <- scoreNum
                            bestArb<-arbResult
                            bestLik <- pc1Sig$opt$lnL}}
                            
                          
                          
                          
                          if(typeof(thisModFit) != "character"){print(paste(k, " is the current model with fit = ", thisModFit, " for ", i, " BestMod is", bestMod))}
                          
                          ##for(j in c(1:6)){if(arbResult$p.values[j] < .05 || is.na(arbResult$p.values[j])){fitScore <- fitScore - 1}}
                          
                        }
                        
                        thisRes <- cbind(i, bestLik, bestMod, bestArb$p.values[1], bestArb$p.values[2], bestArb$p.values[3], bestArb$p.values[4], bestArb$p.values[5], bestArb$p.values[6])
                        fullDataTable <- rbind(fullDataTable, thisRes)
                       
                        
                        if(is.null(bigData))
                            {
                              bigData <- cbind(i, fitMax, bestMod)
                              colnames(bigData) <- c("patchName", "fitMaxMahal", "bestMod")
                              print(bigData)
                        }
                            if(bestScore < cutoff){fitMax <- 0}
                            measure <- cbind(i, fitMax, bestMod)
                            bigData <- rbind(bigData, measure)}
                        
                        
                        
                          #if (fitScore > fitMax)
                          #{
                          #  fitMax <- fitScore
                          #bestMod <- k
                          # print(paste("fitscore was bigger. Best mod is: ", bestMod))
                          #}
                        
                        
                        #if(fitScore > 1)
                        #{
                        #  if(is.null(bigData))
                        #  {
                        #    bigData <- cbind(i, fitScore, bestMod)
                        #    colnames(bigData <- c("patchName", "SigK", "bestMod"))
                        #    print(bigData)
                        #  }
                        #  else{measure <- cbind(i, fitScore, bestMod)
                        #  bigData <- rbind(bigData, measure)}
                        #  print(bigData)
                        #}
                          
                        #if(fitScore < 1)
                        #{
                        #  measure <- cbind(i, 0, "na")
                        #  bigData <- rbind(bigData, measure)
                        #  print(bigData)
                        #}
                        
                        
                        
                    } 
    print(bigData)
    colnames(fullDataTable) <- c("Patch", "bestLik", "bestMod", "m.sig", "c.var", "s.var", "s.asr", "s.hgt", "d.cdf")
    rownames(fullDataTable) <- fullDataTable$Patch
    if(retAICTab)
    {colnames(aicList) <- c("OU", "BM", "White.Noise")
    return(aicList)}
    if(retTab == T)
      {return(fullDataTable)}
    return(bigData)
}
   




                   
                    #  colConv <- cbind(c("BM", "OU"), c(1, 2))
#                    cols <- mapvalues(bigData, from = colConv[,1], to = colConv[,2])
                    
                    ##colVals <- hsv(s = (as.numeric(bigDataTib$SigK)/6), h = cols/3, alpha =.7)
                    ##colVals <- hsv(h = (as.numeric(bigDataTib$SigK)), s = .5, alpha =.8 )
 #                   bigDataCols <- cbind(bigDataTib$patchName, colVals)
#                    colnames(bigDataCols) <- c("patch", "colVals", "colMod")
 #                   colConv <- cbind(modelChoice, c(1:length(modelChoice)))
  #                  cols <- mapvalues(bigData, from = colConv[,1], to = colConv[,2])
   #                 colnames(bigDataCols) <- c("patch", "colVals")
    #                bigDataJoin <- as.data.frame(stringsAsFactors = F, bigDataCols)
     #               joinedSig <- sp::merge(x = patchMapTotal, bigDataJoin, by.x = "PatchName", by.y = "patch", duplicateGeoms = TRUE)
      #              colsSig <- joinedSig@data$colVals
       #             plot(backgroundMap, col = "#d7d9d nd")
        #            plot(add = T, joinedSig, col = colsSig)
                   























                    
                    