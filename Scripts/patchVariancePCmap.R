#PatchVariance
makeModelFitMap <- function(fileDestination = here::here("Output/ModelMaps"), setModel = NULL, setModelPC = 1, modelOptionSet = T)
{
byPatchPCVariance <- function(cleanedDataSet = cleanedData, raw = F, modelOption = modelOptionSet, modelPC = setModelPC, modelPick = setModel)
{
  if(!raw)
  {
  fullPCAData <- prcomp((cleanedDataSet$RawDat[ , c("u", "s", "m", "l")]), center = T, scale. = T)
  patchMeld <- cbind(cleanedDataSet$RawDat, fullPCAData$x[, setModelPC])
  
  }
  if(raw)
  {
    fullPCAData <- prcomp((cleanedDataSet$RawDat[ , c("uvMean", "swMean", "mwMean", "lwMean")]))
  }
  patchMeld <- cbind(cleanedDataSet$RawDat, fullPCAData$x[, c(1, 2)])
  patchList <- data.frame(stringsAsFactors = F)
  datSum <- data.frame(stringsAsFactors = F)
  for(i in unique(cleanedDataSet$RawDat$Patch))
  {
    matNamesBoundFrame <- NULL
    nameBoundPC <- NULL
    matchLocList <- unique(grep(pattern = i, ignore.case = T, value = F, x = cleanedData$RawDat$Patch))
    if(!raw)
    {
     fullPCAData1 <- as.numeric(patchMeld[matchLocList, 24+modelPC])
     print(fullPCAData1)
     fullPCADataX <- as.numeric(scale.default(fullPCAData1, center = min(fullPCAData1), scale = (max(fullPCAData1)-min(fullPCAData1))))

     print(i)
    }
    if(raw)
    {
      ##fullPCAData <- prcomp((cleanedDataSet$RawDat[matchLocList , c("uvMean", "swMean", "mwMean", "lwMean")]))
    }
    
    
    
    ##thisSpaceCol <- na.omit(fullPC[matchLocList, c(25, 26)])
    
    
    
    
    if(modelOption)
    {
    thisSpaceCol <- data.frame(fullPCADataX)
    print("1")
    print(thisSpaceCol)
    print(kurtosis(thisSpaceCol$fullPCADataX))
    print(i)
    nameBoundPC <- cbind(cleanedDataSet$RawDat$mappedVals[matchLocList], thisSpaceCol)
    print("2")
    these.names = make.names(nameBoundPC[, 1], unique = T, allow_ = T)
    matNamesBoundFrame <-  as.numeric(nameBoundPC[ , 2])
    names(matNamesBoundFrame) <- these.names
    print(head(matNamesBoundFrame))
    thisMatch <- match.phylo.data(phy = cleanedData$CleanTree, matNamesBoundFrame)
    print("3")
    ultra <- force.ultrametric(thisMatch$phy)
    dataForOU <- as.numeric(paste(thisMatch$data))
    print(dataForOU)
    names(dataForOU) <- names(thisMatch$data)
    print(dataForOU)
    thisFit <- fitContinuous(ultra, dataForOU, model = modelPick)
    print("6")
    if(modelPick == "delta"){datSum <- cbind(i, thisFit$opt$delta)}
    if(modelPick == "BM"){datSum <- cbind(i, thisFit$opt$sigsq)}
    if(modelPick == "OU"){datSum <- cbind(i, thisFit$opt$alpha)}
    if(modelPick == "lambda"){datSum <- cbind(i, thisFit$opt$lambda)}
    print("7")
    
    } else {
    varPC1 <- var(fullPCADataX)
    datSum <- cbind(i, varPC1)
    }
    patchList <- rbind(patchList, datSum)
  }
  
  return(patchList)
  
}

if(modelOptionSet == T){pcDatCols3 <- (byPatchPCVariance(cleanedDataSet = cleanedData, raw = F, modelOption = modelOptionSet, modelPC = setModelPC, modelPick = setModel))
} else {pcDatCols3 <- (byPatchPCVariance(modelOption = F, modelPC = setModelPC))}

##pcDatCols3 <- matrix(testcols)
if(modelOptionSet == F) {
  pcDatCols2<- (pcDatCols3)
} else {pcDatCols2<- pcDatCols3}

##hist(pcDatCols)

##pcDatCols <- data.frame(stringsAsFactors = F, unique(cleanedData$ColSpace$Patch), pcDatCols2)
colnames(pcDatCols2) <- c("i", "varPC1")
joined <- merge(x = patchMapTotal, pcDatCols2, by.x = "PatchName", by.y =  "i")

print(head(joined@data))

##plots patchmaps and writes to file
setwd(fileDestination)
namePlot = paste(setModel, "pc", setModelPC, ".png", sep = "")
##png(filename= namePlot, units = "in", height = 5, width = 5, res = 300, bg = "transparent")

plot(backgroundMap, border = "black", col = "white")
##coList <- rgb(maxColorValue = 1, green = 1-(.3+(as.numeric(paste(joined@data$pcDatCols)))), red = .2, blue =  .2)

if(!is.null(setModel))
{
if(setModel == "lambda")
{
coList <- rgb(maxColorValue = 1, red = rescale(to = c(0, .7), (log(1+as.numeric(paste(joined@data$varPC1))))), blue = rescale(to = c(0, .7), (log(1+as.numeric(paste(joined@data$varPC1))))), green =  .1)
}
if(setModel == "delta")
{
  coList <- rgb(maxColorValue = 5, blue = 5-rescale(to = c(2, 4), ((as.numeric(paste(joined@data$varPC1))))), green = 4-rescale(to = c(2, 3), ((as.numeric(paste(joined@data$varPC1))))), red =  1)
}
if(setModel == "OU")
{
  coList <- rgb(maxColorValue = 3, green = rescale(to = c(0, 2.5), 1.5+((log(2)/as.numeric(paste(joined@data$varPC1))))), red = 0, blue =  (as.numeric(paste(joined@data$varPC1))))
  print(coList)
}
if(setModel == "BM")
{
  coList <- rgb(maxColorValue = 1, blue = .1, red = rescale(to = c(.2, .8),((as.numeric(paste(joined@data$varPC1))))), green = .15)

  }
} else {
    coList <- rgb(maxColorValue = 1, blue = .1, red = rescale(to = c(.2, .8),((as.numeric(paste(joined@data$varPC1))))), green = rescale(to = c(.2, .6),((as.numeric(paste(joined@data$varPC1))))))
}





plot(add = T, joined, col = coList, border = coList)
  plot(wingMap, add = T, col = alpha("gray", .3))

#legend(cex = .5, title = "Color Variance by Patch",   # location of legend
       #legend = c("Maximum Variance", "Minimum Variance"), # categories or elements to render in
       # the legend
       #fill = range(joined@data$varPC1))

}




#patternDat <- read.delim("C://Users/jmerwin/Downloads/Lory Pattern Data - Side View Full Bird Pattern Data.tsv", header = T, sep = "\t")
#patternDat <- patternDat[-c(1) , ]
#rownames(patternDat) <- patternDat$Species
#patternPC <- prcomp(patternDat[, -c(1)])
#scatterplot3d::scatterplot3d(patternPC$x[ , c(1:3)])




