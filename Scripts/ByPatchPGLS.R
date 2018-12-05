runPGLSanalysis <- function(colPCSelect = 1, brightnessOn = F, bodySize = F, retHist = F, plotHist = F)
{
  bigDatPCA <- prcomp(cleanedData$ColSpace[ ,c(20:23)], center = T, scale. = T, retx = T)
  rownames(bigDatPCA$x) <- make.names(unique = T, names = cleanedData$RawDat$mappedVals, allow_ = T)
  
  
  
byPatchPGLS <- function(patches = NULL, colpcVal = 1, cleanedDataset = cleanedData, climData, brightnessSelect = brightnessOn)
{
  
  pglsModel2 <- NULL
  matchCrown <- NULL
  
  
  patchCols <- cleanedDataset$ColSpace[which(cleanedDataset$ColSpace$Patch == patches), c(20:23)]
  rownames(patchCols) <- make.names(unique = T, names = cleanedDataset$RawDat$mappedVals[which(cleanedDataset$ColSpace$Patch == "Wrist")])
  
  pcaClim <- prcomp((climData[, -c(1)]^10), center = T, scale. = T)
  print(summary(pcaClim))
  matchPGLS <- match.phylo.data(cleanedData$CleanTree, data = pcaClim$x)
  print("ClimTree Matched")
  if(!brightnessSelect)
  {
  crownPC <- as.numeric(scale.default(bigDatPCA$x[which(cleanedDataset$ColSpace$Patch == patches) , colpcVal], center = T, scale = T))
  }
  names(crownPC) <- make.names(unique = T, names = cleanedData$RawDat$mappedVals[which(cleanedDataset$ColSpace$Patch == patches)], allow_ = T)
  
  
  if(brightnessSelect)
  {
    
    patchCols = cleanedData$ColSpace$lumMean[which(cleanedData$ColSpace$Patch == patches)]
    
    names(patchCols) <- make.names(unique = T, names = cleanedData$RawDat$mappedVals[which(cleanedData$ColSpace$Patch == patches)], allow_ = T)
    crownPC <- scale.default(patchCols, center = min(patchCols), scale = max(patchCols) - min(patchCols))
    names(crownPC) <- names(patchCols)
    crownPC <- as.matrix(crownPC)
    matchCrown <- match.phylo.data(phy = matchPGLS$phy, data = (crownPC[ , 1]))
  
    
 
  } else {matchCrown <- match.phylo.data(matchPGLS$phy, crownPC)}
  
  if(bodySize)
  {
    matchCrown <- match.phylo.data(matchPGLS$phy, )
    loryMorphoSheet <- read.csv(file = "~/Lory MorphoMetrics - MorphoSheetForTipMatch.csv", header = T, stringsAsFactors = F)
    cleanDataSetMorpho <- data.frame(stringsAsFactors = F, loryMorphoSheet[ , c(8:10)])
    cleanDataSetMorphoPCA <- prcomp(cleanDataSetMorpho)
    dataforMorpho <- lapply(FUN = scale.default(), X = cleanDataSetMorphoPCA$x, center = min(X), scale = max(X) - min(X))
    matchedMorpho <- match.phylo.data(phy = cleanedData$CleanTree, data = cleanDataSetMorphoPCA$x)
  
  }
  
  
  
  
  #hist(crownPC, main = paste(patches, skewness(crownPC)))
  
  colorDataSet <- cbind(matchCrown$data, matchPGLS$data)
  
  colnames(colorDataSet)[1] <- "V1" 

  print("melded")
  colorDataSetFrame <- as.data.frame(colorDataSet)
  colorDataSetFrame[ , 2] <- colorDataSetFrame[ , 2]
  
  colorDataSetFrameMerge <- cbind(colorDataSetFrame, names(matchCrown$data))
  caperCompare <- caper::comparative.data(phy = nameNodes(tree = force.ultrametric(cleanedData$CleanTree)), data = colorDataSetFrameMerge, names.col = "names(matchCrown$data)")
  ##pglsModel2 <- gls(V1 ~ PC1*PC2*PC3, correlation = corMartins(value = 1, phy = matchPGLS$phy), data = colorDataSetFrame, method = "ML")
  caperModel <- caper::pgls(formula = V1 ~ PC1, data = na.omit(caperCompare))
  
  print("pglsRun")
  return(caperModel)
}


setwd(here::here("Output"))
if(retHist)
{
  returnTab <- list()
}

for(j in c(1:35))
{
  i <- unique(cleanedData$RawDat$Patch)[j]
  print(i)
  
  name <- path.expand(here::here("Output", "PGLS", paste(i, colPCSelect, ".png", sep = "")))
  if(!retHist){
  png(filename = name, res = 200, width = 10, height = 5, units = "in", bg = "white")
  }
  #par(mfrow = c(1, 2))
  thisPGLS <- byPatchPGLS(patch = i, cleanedDataset = cleanedData, climData = climDataImport, colpcVal = colPCSelect)
  
  
  summ <- caper::summary.pgls(thisPGLS)
  
  if(!retHist)
  {
  plot(thisPGLS$data$data$V1~thisPGLS$data$data$PC1, main = paste(i, summ$coefficients[2, 4]), sub = paste( summ$adj.r.squared))
  text(y = thisPGLS$data$data$V1+.5, x = thisPGLS$data$data$PC1+.4, cex = .2)
  flag <- F
  for(i in summ$coefficients[ , 4]){if(i<.001428){flag <- T}}
  red1 <- 0
  if(flag){red1 = 1}
  abline(thisPGLS$model$coef[c(1,2)], col = rgb(red = red1,blue = 0, green = 0)) 
  #dev.off()
  }
  if(retHist)
  {
    if(!is.null(returnTab)){returnTab <- c(returnTab, summ$adj.r.squared)
    print(returnTab)
    } else {
    returnTab <- summ$adj.r.squared}
    print(returnTab)
  }
}
if(plotHist){return(hist(as.numeric(paste(returnTab))))}
if(retHist && !plotHist){(return(as.numeric(paste(returnTab))))}


}

