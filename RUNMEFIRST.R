#dependencies and source other files

require("plyr")
require("rgdal")
require("phytools")
require("geiger")
require("scales")
require("stringi")
require("here")
require("ggtree")
require("picante")
require("surface")



patchMapTotal <- readOGR(dsn = path.expand(here("Shapes")), layer = "PatchNamedPatchMapFeb7")
backgroundMap <- readOGR(dsn = path.expand(here("Shapes")), layer = "patches1")
wingMap <- readOGR(dsn = path.expand(here("Shapes")), layer = "wingLine")
shp <- patchMapTotal


#sourceIncludedMethods
##install.packages("BiocInstaller")
##library(BiocInstaller)
##biocLite("ggtree")

for(i in c(paste(sep = "/", here::here("Scripts"), list.files(here::here("Scripts")))))
{
  cat(i, "\n")
  try(source(i))
}

#Step1: run DataCleaner, Output Trimmed Tree and Cleaned Colorspace. Make Patchmaps and Tree
pdf(file = here("Output/TreeAndColspace.pdf"), bg = "white", height = 4, width = 4)
cleanedData <- dataCleaner(rawColDataLoc = here("Data/FullCompleteMasterSheet.csv"), DataKeyLoc = here("Data/FullMatchTreeD70Full.csv"), treeLoc = here("Data/lcr20_all_treePL.tre"))

#generatespatchmaps from data
makePatchMaps(cleanedData, patchMap = shp, numPatches = 35, fileDestination = here::here("Shapes", "TipNamedPatchMaps"))

#outputPatchMap and Tree. Requires ggtree bioconductor install
patchMapTreePlot(cleanedData)
dev.off()




#Step 2: import blank patchmap, output 4 patchmaps for 2 PCs each, 8 pdf outputs, variance, brownian rate, sigma value, arbutus maps
{
patchMapTotal <- readOGR(dsn = path.expand(here("Shapes")), layer = "PatchNamedPatchMapFeb7")
backgroundMap <- readOGR(dsn = path.expand(here("Shapes")), layer = "patches1")
wingMap <- readOGR(dsn = path.expand(here("Shapes")), layer = "wingLine")
shp <- patchMapTotal

signalMapPC1 <- patchSignalMap(cleanedData, pcVal = 1, cutoff = 4, retTab = F)
fillSkeletonMap(signalMapPC1, modelsPicked = c("OU", "BM", "delta"), title = "PatchSignalMapTestPC1MovedScaled.pdf")
signalTabPC1 <- patchSignalMap(cleanedData, pcVal = 1, cutoff = 4, retTab = T)

signalMapPC2 <- patchSignalMap(cleanedData, pcVal = 2, cutoff = 4)
fillSkeletonMap(signalMapPC2, modelsPicked = c("OU", "BM", "delta"), title = "PatchSignalMapTestPC2MovedScaled.pdf")
signalTabPC2 <- patchSignalMap(cleanedData, pcVal = 2, cutoff = 4, retTab = T)

##AIC comparison 
fitTestAIC <- relativeFitTable(pcVal = 2)
colnames(fitTestAIC) <- unique(cleanedData$ColSpace$Patch)
barplot(fitTestAIC, col = c("coral4", "darkgoldenrod3", "skyblue3", "white"), legend.text = c(rownames(fitTestAIC)), main = "Relative AICc Weights by Patch", pch = .4, cex.names = .8, las = 2, space = .5)

                       
png(filename = here::here("Output/columnMap.png"), res = 500, width = 7, height = 7, bg = "transparent", units = "in")
par(mfrow = c(6, 2))
par(mar = c(0, 0, 0, 0))

makeModelFitMap(modelOptionSet = F, setModelPC = 1)
makeModelFitMap(modelOptionSet = F, setModelPC = 2)

{
  makeModelFitMap(setModel = "lambda")
  makeModelFitMap(setModel = "lambda", setModelPC = 2)
}

{
  makeModelFitMap(setModel = "BM")
  makeModelFitMap(setModel = "BM", setModelPC = 2)
}

{
makeModelFitMap(setModel = "delta")
makeModelFitMap(setModel = "delta", setModelPC = 2)
}

{
makeModelFitMap(setModel = "OU")
makeModelFitMap(setModel = "OU", setModelPC = 2)
}
dev.off()

}


##AIC comparison 
fitTestAIC <- relativeFitTable(pcVal = 2)
colnames(fitTestAIC) <- unique(cleanedData$ColSpace$Patch)
barplot(fitTestAIC, col = c("coral4", "darkgoldenrod3", "skyblue3", "white"), legend.text = c(rownames(fitTestAIC)), main = "Relative AICc Weights by Patch", pch = .4, cex.names = .8, las = 2, space = .5)


#Step 3: Run Climate Analysis and out put PGLsplots

#{
  cleanedDataSet <- cleanedData
  cleanedDataset <- cleanedData
  climDataImport <- read.delim(sep = ",", file = here::here("Data/climData.csv"), header = T, stringsAsFactors = F)
  rownames(climDataImport) <- make.names(climDataImport$tipNames, unique = T)
  climDataImport <- climDataImport[, -c(1)]
  runPGLSanalysis(colPCSelect = 1, brightnessOn = T)
  runPGLSanalysis(colPCSelect = 2, brightnessOn = T)
#}
 


  
  
#ArbutusTable and Color

#SupplementaryAnalyses
  
  par(mfrow = c(5, 5))
  for (i in unique(cleanDataSet$RawDat$Patch))
  {
    aicListTMP <- data.frame(stringsAsFactors = F)
    if(stri_length(i) > 2)
    {
      ##pcaList <- pcaData$x[which(tipLabeledData$Patch == i), pcVal]
      pcData <- cleanedData$ColSpace[which(tipLabeledData$Patch == i), c(20:23)]
      pcaComp<- prcomp(pcData, center = T, scale. = T)  
      plot(pcaComp$x[, c(1, 2)])
      print(summary(pcaComp))
      print(head(pcaComp$rotation))
    }
  }
  

  #Continuous Trait Mapping
  png(filename = here::here("contmappingClearUncertainty2.png"), bg = "transparent", width = 8, height = 4, units = "in", res = 500)
  par(mfrow = c(1, 4))
  contMapForColorVsBodySize(patchSelect = "Crown")
  contMapForColorVsBodySize(patchSelect = "LowerFlank")
  contMapForColorVsBodySize(patchSelect = "Lores")
  contMapForColorVsBodySize(body = T)
  dev.off()
  

  
##supplemental figs
if(supp = T)
{
#fitTables, relative AIC weight

  signalTabPC1 <- patchSignalMap(cleanedData, pcVal = 1, cutoff = 4, retTab = T)
  signalTabPC2 <- patchSignalMap(cleanedData, pcVal = 2, cutoff = 4, retTab = T)
  write.csv(signalTabPC1, file = here::here("SignalTabPC1.csv"))
  write.csv(signalTabPC2, file = here::here("SignalTabPC2.csv"))
  fitTestAIC <- relativeFitTable(pcVal = 2)
  colnames(fitTestAIC) <- unique(cleanedData$ColSpace$Patch)
  barplot(fitTestAIC, col = c("coral4", "darkgoldenrod3", "skyblue3", "white"), legend.text = c(rownames(fitTestAIC)), main = "Relative AICc Weights by Patch", pch = .4, cex.names = .8, las = 2, space = .5)
  
  
#contmap with error bars?
  


#all 35 

for(p in unique(cleanedData$ColSpace$Patch))
{
  par(mfrow = c(1, 1))
  thisNom <- paste("Output/Supplement/contmappingAll", p, ".png", sep = "")
  png(filename = here::here(thisNom), bg = "transparent", width = 2, height = 2, units = "in", res = 300)
  contMapForColorVsBodySize(patchSelect = p)
  dev.off()
}

#PGLS results, PGLS brightness results
  #Step 3: Run Climate Analysis and out put PGLsplots
  
  #{
  cleanedDataSet <- cleanedData
  cleanedDataset <- cleanedData
  climDataImport <- read.delim(sep = ",", file = here::here("Data/climData.csv"), header = T, stringsAsFactors = F)
  rownames(climDataImport) <- make.names(climDataImport$tipNames, unique = T)
  climDataImport <- climDataImport[, -c(1)]
  par(mfrow = c(1, 2))
  PGLSlist <- runPGLSanalysis(colPCSelect = 1, brightnessOn = F, retHist = T)
  hist(PGLSlist)
  patchList <- as.data.frame(cbind(unique(cleanedData$RawDat$Patch), as.numeric(paste(PGLSlist))))
  colnames(patchList) <- c("i", "varPC1")
  patchList[, 2] <- as.numeric(paste(patchList[, 2]))
  joined <- merge(x = patchMapTotal, patchList, by.x = "PatchName", by.y =  "i")
  plot(backgroundMap, border = "black", col = "white")
  coList <- rgb(maxColorValue = 1, red = rescale(to = c(.1, .6), (log(1+as.numeric(paste(joined@data$varPC1))))), blue = rescale(to = c(.1, .6), (log(1+as.numeric(paste(joined@data$varPC1))))), green =  .1)
  plot(add = T, joined, col = coList, border = coList)
  plot(wingMap, add = T, col = alpha("gray", .3))
  
  #}
  
  


}
  
  
  
  
  
  



