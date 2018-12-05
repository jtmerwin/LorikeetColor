patchVarSkeletonFill <- function(patchVarOutput)
{
  
  pcDatCols3<- matrix(byPatchPCVariance(cleanedData, raw = F, ouOption = F, modelPC = 1))
  pcDatCols2<- ((as.numeric(paste(pcDatCols3[[2]]))))
  
  pcDatCols <- data.frame(levels(varList$i), pcDatCols2)
  joined <- merge(x = patchMapTotal, pcDatCols, by.x = "PatchName", by.y =  "levels.varList.i.")
  
  ##plots patchmaps and writes to file
  
  namePlot = paste(toString(cleanDataSet$RawDat[((chunk)*numPatches)-1, 1]), ".png", sep = "")
  png(filename= paste(here::here("Output"), sep = "/", namePlot), units = "in", height = 5, width = 5, res = 200, bg = "transparent")
  plot(backgroundMap, border = "black", col = "white")
  coList <- rgb(maxColorValue = 5, red = 5-(as.numeric(paste(joined@data$pcDatCols))), green = 2, blue =  2)
  plot(add = T, joined, col = coList, border = coList)
  plot(wingMap, add = T, col = alpha("gray", .3))
  
  legend(cex = .5, title = "Color Variance by Patch",   # location of legend
  legend = c("Maximum Variance", "Minimum Variance"), # categories or elements to render in the legend
  fill = range(joined@data$varPC1))
  
}