patchMapTreePlot <- function(cleanDataSet, dirLoad = here::here("Shapes/TipNamedPatchMaps/"))
{
  library(ggtree)
  q <- ggtree(cleanDataSet$CleanTree, layout = "rectangular", branch.length = "none")
  joinList <- data.frame(stringsAsFactors = F)
  joinList <- as.data.frame(paste(dirLoad, "/", cleanDataSet$CleanTree$tip.label, ".png", sep = ""), stringsAsFactors = F)
  joinList <- cbind(as.character(cleanDataSet$MetaDat$TreeLabel), joinList, as.data.frame(log10(as.numeric(cleanedData$MetaDat$Body.Size))))
  rownames(joinList) <- cleanDataSet$MetaDat$TreeLabel
  colnames(joinList) <- c("tipNames", "fileExt", "bodySize")
  q<- (q %<+% data.frame(stringsAsFactors = F, joinList))
  namePlot <- here::here("ParrotGGTreeOutput18.png")
  png(width = 10, height = 10, units = "in", res = 1000, filename = namePlot)
  par(mar=c(0,2,0,2))
  ##q + xlim_tree(4) +  xlim(NA, 400) + geom_tiplab(aes(image = fileExt),  size = rescale(as.numeric(as.character(q$data$bodySize[c(1:98)])), to = c(.8, 2))/20 , geom="image", align=T)
  q + xlim_tree(4) +  xlim(NA, 400) + geom_tiplab(aes(image = fileExt),  size = log(3+as.numeric(as.character(q$data$bodySize[c(1:98)])))/20 , geom="image", align=T)
}
