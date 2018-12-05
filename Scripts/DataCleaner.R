##This function combines raw input data into a linked data frame

dataCleaner <- function(rawColDataLoc, DataKeyLoc, treeLoc)
{
  require(pavo)
  require(picante)
  require(plyr)
  
  rawColData <- read.csv(rawColDataLoc, header = T, stringsAsFactors = F)
  DataKey <- read.csv(DataKeyLoc, header = T, stringsAsFactors = F)
  unrootedTree <- read.tree(treeLoc)
  
  #transform into color space
  newBigSpace <- colspace(rawColData[ , c(10, 12, 14, 16)], space = "tcs", qcatch = "Qi")
  colList <- rgb(red = newBigSpace$l, green = newBigSpace$m, blue = newBigSpace$s)
  
  #attach color space to rawData sheet
  spaceColData <- cbind.data.frame(rawColData, newBigSpace[, c(1:4)])
  
  #tip label the data
  mappedVals <- mapvalues(x = spaceColData$TaxaName, from = DataKey$NewTreeLabel, to = DataKey$TreeLabel)
  tipLabeledData <- cbind.data.frame(stringsAsFactors = F, mappedVals, spaceColData)
  
  #root and trim tree
  unrootedTreeTrim <- ape::drop.tip(phy = unrootedTree, tip = "Psittaculirostris_edwardsii_DOT7827")
  rootedTree <- root(unrootedTreeTrim, outgroup = "Melopsittacus_undulatus_DOT9145", resolve.root = TRUE)
  taxaKey <- DataKey
  rownames(taxaKey) <- taxaKey$TreeLabel
  prunedTree <- match.phylo.data(phy = rootedTree, data = taxaKey)
  trimmedRootedTree <- prunedTree$phy
  par(mfrow = c(1, 2))
  ##pdf(file = "~/LoryProject/FinalManuscriptFolderLoryProject/Figures/FinalLayout/ColSpaceBig48.pdf", bg = "white")
  plot(trimmedRootedTree, cex = .1)
  plot(newBigSpace, col = colList, theta = 0, phi = 14, out.lwd = 2, vert.cex = 3) 
  dev.off()
  resultsList <- (list(tipLabeledData, spaceColData, trimmedRootedTree, prunedTree$data))
  names(resultsList) <- c("RawDat", "ColSpace", "CleanTree", "MetaDat")
  return(resultsList)
}

