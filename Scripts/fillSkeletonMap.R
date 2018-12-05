fillSkeletonMap <- function(mapFilledList, modelsPicked = c("OU", "BM", "white"), title = NULL)
{
  bigDataCols <- cbind.data.frame(stringsAsFactors = F, mapFilledList[-c(1), 1], as.numeric(mapFilledList[-c(1), 2]), mapFilledList[-c(1), 3])
  modelChoice <- modelsPicked
  colnames(bigDataCols) <- c("patch", "colVals", "colMod")
  colConv <- cbind(modelChoice, c("coral4", "darkgoldenrod3", "skyblue3"))
  
  ##filterList <- grep(pattern = "0", x = bigDataCols[ , 2], value = F, )
  filterList <- which(as.numeric(bigDataCols$colVals) == 0)
  
  cols <- mapvalues(bigDataCols[ , 3], from = colConv[,1], to = colConv[,2])
  cols[filterList] <- "gray"
  
  colnames(bigDataCols) <- c("patch", "colVals")
  bigDataJoin <- cbind.data.frame(stringsAsFactors = F, bigDataCols, cols)
  joinedSig <- sp::merge(x = patchMapTotal, bigDataJoin, by.x = "PatchName", by.y = "patch", duplicateGeoms = TRUE)
  colsSig <- joinedSig@data$colVals
  pdf(file = here("Output", title), bg = "white")
  
  plot(backgroundMap, col = "white", border = "black")
  list <- cols
  ##supCols <- (as.numeric(cols))
  plot(add = T, joinedSig, col = joinedSig@data$cols, border = joinedSig@data$cols)
  plot(wingMap, add = T, col = alpha("white", .5))
  
  legend(cex = .5, title = "OU Vs. BM, FILTERED", "bottom",   # location of legend
         legend = colConv[c(1:3), 1], # categories or elements to render in
         # the legend
         fill = colConv[c(1:3) , 2])
  dev.off()
}