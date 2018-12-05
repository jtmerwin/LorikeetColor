#' @title makePatchMaps of Color Data.
#'
#' @description
#' \code{makePatchMaps} writes patchmaps to file.
#'
#' @details
#' ColorData must be a 7 column data frame object with Species names (or filenames) in the
#' first column, "PatchName" in the second column. columns 3, 4, 5, and 6
#' must be UV, Short, Medium, and Longwave reflectance data in that order.
#' All data must normalize to 1. Ex: .20, .23, .27, .30 for u, s, m, l. 
#' Patchmap must be a shapefile with internal polygons with values for a "PatchName" in a VAT.
#' Number of polygons must be consistent across different specimens.
#' Numpatches is the number of patches for which you have data.
#' You must have equal numbers of patches for both colorData and patchMap
#' fileDestination is the address for the folder where you wish to save the generated PNGs.
makePatchMaps <- function(cleanDataSet, patchMap, numPatches, fileDestination)
{
  require("rgdal")
  require("png")
  
  patchMapTotal <- patchMap
  chunkSize <- numPatches
  
  #sets up loop to loop through all entries in list, chunking by numPatches
  done <- FALSE
  chunk <- 1
  while(!done)
  {
    firstColorSetNamed <- cbind.data.frame(stringsAsFactors = F, cleanDataSet$RawDat$Patch[c(1:35)], cleanDataSet$RawDat[c((35*chunk-34):(35*chunk)), c(21:24)])
    joined <- merge(x = patchMapTotal, firstColorSetNamed, by.x = "PatchName", by.y = "cleanDataSet$RawDat$Patch[c(1:35)]")
    
    #plots patchmaps and writes to file
    setwd(fileDestination)
    namePlot = paste(toString(cleanDataSet$RawDat[((chunk)*numPatches)-1, 1]), ".png", sep = "")
    png(filename= namePlot, units = "in", height = 5, width = 5, res = 200, bg = "transparent")
    plot(backgroundMap, border = "black", col = "black")
    plot(add = T, joined, col = rgb(red = joined@data$l, blue = joined@data$s, green = joined@data$m), border = rgb(red = joined@data$l, blue = joined@data$s, green = joined@data$m))
    plot(wingMap, add = T, col = alpha("gray", .3))    
    dev.off()
    
    #checks if loop is done, notify user if done
    chunk <- (chunk + 1)
    if((chunk*numPatches) > NROW(cleanDataSet$RawDat))
    {
      done = TRUE
      print("done")
    }
  }
}





