substr(a$data,1,nchar(a$data)-3)


testPatternList <- makeList(prepath = "~/PatternizeTest/", IDlist = list.files(path = "~/PatternizeTest/"), type = "image")

ID <- list.files("~/PatternizeTest/SomuLightenTest/")
ID <- c("Lorius_lory_somu802073_0000", "Lorius_lory_somu802074_0000", "Lorius_lory_somu802076_0000", "Lorius_lory_somu802077_0000", "Lorius_lory_somu802079_0000")  
testPatternList1 <- makeList(prepath = "~/PatternizeTest/SomuLightenTest/", IDlist = ID, type = 'image', extension = ".tiff")

testRGBRasterList <- patRegK(maskOutline = thisOutLine, sampleList = testPatternList1, target = testPatternList1[[1]], plot = T, resampleFactor = 2, k = 5)
rastSum <- sumRaster(rList = testRGBRasterList, IDlist = ID, type = "k")
plotHeat(summedRaster = rastSum, IDlist = ID)


colnames(thisOutLine) <- c("x", "y")

testRGBRasterList <- patRegK(sampleList = testPatternList, target = testPatternList1[[1]], plot = T, maskOutline = thisOutLine, crop = c(0, 8, 0, 6), k = 3)
rastSum <- sumRaster(rList = testRGBRasterList, IDlist = ID, type = "k")
plotHeat(testRGBRasterList, IDlist = ID)


thisOutLine <- read.delim(sep = "\t", header = F, file = "~/PatternizeTest/JustOne/Chalcopsitta_cardinalis217110_0000.txt")

colnames(thisOutLine) <- c("x", "y")
thisOutLine <- as.matrix(thisOutLine)
thisOutline <- maskOutline(RasterStack = testRGBRasterList[[1]], outline = thisOutLine)


?createTarget
testRaster <- sumRaster(rasterList_regK, IDlist = IDlist)
plotHeat(rasterList_regK$BC0077)
plot(landmarkArray)


testPatternList <- makeList(prepath = "~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/", IDlist = list.files(path = "~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/"), type = "image")

ID <- list.files(path = "~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/")

testRGBRasterList <- patRegK(useBlockPercentage = 75, sampleList = testPatternList, plot = T, k = 4, target = testPatternList$Chalcopsitta_atra_atra_AMNH233675.png)
rastSum <- sumRaster(rList = testPatternList[c(1:4)], IDlist = ID, type = "k")
plotHeat(summedRaster = rastSum, IDlist = ID)



x = raster("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Chalcopsitta_atra_atra_AMNH233675.png")
rast1 <- raster.gaussian.smooth(x = raster("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Chalcopsitta_atra_atra_AMNH233675.png"), sigma = 100, n = 21)
rast2 <- raster.gaussian.smooth(x = raster("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Chalcopsitta_duivenbodei_duivenbodei_DOT13128.png"), sigma = 100, n = 21)
par(mfrow = c(1, 3))
plot(rast1)
plot(rast2)
plot(rast1-rast2)

testJPG <- magick::image_read(path = "~/PatternizeTest/JPG/1x/Asset 1-100.jpg")
testSquid <- magick::image_read(path = "~/dots (2).png")

testBird1 <- image_read("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Chalcopsitta_duivenbodei_duivenbodei_DOT13128.png")
testBird2 <- image_read("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Chalcopsitta_duivenbodei_syringanuchalis_AMNH828937.png")
testBird3 <- image_read("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Trichoglossus_flavoviridis_flavoviridis_AMNH618606.png")
testBird4 <- image_read("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/Vini_peruviana_DOT13079.png")

listBords <- list.files("~/LoryProject/FinalPackagedScripts/Shapes/TipNamedPatchMaps/", full.names = T)



testBords <- image_read(path = )
image_flatten(operator = "Minus", image_read(listBords))

bords <- NULL
bords <- image_animate(image_morph(frames = 5, image_resize(image_read(listBords), geometry = "1000x1000>")))
image_write(bords, "~/allbirdanim.gif")

morphBords<- image_morph(bords, frames = 5)


animation3 <- image_animate(image_morph(image_scale(bords), frames = 8))
image_write(animation3, "~/allbirdanim.gif")

image_morph(bords)


bords <- lapply(listBords, magick::image_read)
??ImageList

animation3 <- image_animate((bords))
image_write(animation3, "~/allbirdanim.gif")


plot((testSquid))
plot(image_edge(testJPG))
plot(image_charcoal(testSquid))

squids <- c(testSquid, image_implode(testSquid))
bords <- c(testBird1, testBird2, testBird3, testBird4, testBird4, testBird3, testBird2, testBird1)
(animation3 <- image_animate(image_morph(unlist(bords), frames = 10)))
image_write(animation3, "~/birdanim.gif")


# Morph effect  <-- result of this is shown
(animation2 <- image_animate(image_morph(squids, frames = 10)))
image_write(animation2, "~/anim2.gif")

colHistList <- NULL



list <- sapply(1:length(cleanedData$RawDat$mappedVals), function(x, y) which.max(y))





