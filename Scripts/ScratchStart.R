##scratch
##CodetogetStarted
##FInal Analysis

       

##clean data
library(plyr)
library(rgdal)
library(phytools)
library(geiger)
library(scales)
library(stringi)
library(here)
library(e1071)


matchedTree <- cleanedData$CleanTree
patchMapTotal <- readOGR(dsn = path.expand(here("Shapes")), layer = "PatchNamedPatchMapFeb7")
backgroundMap <- readOGR(dsn = path.expand(here("Shapes")), layer = "patches1")
wingMap <- readOGR(dsn = path.expand("C://Users/jmerwin/Documents/LoryProject"), layer = "wingLine")
shp <- patchMapTotal






#Input and Clean Data
cleanedData <- dataCleaner(rawColDataLoc = here("Data/FullCompleteMasterSheet.csv"), DataKeyLoc = here("Data/FullMatchTreeD70Full.csv"), treeLoc = here("Data/lcr20_all_treePL.tre"))
signalPatchMap <- patchSignalMap(cleanedData, retTab = F, pcVal = 2)


    


{
     bigDataCols <- cbind.data.frame(stringsAsFactors = F, signalPatchMap[-c(1), 1], as.numeric(signalPatchMap[-c(1), 2]), signalPatchMap[-c(1), 3])
     ##modelChoice <- c("OU", "BM", "white")   
    modelChoice <- c("OU", "BM", "white")
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
           pdf(file = here("Output/PatchSignalMapTestPC2.pdf"), bg = "white")
           
           plot(backgroundMap, col = "black", border = "black")
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
           
        
           
## careful this takes 30 mins.
crownCols <- birdPatchSpaceStats(cleanDataSet, patch = "Crown")
matchup <- match.phylo.data(matchup$phy, crownCols[ , c(21:24)])
surfaceRunCrown <- runSurface(tree =nameNodes(matchup$phy), as.data.frame(matchup$data, stringsAsFactors = F), verbose = TRUE, error_skip = T)    

par(mfrow = c(1, 1))
surfaceSumm <- surfaceSummary(fwd = surfaceRunLores$fwd, bwd = surfaceRunLores$bwd)
colsList <- rgb(blue = surfaceSumm$theta[ , 2], green = surfaceSumm$theta[ , 3], red = surfaceSumm$theta[ , 4], alpha = .97)
surfaceTreePlot(tree = nameNodes(matchup$phy), hansenfit = surfaceRunLores$bwd[[9]], cex = .7, convcol = F,
                  col = colsList, edge.width = 7) 

surfaceTraitPlot(dat = matchup$data, surfaceRunLores$bwd[[9]], convcol = T, optellipses = F)
surfaceAICPlot(out = surfaceRunFirstPass, traitplot = "aic")
surfaceAICPlot(out = surfaceRunFirstPass, traitplot = "dev")
surface::surfaceSummary(fwd = surfaceRunFirstPass$fwd, bwd = surfaceRunFirstPass$bwd)
           

namedDat <- (matchup$data$centroid.u)
names(namedDat) <- rownames(matchup$data)
contMap(res = 1000, matchup$phy, namedDat,ftype = "reg", fsize = .25, lwd = 2, outline = F)
dttTry <- dtt(phy = matchup$phy, namedDat, plot = T, nsim = 1000, calculateMDIp = T)

    climDataImport <- read.csv(here::here("Data/climData.csv"), head = T)
    rownames(climDataImport) <- climDataImport$tipNames
    climDataImport <- climData[ , -c(1)]
    
    setwd(here::here("Output"))
    for(j in c(1:35))
    {
    i <- unique(cleanedData$RawDat$Patch)[j]
    print(i)
    
    name <- path.expand(here::here("Output", "PGLS", paste(i, ".png", sep = "")))
    png(filename = name, res = 200, width = 10, height = 5, units = "in", bg = "white")
    par(mfrow = c(1, 2))
    thisPGLS <- byPatchPGLS(patch = i, cleanedDataset = cleanedData, climData = climDataImport, colpcVal = 1)
    
    summ <- caper::summary.pgls(thisPGLS)
    plot(thisPGLS$data$data$V1~thisPGLS$data$data$PC1, main = paste(i, summ$coefficients[2, 4]), sub = paste( summ$adj.r.squared))
    text(y = thisPGLS$data$data$V1+.5, x = thisPGLS$data$data$PC1+.4, cex = .2)
    flag <- F
    for(i in summ$coefficients[ , 4]){if(i<.001428){flag <- T}}
    red1 <- 0
    if(flag){red1 = 1}
    abline(thisPGLS$model$coef[c(1,2)], col = rgb(red = red1,blue = 0, green = 0)) 
    dev.off()
    }

#dependencies and source other files




#Step1: run DataCleaner, Output Trimmed Tree and Cleaned Colorspace. 




#Step 2: import blank patchmap, output 4 patchmaps for 2 PCs each, 8 pdf outputs, variance, brownian rate, sigma value, arbutus maps





#Step 3: Run Climate Analysis and out put PGLs maps
















           