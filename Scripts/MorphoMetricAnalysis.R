#Lory Morphometrics
runMorphoAnalysis <- function()
{
loryMorphoSheet <- read.csv("~/Lory MorphoMetrics - MorphoSheetForTipMatch.csv", header = T,stringsAsFactors = F)

hist(loryMorphoSheet$Tail.length..mm.)

plot(loryMorphoSheet[ , c(4:11)])


loryMorphoSheet$TipLabel

loryMorphoSheet[,5]

list <- c(3, 6, 3, 8, 4, 6)

list <- c(1:100)

plot(cleanedData$CleanTree, pch = .5, cex = .5)

cleanDataSetMorpho <- data.frame(stringsAsFactors = F, loryMorphoSheet[ , c(8:10)])
rownames(cleanDataSetMorpho) <- loryMorphoSheet$TipLabel
cleanDataSetMorphoPCA <- prcomp(cleanDataSetMorpho)
summary(cleanDataSetMorphoPCA)
biplot(cleanDataSetMorphoPCA, cex = .5)



matchedMorpho <- match.phylo.data(phy = cleanedData$CleanTree, data = cleanDataSetMorphoPCA$x)

##phenogram(tree = matchedMorpho$phy, x = matchedMorpho$data[ , 1])
contMap(tree = matchedMorpho$phy, x = matchedMorpho$data[, 1])
plot()
ouMorphoFit <- fitContinuous(phy = force.ultrametric(matchedMorpho$phy), dat = matchedMorpho$data, model = "OU")
bmMorphoFit <-fitContinuous(phy = force.ultrametric(matchedMorpho$phy), dat = matchedMorpho$data, model = "BM")
whiteMorphoFit <-fitContinuous(phy = force.ultrametric(matchedMorpho$phy), dat = matchedMorpho$data, model = "white")
kappaMorphoFit <-fitContinuous(phy = force.ultrametric(matchedMorpho$phy), dat = matchedMorpho$data, model = "delta")
print(c(ouMorphoFit$opt$aic, bmMorphoFit$opt$aicc, whiteMorphoFit$opt$aicc, kappaMorphoFit$opt$aicc))
modScores <- c(ouMorphoFit, bmMorphofit, whiteMorphoFit, kappaMorphoFit)



patternDat <- read.delim("C://Users/jmerwin/Downloads/Lory Pattern Data - Side View Full Bird Pattern Data.csv", header = T, sep = ",")
patternDat <- patternDat[-c(1) , ]
patNames <- make.names(patternDat$TipName, unique = T)
rownames(patternDat) <- patNames
patternPC <- prcomp(patternDat[, -c(1, 2)], center = T, scale. = T)
biplot(patternPC, cex = .5)
dm <- dist(x = patternPC$x[ , c(1, 2, 3)], method = "euclidean")
upgmaTree <- upgma(D = dm)
plot(upgmaTree, cex = .5)

plot3d(patternPC$x)
hist(log(patternPC$x[ , 2]+1000))

plot3d(log(patternPC$x+1000), rgb(maxColorValue = 7, patternPC)

patternPCTransformed <- log(patternPC$x[ , 1]+1000)
matchedMorpho <- match.phylo.data(phy = cleanedData$CleanTree, data = patternPCTransformed)
matchedPattern <- match.phylo.data(phy = upgmaTree, data = matchedMorpho$data)
contMap(tree = matchedMorpho$phy, x = matchedMorpho$data)

plot(matchedMorpho$data[ , 2])

compareTrees <- cbind((matchedMorpho$phy$tip.label), (matchedMorpho$phy$tip.label))

plot(cophylo(tr1 = matchedMorpho$phy, tr2 = matchedPattern$phy,rotate = T))





joinList <- data.frame(stringsAsFactors = F)
joinList <- as.data.frame(paste(dirLoad, cleanDataSet$CleanTree$tip.label, ".png", sep = ""), stringsAsFactors = F)
joinList <- cbind.data.frame(cleanDataSet$MetaDat$TreeLabel, joinList, as.data.frame(log10(as.numeric(cleanDataSet$MetaDat$Body.Size)), stringsAsFactors = F), stringsAsFactors = F)
rownames(joinList) <- cleanDataSet$MetaDat$TreeLabel


joinList <- data.frame(stringsAsFactors = F)
joinList <- as.data.frame(paste(dirLoad, cleanDataSet$CleanTree$tip.label, ".png", sep = ""), stringsAsFactors = F)
joinList <- cbind.data.frame(cleanDataSet$MetaDat$TreeLabel, joinList, as.data.frame(log10(as.numeric(cleanDataSet$MetaDat$Body.Size)), stringsAsFactors = F), stringsAsFactors = F)
rownames(joinList) <- cleanDataSet$MetaDat$TreeLabel

}

