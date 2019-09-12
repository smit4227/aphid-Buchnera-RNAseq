# R script for WGCNA Network Construction
# written by Thomas E. Smith
# last updated Sep 2019
# for use in RStudio, or remove comment from line 13
# based on instructions in https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

# setting the stage
wkdir <- "/Users/User/Work"
setwd(wkdir)
library("WGCNA")
library("flashClust")
options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
lnames <- load(file = "DataInput_rld_noLSR1.RData")
lnames

# selection of soft-thresholiding power by analysis of expression network fit 
# with scale-free topology model
samp <- "Ap"
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExprAp, powerVector = powers, verbose = 5)
pdf(paste0("softThresholdPlots_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2))
cex1 = 0.8
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# the next line adds a line at an R^2 cutoff value of 0.8 
abline(h=0.90,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# same for Buchnera
samp <- "Ba"
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExprBa, powerVector = powers, verbose = 5)
pdf(paste0("softThresholdPlots_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
par(mfrow = c(1,2))
cex1 = 0.8
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# the next line adds a line at an R^2 cutoff value of 0.8 
abline(h=0.80,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# look at 'softThresholdPlots.pdf' and decide on the smallest power level that 
# gets above R^2 = 0.8
# set soft thresholding power and calculate adjacencies
softPowerAp <- 4
adjacencyAp <- adjacency(datExprAp, power = softPowerAp)
softPowerBa <- 9
adjacencyBa <- adjacency(datExprBa, power = softPowerBa)

# transform adjacencies into a Topological Overlap Matrix, then calculate dissimilarity
TOMAp = TOMsimilarity(adjacencyAp)
dissTOMAp = 1-TOMAp
TOMBa = TOMsimilarity(adjacencyBa)
dissTOMBa = 1-TOMBa

# clustering using TOM
geneTreeAp = hclust(as.dist(dissTOMAp), method = "average")
geneTreeBa = hclust(as.dist(dissTOMBa), method = "average")

# If "Error: vector memory exhausted (limit reached?)" is obtained, open Terminal and enter the following commands
# cd ~
# touch .Renviron
# open .Renviron
# as the first line in .Renviron, enter "R_MAX_VSIZE=100Gb" and save 

# plot the resulting gene clustering dendrograms
samp <- "Ap"
pdf(paste0("geneTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
print(plot(geneTreeAp, xlab="", sub="", 
           main = "Gene clustering on TOM-based dissimilarity, power level = 4",
           labels = FALSE, hang = 0.04))
dev.off()
# and for Buchnera
samp <- "Ba"
pdf(paste0("geneTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
print(plot(geneTreeBa, xlab="", sub="", 
           main = "Gene clustering on TOM-based dissimilarity, power level = 9",
           labels = FALSE, hang = 0.04))
dev.off()

# set the minimum module size relatively high
minModuleSizeAp = 30
minModuleSizeBa = 10

# module identification using dynamic tree cut:
dynamicModsAp = cutreeDynamic(dendro = geneTreeAp, distM = dissTOMAp,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSizeAp)
table(dynamicModsAp)
# for Buchnera
dynamicModsBa = cutreeDynamic(dendro = geneTreeBa, distM = dissTOMBa,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSizeBa)
table(dynamicModsBa)

# convert numeric lables into colors
dynamicColorsAp = labels2colors(dynamicModsAp)
table(dynamicColorsAp)
dynamicColorsBa = labels2colors(dynamicModsBa)
table(dynamicColorsBa)

# plot the dendrogram and colors underneath
samp <- "Ap"
pdf(paste0("geneTreeModules_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
print(plotDendroAndColors(geneTreeAp, dynamicColorsAp, "Dynamic Tree Cut",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = "Gene dendrogram and module colors, power level = 4"))
dev.off()
# plot the dendrogram and colors underneath
samp <- "Ba"
pdf(paste0("geneTreeModules_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
print(plotDendroAndColors(geneTreeBa, dynamicColorsBa, "Dynamic Tree Cut",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = "Gene dendrogram and module colors, power level = 9"))
dev.off()

# merge modules with similar expression profiles
# calculate eigengenes
samp <- "Ap"
MEList = moduleEigengenes(datExprAp, colors = dynamicColorsAp)
MEsAp = MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEsAp)
# cluster module eigengenes
METreeAp = hclust(as.dist(MEDiss), method = "average")
MEDissThresAp = 0.3
# plot results
pdf(paste0("METree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
print(plot(METreeAp, main = "Clustering of module eigengenes",
           xlab = "", ylab="dissimilarity", sub = ""))
abline(h=MEDissThresAp, col = "red")
dev.off()
# look at METree.pdf and choose a cutoff height, below which modules will be merged
MEDissThresAp = 0.37
# go back and plot again, changes in MEDissThres will move the abline in the plot

# for Buchnera, too
samp <- "Ba"
MEList = moduleEigengenes(datExprBa, colors = dynamicColorsBa)
MEsBa = MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEsBa)
# cluster module eigengenes
METreeBa = hclust(as.dist(MEDiss), method = "average")
MEDissThresBa = 0.3
# plot results
pdf(paste0("METree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
print(plot(METreeBa, main = "Clustering of module eigengenes",
           xlab = "", ylab="dissimilarity", sub = ""))
abline(h=MEDissThresBa, col = "red")
dev.off()
# look at METree.pdf and choose a cutoff height, below which modules will be merged
MEDissThresBa = 0.2
# go back and plot again, changes in MEDissThres will move the abline in the plot


# merge modules and assign new colors to merged modules
mergeAp = mergeCloseModules(datExprAp, dynamicColorsAp, cutHeight = MEDissThresAp, 
                            verbose = 3)
mergedColorsAp = mergeAp$colors
# eigengenes of the new merged modules:
mergedMEsAp = mergeAp$newMEs
# and for Buchnera
mergeBa = mergeCloseModules(datExprBa, dynamicColorsBa, cutHeight = MEDissThresBa, 
                            verbose = 3)
mergedColorsBa = mergeBa$colors
mergedMEsBa = mergeBa$newMEs

# re-plot results with merged modules      ##### MEDissThresAp from before is 0.37
samp <- "Ap"
pdf(paste0("geneTreeMergedModules_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTreeAp, cbind(dynamicColorsAp, mergedColorsAp),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# and for Buchnera
samp <- "Ba"
pdf(paste0("geneTreeMergedModules_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTreeBa, cbind(dynamicColorsBa, mergedColorsBa),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# rename objects for next step and save R environment
moduleColorsAp = mergedColorsAp
moduleColorsBa = mergedColorsBa
# construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabelsAp = match(moduleColorsAp, colorOrder)-1
moduleLabelsBa = match(moduleColorsBa, colorOrder)-1
MEsAp = mergedMEsAp
MEsBa = mergedMEsBa
save(MEsAp, MEsBa, moduleLabelsAp, moduleLabelsBa, moduleColorsAp, 
     moduleColorsBa, geneTreeAp, geneTreeBa, 
     file = paste0("NetworkConstruction_rld_noLSR1.RData"))
