# R script for WGCNA Data Input and Cleaning
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
lnames <- load(file = "DESeq2ResultsAndGeneInfo_noLSR1.RData")
lnames

# load log-transformed counts from DESeq2
# must be read from .csv file, rld is an incompatible object type
ApData <- read.csv("RLD_deseq2_Ap_gene_LRT_rldAp_noLSR1.csv")
datExprAp0 <- as.data.frame(t(ApData[, -1]))
names(datExprAp0) <- ApData$X
rownames(datExprAp0) <- names(ApData[, -1])
# for Buchnera, too
BaData <- read.csv("RLD_deseq2_Ba_gene_LRT_rldBa_noLSR1.csv")
datExprBa0 <- as.data.frame(t(BaData[, -1]))
names(datExprBa0) <- BaData$X
rownames(datExprBa0) <- names(BaData[, -1])

# check for columns with excessive missing values ("NA"s) and identify outliers
gsg <- goodSamplesGenes(datExprAp0, verbose=3)
gsg$allOK
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprAp0)[!gsg$goodGenes], 
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprAp0)[!gsg$goodSamples],
                                                collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprAp0 = datExprAp0[gsg$goodSamples, gsg$goodGenes]
}
# now for Buchnera
gsg <- goodSamplesGenes(datExprBa0, verbose=3)
gsg$allOK
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprBa0)[!gsg$goodGenes], 
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprBa0)[!gsg$goodSamples],
                                                collapse = ", ")));
  datExprBa0 = datExprBa0[gsg$goodSamples, gsg$goodGenes]
}

# cluster samples by expression
samp <- "Ap"
sampleTreeAp <- hclust(dist(datExprAp0), method = "average")
pdf(paste0("sampleTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeAp, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# look at pdf, choose a height to remove any outliers, and enter this value as cutHeightAp
# replot with final cutHeightAp value
cutHeightAp <- 40
abline(h=cutHeightAp, col="red")
dev.off()
# cut the tree to remove the outlier
clustAp <- cutreeStatic(sampleTreeAp, cutHeight = cutHeightAp, minSize = 10)
table(clustAp)
# shows the clusters formed, clust 1 contains the desired samples
keepSamples <- (clustAp==1)
datExprAp <- datExprAp0[keepSamples, ]
nGenesAp <- ncol(datExprAp)
nSamplesAp <- nrow(datExprAp)

# repeat all of this for Buchnera
samp <- "Ba"
sampleTreeBa <- hclust(dist(datExprBa0), method = "average")
pdf(paste0("sampleTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeBa, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# look at pdf, choose a height to remove any outliers, and enter this 
# value as cutHeightAp
# replot with final cutHeightAp value
cutHeightBa <- 5
abline(h=cutHeightBa, col="red")
dev.off()
# cut the tree to remove the outlier
clustBa <- cutreeStatic(sampleTreeBa, cutHeight = cutHeightBa, minSize = 10)
table(clustBa)
# shows the clusters formed, clust 1 contains the desired samples
keepSamples <- (clustBa==1)
datExprBa <- datExprBa0[keepSamples, ]
nGenesBa <- ncol(datExprBa)
nSamplesBa <- nrow(datExprBa)

# format trait data for WGCNA
sampleColor <- c("green", "green", "green", "red", "red", "red", 
                 "red", "red", "red", "green", "green", "green", 
                 "green", "green", "green", "green", "green", "green", 
                 "green", "green", "green")
sampleParentage <- c("A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "T", "T",
                     "T", "T", "T", "T", "T", "T", "A", "A", "A")
sampleTiter <- c(42, 42, 42, 82, 82, 82, 62, 62, 62, 25, 25, 25, 15, 15, 15, 
                 40, 40, 40, 35, 35, 35)
# apply to A. pisum data
names(sampleTiter) <- names(ApData[, 2:22])
names(sampleColor) <- names(ApData[, 2:22])
names(sampleParentage) <- names(ApData[, 2:22])
datTraitsAp0 <- as.data.frame(cbind(sampleTiter, sampleColor, sampleParentage))
keepColumns <- (rownames(datTraitsAp0) %in% rownames(datExprAp))
datTraitsAp <- as.data.frame(datTraitsAp0[keepColumns,])
# binarize sampleColor and sampleParentage, "red" is 1 and "green" is 0, 
# "T" is 1 and "A" is 0
datTraitsAp <- binarizeCategoricalColumns(datTraitsAp, convertColumns = 
                                            c("sampleColor", "sampleParentage"))
# invert sampleParentage so that "T" is 0 and "A" is 1
datTraitsAp$sampleParentage.T.vs.all <- abs(datTraitsAp$sampleParentage.T.vs.all - 1)
colnames(datTraitsAp) <- colnames(datTraitsAp0)

# now for Buchnera
names(sampleTiter) <- names(BaData[, 2:22])
names(sampleColor) <- names(BaData[, 2:22])
names(sampleParentage) <- names(BaData[, 2:22])
datTraitsBa0 <- as.data.frame(cbind(sampleTiter, sampleColor, sampleParentage))
keepColumns <- (rownames(datTraitsBa0) %in% rownames(datExprBa))
datTraitsBa <- as.data.frame(datTraitsBa0[keepColumns,])
# binarize sampleColor and sampleParentage, "red" is 1 and "green" is 0, 
# "T" is 1 and "A" is 0
datTraitsBa <- binarizeCategoricalColumns(datTraitsBa, convertColumns = 
                                            c("sampleColor", "sampleParentage"))
# invert sampleParentage so that "T" is 0 and "A" is 1
datTraitsBa$sampleParentage.T.vs.all <- abs(datTraitsBa$sampleParentage.T.vs.all - 1)
colnames(datTraitsBa) <- colnames(datTraitsBa0)

# cluster samples again with outlier sample removed
samp <- "Ap"
sampleTreeAp2 <- hclust(dist(datExprAp), method = "average")
pdf(paste0("sampleTree2_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeAp2, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
# and for Buchnera
samp <- "Ba"
sampleTreeBa2 <- hclust(dist(datExprBa), method = "average")
pdf(paste0("sampleTree2_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeBa2, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# below lines can be commented out if only one trait of interest is present
samp <- "Ap"
titerColors <- numbers2colors(as.numeric(datTraitsAp$sampleTiter), signed=FALSE)
colorColors <- numbers2colors(as.numeric(datTraitsAp$sampleColor), signed=FALSE)
parentageColors <- numbers2colors(as.numeric(datTraitsAp$sampleParentage), signed=FALSE)
traitColors <- as.matrix(cbind(titerColors, colorColors, parentageColors))
pdf(paste0("sampleTree2andTiter_rld", samp, "_noLSR1.pdf"), width=12, height=9)
plotDendroAndColors(sampleTreeAp2, traitColors, 
                    groupLabels = names(datTraitsAp), 
                    main = "Sample dendogram and trait heatmap")
dev.off()
# and for Buchnera
samp <- "Ba"
titerColors <- numbers2colors(as.numeric(datTraitsBa$sampleTiter), signed=FALSE)
colorColors <- numbers2colors(as.numeric(datTraitsBa$sampleColor), signed=FALSE)
parentageColors <- numbers2colors(as.numeric(datTraitsBa$sampleParentage), signed=FALSE)
traitColors <- as.matrix(cbind(titerColors, colorColors, parentageColors))
pdf(paste0("sampleTree2andTiter_rld", samp, "_noLSR1.pdf"), width=12, height=9)
plotDendroAndColors(sampleTreeBa2, traitColors, 
                    groupLabels = names(datTraitsBa), 
                    main = "Sample dendogram and trait heatmap")
dev.off()

# save objects and R environment 
save(datExprAp, datExprBa, datTraitsAp, datTraitsBa,
     file = paste0("DataInput_rld_noLSR1.RData"))

