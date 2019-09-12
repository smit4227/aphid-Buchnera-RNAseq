# R script for correlation analysis of host and symbiont DESeq2 log-transformed counts 
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"
setwd(wkdir)
lnames <- load(file = "DESeq2ResultsAndGeneInfo_noLSR1.RData")
lnames

## correlation of aphid and Buchnera genes

## identify and remove outlier samples from dataset (from WGCNA tutorial)
library("WGCNA")
library("flashClust")
# load rld data
samp <- "Ap"
rldAp00 <- read.table("RLD_deseq2_Ap_gene_LRT_rldAp_noLSR1.csv", sep = ",", header=T, 
                      stringsAsFactors = F)
rownames(rldAp00) <- rldAp00[,1]
rldAp00[,1] <- NULL
datExprAp <- as.data.frame(t(rldAp00))
# cluster samples with similar patterns of gene expression
# look at resulting pdf, choose a height to remove any outliers, and enter this value as cutHeightAp
# then replot with final cutHeightAp value
sampleTreeAp <- hclust(dist(datExprAp), method = "average")
pdf(paste0("sampleTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeAp, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
cutHeightAp <- 40
abline(h=cutHeightAp, col="red")
dev.off()
# which sample was/were outlier(s)?
outlierAp <- "AL4.2"

## do the same thing for Buchnera
# load rld data
samp <- "Ba"
rldBa00 <- read.table("RLD_deseq2_Ba_gene_LRT_rldBa_noLSR1.csv", sep = ",", header=T, 
                      stringsAsFactors = F)
rownames(rldBa00) <- rldBa00[,1]
rldBa00[,1] <- NULL
datExprBa <- as.data.frame(t(rldBa00))
# cluster samples with similar patterns of gene expression
# look at resulting pdf, choose a height to remove any outliers, and enter this value as cutHeightBa
# then replot with final cutHeightBa value
sampleTreeBa <- hclust(dist(datExprBa), method = "average")
pdf(paste0("sampleTree_rld", samp, "_noLSR1.pdf"), width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTreeBa, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
cutHeightBa <- 5
abline(h=cutHeightBa, col="red")
dev.off()
# which sample was/were outlier(s)?
outlierBa <- "TUCc.3"

## remove from both data sets
rldAp0 <- rldAp00[, -(grep(outlierAp, names(rldAp00)))]
rldAp <- rldAp0[, -(grep(outlierBa, names(rldAp0)))]
rldBa0 <- rldBa00[, -(grep(outlierAp, names(rldBa00)))]
rldBa <- rldBa0[, -(grep(outlierBa, names(rldBa0)))]




## calculate the inner quartile range (IQR) of each gene to use 
# as a measure of variation among samples
samp <- "Ap"
# calculate the IQR for each gene based on rld data
rldApIqr <- apply(as.data.frame(rldAp), 1, IQR)
range(rldApIqr)
mean(rldApIqr)
pdf(paste0("IQR_rld", samp, "_noLSR1.pdf"), width=12, height=9)
hist(rldApIqr[rldApIqr < 2], xlim = c(0,2), breaks = 40, labels=T)
dev.off()
iqrCutoffAp <- 0.25
IQRselectionAp <- rldApIqr > iqrCutoffAp
# what percentage of the total genes does this cutoff represent?
print(paste0("This represents ", ((sum(IQRselectionAp)/length(IQRselectionAp))*100), 
             " % of the total genes")) 

## do the same for Buchnera
samp <- "Ba"
# calculate the IQR for each gene based on rld data
rldBaIqr <- apply(as.data.frame(rldBa), 1, IQR)
range(rldBaIqr)
mean(rldBaIqr)
pdf(paste0("IQR_rld", samp, "_noLSR1.pdf"), width=12, height=9)
hist(rldBaIqr[rldBaIqr < 0.4], xlim = c(0,0.4), breaks = 40, labels=T)
dev.off()
iqrCutoffBa <- 0.1
IQRselectionBa <- rldBaIqr > iqrCutoffBa
# what percentage of the total genes does this cutoff represent?
print(paste0("This represents ", ((sum(IQRselectionBa)/length(IQRselectionBa))*100),
             " % of the total genes")) 



## include in this analysis only those genes with padj < 0.05 from deseq2 
## or exclude those with IQR below the IQR treshhold to remove genes that 
## do not vary between samples or are not significant from the gene Universe

# first load deseq2 results table to obtain padj values
resAp <- read.table("deseq2_Ap_gene_LRT_resOrderedAp_noLSR1.csv", 
                    sep = ",", header=T, stringsAsFactors = F)
rownames(resAp) <- resAp[,1]
resBa <- read.table("deseq2_Ba_gene_LRT_resOrderedBa_noLSR1.csv", 
                    sep = ",", header=T, stringsAsFactors = F)
rownames(resBa) <- resBa[,1]


## now remove genes that don't meet the above criteria
species <- c("Ap", "Ba")
rldList <- list(rldAp, rldBa)
resList <- list(resAp, resBa)
iqrList <- list(rldApIqr, rldBaIqr)
iqrCutoffList <- list(iqrCutoffAp, iqrCutoffBa)
for (i in 1:length(rldList)) {
  rld <- rldList[[i]]
  res <- resList[[i]]
  rldIQR <- iqrList[[i]]
  cutoff <- iqrCutoffList[[i]]
  # ensure row order of res matches that of rld
  res <- res[match(row.names(rld), res$entrez), ]
  # create empty data.frames to fill in
  rldFilt <- data.frame()
  resFilt <- data.frame()
  if (nrow(res) == nrow(rld) & nrow(res) == length(rldIQR)) {
    for (j in 1:nrow(rld)) {
      if (!is.na(res$padj[j])) {
        if (rldIQR[j] > cutoff | res$padj[j] < 0.05) {
          rldFilt <- rbind(rldFilt, rld[j, ])
          resFilt <- rbind(resFilt, res[j, ])
        }
      }
      assign(paste0("rld", species[i], "Filt"), rldFilt)
      assign(paste0("res", species[i], "Filt"), resFilt)
    }
  }
}



## create a matrix of correlations between counts of each pair of aphid and 
## Buchnera genes. This will take a while...
c1 <- rldApFilt
c2 <- rldBaFilt
corMatrixNames <- list(row.names(c1), row.names(c2))
corMatrix <- matrix(nrow = nrow(c1), ncol = nrow(c2), dimnames = corMatrixNames)
for (i in 1:nrow(c2)) {
  BaGene <- row.names(c2)[i]
  for (j in 1:nrow(c1)) {
    ApGene <- row.names(c1)[j]
    countsComp <- data.frame(as.numeric(c1[ApGene,]), as.numeric(c2[BaGene,]))
    countsCor <- cor(countsComp[,1], countsComp[,2], method = "pearson")
    corMatrix[ApGene, BaGene] <- countsCor
  }
}
nSamples <- ncol(rldApFilt)
corMatrixPval = corPvalueStudent(corMatrix, nSamples)

# save the matrix so it can be accessed later
write.csv(as.data.frame(corMatrix), 
          file = paste0("GeneCorrelationMatrix_corMatrix_rld.csv"))
write.csv(as.data.frame(corMatrixPval), 
          file = paste0("GeneCorrelationMatrixPval_corMatrixPval_rld.csv"))

# look at the distribution of correlation coefficients in the matrix
pdf(paste0("ApBaCorrelations_corMatrix_noLSR1.pdf"), width=12, height=9)
hist(corMatrix, xlim = c(-1,1), ylim = c(0, 160000), breaks = 20, labels=T)
dev.off()
# look at the distribution of p-values in the matrix
pdf(paste0("ApBaCorrelationPvalues_corMatrixPval_noLSR1.pdf"), width=12, height=9)
hist(corMatrixPval, xlim = c(0,1), ylim = c(0, 200000), breaks = 20, labels=T)
dev.off()

#### make a heatmap of the correlation matrix, with genes clustered by similar 
#### patterns of gene correlation

## cluster rows and columns by similar correlation profiles
## first for Buchnera
library("stats")
library("flashClust")
library("WGCNA")
library("dendextend")
colDist <- dist(t(corMatrix))
colClust <- hclust(colDist)
pdf("Ba_cor_tree.pdf")
plot(colClust, main = "Ba correlation dendrogram", xlab = "", sub = "")
dev.off()
# combine similar branches and group genes together
colClustCut <- cutree(colClust, k = 12)
table(colClustCut)
colClustCutCol <- labels2colors(colClustCut)
table(colClustCutCol)
names(colClustCutCol) <- colClust[["labels"]]
pdf("Ba_cor_tree_modules.pdf")
print(plotDendroAndColors(colClust, colClustCutCol, "Dynamic Tree Cut",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = "Ba correlation dendrogram with clusters"))
dev.off()


## now clustering for A. pisum
rowDist <- dist(corMatrix)
rowClust <- hclust(rowDist)
pdf("Ap_cor_tree.pdf")
plot(rowClust, main = "Ap correlation dendrogram", xlab = "", sub = "")
dev.off()
# combine similar branches and group genes together
rowClustCut <- cutree(rowClust, k=16)
table(rowClustCut)
rowClustCutCol <- labels2colors(rowClustCut)
table(rowClustCutCol)
names(rowClustCutCol) <- rowClust[["labels"]]
pdf("Ap_cor_tree_modules.pdf")
print(plotDendroAndColors(rowClust, rowClustCutCol, "Dynamic Tree Cut",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = "Ap correlation dendrogram with clusters"))
dev.off()

## generate a heatmap for correlations between individual genes
# function to convert color name to hex code
col2hex <- function(color) {
  df <- col2rgb(color)
  hex <- rgb(df[1,], df[2,], df[3,], max = 255)
  return(hex)
}
# set colors for heatmap
library("pheatmap")
# assign colors to each gene entrez ID for both Ap and Ba
rowAnnotation <- data.frame(ApModule = rowClustCutCol)
rownames(rowAnnotation) <- rowClust[["labels"]]
colAnnotation <- data.frame(BaModule = colClustCutCol)
rownames(colAnnotation) <- colClust[["labels"]]
# define color to be assigned to each module name and arrange it in a list
rowCol <- unique(rowClustCutCol) %>% col2hex() # not neccesary to convert to hex code but I wrote it and it's mine
names(rowCol) <- unique(rowClustCutCol)
colCol <- unique(colClustCutCol) %>% col2hex() # not neccesary to convert to hex code but I wrote it and it's mine
names(colCol) <- unique(colClustCutCol)
annotationColors <- list(ApModule = rowCol, BaModule = colCol)
# set up a diverging color scheme for correlation data so we can see both neg and pos correlations
library("RColorBrewer")
heatmapColors <- colorRampPalette(brewer.pal(8, "Spectral"))(20) %>% rev()   # order is blue to red (neg to pos)
# define legend breakpoints (for the 21 colors in heatmapColors)
breakList <- seq(-1, 1, by = 0.1)
# make the gene correlation heatmap using pretty heatmap package

# first save the heatmap as an object, to extract the gene order later
corMap <- pheatmap(corMatrix, color = heatmapColors, scale = "none",
                   main = "Gene correlation heatmap",
                   cutree_rows = length(rowCol), cutree_cols = length(colCol),
                   legend_breaks = breakList, 
                   fontsize = 6, fontsize_row = 1, fontsize_col = 1,
                   annotation_row = rowAnnotation, annotation_col = colAnnotation, 
                   annotation_names_row = F, annotation_names_col = F,
                   annotation_colors = annotationColors, annotation_legend = T,
                   breaks = breakList,
                   angle_col = "45")
dev.off()

# then save as image files. First, a pdf you can zoom in on to enable you to look 
# at row/column labels
pdf(paste0("gene_rlog_correlation_heatmap.pdf"), width = 8.5, height = 45)
pheatmap(corMatrix, color = heatmapColors, scale = "none",
         main = "Gene correlation heatmap",
         cutree_rows = length(rowCol), cutree_cols = length(colCol),
         legend_breaks = breakList, 
         fontsize = 6, fontsize_row = 1, fontsize_col = 1,
         annotation_row = rowAnnotation, annotation_col = colAnnotation, 
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = annotationColors, annotation_legend = T,
         breaks = breakList,
         angle_col = "45")
dev.off()

# then as a png for smaller image size and to use as an actual figure
png(paste0("gene_rlog_correlation_heatmap.png"), width = 800, height = 1000, res=100)
pheatmap(corMatrix, color = heatmapColors, scale = "none",
         main = "Gene correlation heatmap",
         cutree_rows = length(rowCol), cutree_cols = length(colCol),
         legend_breaks = breakList, 
         fontsize = 10, fontsize_row = 0.1, fontsize_col = 0.1,
         annotation_row = rowAnnotation, annotation_col = colAnnotation, 
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = annotationColors, annotation_legend = T,
         breaks = breakList,
         angle_col = "45")
dev.off()



## generate matrix figures for p-value matrix, too

# obtain gene order of corMap object
rowOrder <- corMap$tree_row[["order"]]
colOrder <- corMap$tree_col[["order"]]
# arrange p-value matrix accordingly
corMatrixPvalOrdered <- corMatrixPval[rowOrder, colOrder]
# new unidirectional color scheme and break numbers
heatmapColors <- colorRampPalette(brewer.pal(8, "YlOrRd"))(20) %>% rev()    # significant p-values are red
breakList <- seq(0, 1, by = 0.05)
# create vectors to specify where to introduce gaps between gene modules
# first, arrange gene colors in the same order as corMap
rowCCCordered <- rowClustCutCol[rowOrder]
colCCCordered <- colClustCutCol[colOrder]
# establish the order of module colors in corMap
rowColOrdered <- unique(rowCCCordered)
colColOrdered <- unique(colCCCordered)
# how many genes for each module color?
rowColCount <- table(rowCCCordered)[rowColOrdered]
colColCount <- table(colCCCordered)[colColOrdered]
# the vector to specify gap positions should add the number of genes in the next module to the 
# number of genes in the previous
rowGap <- c(0)
for (i in 1:length(rowColCount)) {
  previous <- rowGap[i]
  current <- rowColCount[i] + previous
  rowGap <- c(rowGap, current)  
}
colGap <- c(0)
for (i in 1:length(colColCount)) {
  previous <- colGap[i]
  current <- colColCount[i] + previous
  colGap <- c(colGap, current)  
}

# heatmap object
corPvalMap <- pheatmap(corMatrixPval, color = heatmapColors, scale = "none",
                       main = "Gene correlation p-value heatmap",
                       cluster_cols = F, cluster_rows = F,
                       gaps_row = rowGap, gaps_col = colGap,
                       legend_breaks = breakList, 
                       fontsize = 6, fontsize_row = 1, fontsize_col = 1,
                       annotation_row = rowAnnotation, annotation_col = colAnnotation, 
                       annotation_names_row = F, annotation_names_col = F,
                       annotation_colors = annotationColors, annotation_legend = T,
                       breaks = breakList,
                       angle_col = "45")
dev.off()

png(paste0("gene_rlog_correlation_pval_heatmap.png"), width = 800, height = 1000, 
    res=100)
pheatmap(corMatrixPvalOrdered, color = heatmapColors, scale = "none",
         main = "Gene correlation p-value heatmap",
         cluster_cols = F, cluster_rows = F,
         gaps_row = rowGap, gaps_col = colGap,
         legend_breaks = breakList, 
         fontsize = 10, fontsize_row = 0.1, fontsize_col = 0.1,
         annotation_row = rowAnnotation, annotation_col = colAnnotation, 
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = annotationColors, annotation_legend = T,
         breaks = breakList,
         angle_col = "45")
dev.off()


# save key objects to use later
save(rldApFilt, rldBaFilt, resApFilt, resBaFilt, 
     corMatrix, corMatrixPval, corMap, corPvalMap,
     rowClustCutCol, colClustCutCol,
     file = paste0("CorrelationAnalysis_rld_noLSR1.RData"))
