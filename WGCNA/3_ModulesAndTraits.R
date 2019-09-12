# R script for WGCNA Modules and Traits
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
lnames <- load(file = "NetworkConstruction_rld_noLSR1.RData")
lnames

nGenesAp = ncol(datExprAp)
nSamplesAp = nrow(datExprAp)
# recalculate MEs with color labels
MEsAp = moduleEigengenes(datExprAp, moduleColorsAp)$eigengenes
MEsAp = orderMEs(MEsAp)
moduleTraitCorAp = cor(MEsAp, datTraitsAp, use = "p")
moduleTraitPvalueAp = corPvalueStudent(moduleTraitCorAp, nSamplesAp)

# and for Buchnera
nGenesBa = ncol(datExprBa)
nSamplesBa = nrow(datExprBa)
MEsBa = moduleEigengenes(datExprBa, moduleColorsBa)$eigengenes
MEsBa = orderMEs(MEsBa)
moduleTraitCorBa = cor(MEsBa, datTraitsBa, use = "p")
moduleTraitPvalueBa = corPvalueStudent(moduleTraitCorBa, nSamplesBa)

# display correlations and their p-values in a heat map
samp <- "Ap"
pdf(paste0("ModuleTraitCor_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
textMatrix = paste(signif(moduleTraitCorAp, 2), "\n(",
                   signif(moduleTraitPvalueAp, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorAp)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCorAp,
               xLabels = names(datTraitsAp),
               yLabels = names(MEsAp),
               ySymbols = names(MEsAp),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, power level 4"))
dev.off()

# and for Buchnera
samp <- "Ba"
pdf(paste0("ModuleTraitCor_rld", samp, "_noLSR1.pdf"), width=12, height=9)
sizeGrWindow(12,9)
textMatrix = paste(signif(moduleTraitCorBa, 2), "\n(",
                   signif(moduleTraitPvalueBa, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorBa)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCorBa,
               xLabels = names(datTraitsBa),
               yLabels = names(MEsBa),
               ySymbols = names(MEsBa),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, power level 9"))
dev.off()

# calculate gene significance and module membership
modNamesAp = substring(names(MEsAp), 3)
geneModuleMembershipAp = as.data.frame(cor(datExprAp, MEsAp, use = "p"))
MMPvalueAp = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipAp), 
                                            nSamplesAp))
names(geneModuleMembershipAp) = paste("MM", modNamesAp, sep="")
names(MMPvalueAp) = paste("p.MM", modNamesAp, sep="")
geneTraitSignificanceAp = as.data.frame(cor(datExprAp, datTraitsAp, use = "p"))
GSPvalueAp = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceAp), 
                                            nSamplesAp))
names(geneTraitSignificanceAp) = paste("GS.", names(datTraitsAp), sep="")
names(GSPvalueAp) = paste("p.GS.", names(datTraitsAp), sep="")

# do it for Buchnera
modNamesBa = substring(names(MEsBa), 3)
geneModuleMembershipBa = as.data.frame(cor(datExprBa, MEsBa, use = "p"))
MMPvalueBa = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipBa), 
                                            nSamplesBa))
names(geneModuleMembershipBa) = paste("MM", modNamesBa, sep="")
names(MMPvalueBa) = paste("p.MM", modNamesBa, sep="")
geneTraitSignificanceBa = as.data.frame(cor(datExprBa, datTraitsBa, use = "p"))
GSPvalueBa = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceBa), 
                                            nSamplesBa))
names(geneTraitSignificanceBa) = paste("GS.", names(datTraitsBa), sep="")
names(GSPvalueBa) = paste("p.GS.", names(datTraitsBa), sep="")

# identify genes with high significance (GS) for titer and high module membership (MM)
# for a module of interest
# "genes highly significantly associated with a trait are often also the most
# important (central) elements of modules associated with the trait"

# using module "blue" as an example for A. pisum
module = "blue"
column = match(module, modNamesAp)
moduleGenesAp = (moduleColorsAp==module)

# plot MM vs. GS for genes from this module
samp <- "Ap"
pdf(paste0("MMvGS_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipAp[moduleGenesAp, column]),
                   abs(geneTraitSignificanceAp[moduleGenesAp, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for titer",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
pdf(paste0("MM_hist_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
hist(abs(geneModuleMembershipAp[moduleGenesAp, column]),
     breaks = 20, xlim = c(0, 1), ylim = c(0, 250),
     xlab = paste("Module Membership"),
     main = paste("Module Membership distribution in", module, "module"))
dev.off()
pdf(paste0("GS_hist_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
hist(abs(geneTraitSignificanceAp[moduleGenesAp, 1]),
     breaks = 20, xlim = c(0, 1), ylim = c(0, 200),
     xlab = paste("Gene Significance for titer"),
     main = paste("Gene Significance for titer distribution in", module, "module"))
dev.off()

# now for Buchnera, using most significantly correlated module (negatively correlated with Tucson parentage)
module = "tan"
column = match(module, modNamesBa)
moduleGenesBa = (moduleColorsBa==module)
samp <- "Ba"
pdf(paste0("MMvGS_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembershipBa[moduleGenesBa, column]),
                   abs(geneTraitSignificanceBa[moduleGenesBa, 3]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for titer",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
pdf(paste0("MM_hist_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
hist(abs(geneModuleMembershipBa[moduleGenesBa, column]),
     breaks = 20, xlim = c(0, 1), ylim = c(0, 8),
     xlab = paste("Module Membership"),
     main = paste("Module Membership distribution in", module, "module"))
dev.off()
pdf(paste0("GS_hist_", module, "_rld", samp, "_noLSR1.pdf"), width=7, height=7)
hist(abs(geneTraitSignificanceBa[moduleGenesBa, 3]),
     breaks = 20, xlim = c(0, 1), ylim = c(0, 8),
     xlab = paste("Gene Significance for parentage"),
     main = paste("Gene Significance for parentage distribution in", module, "module"))
dev.off()

# prepare summary tables of results
# load deseq2 summary tables
resOrderedAp <- read.csv(file = "deseq2_Ap_gene_LRT_resOrderedAp_noLSR1.csv", 
                         header = T, stringsAsFactors = F)
resOrderedBa <- read.csv(file = "deseq2_Ba_gene_LRT_resOrderedBa_noLSR1.csv", 
                         header = T, stringsAsFactors = F)

# obtain the gene order from the datExpr objects 
geneOrderAp <- colnames(datExprAp)
geneOrderBa <- colnames(datExprBa)

# set up tables combining WGCNA results with deseq2 results
resOrderedAp <- resOrderedAp[match(geneOrderAp, resOrderedAp$entrez), ]
geneInfoAp0 <- data.frame(resOrderedAp, moduleColor = moduleColorsAp, 
                          geneTraitSignificanceAp, GSPvalueAp)
resOrderedBa <- resOrderedBa[match(geneOrderBa, resOrderedBa$entrez), ]
geneInfoBa0 <- data.frame(resOrderedBa, moduleColor = moduleColorsBa, 
                          geneTraitSignificanceBa, GSPvalueBa)

# order modules by their significance for titer
sampleTiterAp <- datTraitsAp$sampleTiter
modOrderAp = order(-abs(cor(MEsAp, sampleTiterAp, use = "p")))
# and order Buchnera modules by parentage
sampleTiterBa <- datTraitsBa$sampleParentage
modOrderBa = order(-abs(cor(MEsBa, sampleTiterBa, use = "p")))

# add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembershipAp)) {
  oldNames = names(geneInfoAp0)
  geneInfoAp0 = data.frame(geneInfoAp0, geneModuleMembershipAp[, modOrderAp[mod]],
                           MMPvalueAp[, modOrderAp[mod]])
  names(geneInfoAp0) = c(oldNames, paste("MM.", modNamesAp[modOrderAp[mod]], sep=""),
                         paste("p.MM.", modNamesAp[modOrderAp[mod]], sep=""))
}
# and for Buchnera
for (mod in 1:ncol(geneModuleMembershipBa)) {
  oldNames = names(geneInfoBa0)
  geneInfoBa0 = data.frame(geneInfoBa0, geneModuleMembershipBa[, modOrderBa[mod]],
                           MMPvalueBa[, modOrderBa[mod]])
  names(geneInfoBa0) = c(oldNames, paste("MM.", modNamesBa[modOrderBa[mod]], sep=""),
                         paste("p.MM.", modNamesBa[modOrderBa[mod]], sep=""))
}


# order the genes in the geneInfo variable by deseq2 adj p-value, as in DESeq2 results table
geneOrderAp = order(geneInfoAp0$padj, decreasing = F)
geneOrderBa = order(geneInfoBa0$padj, decreasing = F)

# rename tables, keeping them in order of deseq2 results 
WGCNAresAp = geneInfoAp0[geneOrderAp, ]
WGCNAresBa = geneInfoBa0[geneOrderBa, ]

# save as .csv files
write.csv(WGCNAresAp, file = paste0("WGCNAresAp_noLSR1.csv"), 
          row.names = FALSE)
write.csv(WGCNAresBa, file = paste0("WGCNAresBa_noLSR1.csv"), 
          row.names = FALSE)

# save Rdata for later
save(WGCNAresAp, WGCNAresBa, 
     file = paste0("ModulesAndTraits_rld_noLSR1.RData"))

### module-trait relationship figures
library("flashClust")
library("ape")
library("ggtree")
library("RColorBrewer")
library("pheatmap")

## for both Ap and Ba, make a dendrogram and heatmap of eigengenes relationship to traits
species <- c("Ap", "Ba")
for (i in 1:length(species)) {
  samp <- species[i]
  if (samp == "Ap") {
    MEs <- MEsAp
    moduleTraitCor <- moduleTraitCorAp
    lineSize = 2
  }
  if (samp == "Ba") {
    MEs <- MEsBa
    moduleTraitCor <- moduleTraitCorBa
    lineSize = 3
  }
  # cluster modules by their eigengenes
  eigenDist <- dist(t(MEs))
  eigenClust <- hclust(eigenDist)
  # convert hclust object to phylogenetic tree
  eigenPhylo <- as.phylo(eigenClust)
  
  # plot and save the dendrogram
  png(paste0("ModuleTraitCor_dendrogram_", samp, ".png"), width = 800, 
      height = 1000, res=70)
  plot(ggtree(eigenPhylo, size = lineSize) + 
         theme_tree2() +
         geom_treescale(fontsize = 12, linesize = 4, offset = -2) )
  dev.off()
  
  # convert the phylo object into a ggtree object, from which the order of 
  # branches can be obtained
  tree <- ggtree(eigenPhylo)
  # record the order of branches from the tree
  df <- tree[["data"]]
  df <- df[!is.na(df$label), ]
  df <- df[order(df$y, decreasing = T), ]
  treeOrder <- df$label
  
  # re-order the moduleTraitCor matrix by the tree order
  eigenMatrix <- moduleTraitCor[treeOrder, ]
  
  ## now we make a heatmap
  
  # assign colors for heatmap
  heatmapColors <- colorRampPalette(brewer.pal(8, "RdBu"))(20) %>% rev()
  
  # create vectors with row and column annotations
  rowAnnot <- data.frame(Module = gsub("ME", "", rownames(eigenMatrix)))
  rownames(rowAnnot) <- rownames(eigenMatrix)
  colAnnot <- data.frame(Trait = gsub("sample", "", colnames(eigenMatrix)))
  rownames(colAnnot) <- colnames(eigenMatrix)
  
  # assign colors for row annotations (modules) and column annotations (traits)
  # and organize into a list
  rowColors <- gsub("ME", "", rownames(eigenMatrix))
  names(rowColors) <- gsub("ME", "", rownames(eigenMatrix))
  annotCol <- list(Module = rowColors)
  
  # specify specific points for the heatmap scale
  breaksList <- seq(-1, 1, by = 0.1)
  
  # plot and save the heatmap to a png file
  png(paste0("ModuleTraitCor_pheatmap_", samp, ".png"), width = 800, 
      height = 1000, res=70)
  pheatmap(eigenMatrix, color = heatmapColors, scale = "none",
           cluster_rows = F, cluster_cols = F,
           fontsize = 12, fontsize_row = 12, fontsize_col = 12,
           breaks = breaksList,
           cellwidth = 40,
           cellheight = 30,
           annotation_row = rowAnnot, 
           annotation_col = colAnnot, 
           annotation_names_row = F, annotation_names_col = F,
           show_rownames = F,
           annotation_colors = annotCol, annotation_legend = T,
           angle_col = "45",
           border_color = NA,
           display_numbers = T,          # may want to print one with and one without numbers
           number_color = "black",
           fontsize_number = "12")                 
  dev.off()
}

# number of genes per module
sapply(unique(WGCNAresAp$moduleColor), function(x) 
  { sum(grepl(x, WGCNAresAp$moduleColor)) } )
sapply(unique(WGCNAresBa$moduleColor), function(x) 
  { sum(grepl(x, WGCNAresBa$moduleColor)) } )

# are there any module colors unique to Ba results/ not used for Ap modules?
unique(WGCNAresBa$moduleColor)[!(unique(WGCNAresBa$moduleColor) %in% unique(WGCNAresAp$moduleColor))]




