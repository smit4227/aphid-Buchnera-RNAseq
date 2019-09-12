# R script for calculating eigengenes and module membership of correlation analysis modules
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"
setwd(wkdir)
lnames <- load(file = "CorrelationAnalysis_rld_noLSR1.RData")
lnames

## create a correlation matrix of module eigengenes between Ap and Ba modules

# calculate module eignegenes for Ba
library("WGCNA")
MEListBa <- moduleEigengenes(t(rldBaFilt), colors = colClustCutCol)
MEsBa <- MEListBa$eigengenes
rownames(MEsBa) <- colnames(rldBaFilt)
# calculate module eignegenes for Ap
MEListAp <- moduleEigengenes(t(rldApFilt), colors = rowClustCutCol)
MEsAp <- MEListAp$eigengenes
rownames(MEsAp) <- colnames(rldApFilt)
# order the modules eigengenes to match the order of the gene correlation heatmap
rowGeneOrder <- rowClustCutCol[ corMap$tree_row[["order"]] ] # the order of Ap genes with module colors
rowColorOrder <- unique(rowGeneOrder) # the order of the Ap module colors
colGeneOrder <- colClustCutCol[ corMap$tree_col[["order"]] ] # the order of Ba genes with module colors
colColorOrder <- unique(colGeneOrder) # the order of the Ba module colors
MEsAp <- MEsAp[, match(rowColorOrder, gsub("ME", "", colnames(MEsAp)))]
MEsBa <- MEsBa[, match(colColorOrder, gsub("ME", "", colnames(MEsBa)))]
# calculate percentage of variability explained by each eigengene
MEsVarAp <- propVarExplained(t(rldApFilt), rowClustCutCol, MEsAp, corFnc = "cor", corOptions = "use = 'p'")
MEsVarAp <- (MEsVarAp * 100)
MEsVarBa <- propVarExplained(t(rldBaFilt), colClustCutCol, MEsBa, corFnc = "cor", corOptions = "use = 'p'")
MEsVarBa <- (MEsVarBa * 100)
# create eigengene correlation matrix
e1 <- t(MEsAp)
e2 <- t(MEsBa)
eMatrixNames <- list(row.names(e1), row.names(e2))
eMatrix <- matrix(nrow = nrow(e1), ncol = nrow(e2), dimnames = eMatrixNames)
for (i in 1:nrow(e2)) {
  BaMod <- row.names(e2)[i]
  for (j in 1:nrow(e1)) {
    ApMod <- row.names(e1)[j]
    eComp <- data.frame(as.numeric(e1[ApMod,]), as.numeric(e2[BaMod,]))
    eCor <- cor(eComp[,1], eComp[,2], method = "pearson")
    eMatrix[ApMod, BaMod] <- eCor
  }
}
# calculate p-values for gene correlations
nSamples <- nrow(MEsAp)                    ######## should be 19 aphid samples
eMatrixPval = corPvalueStudent(eMatrix, nSamples)
# set up annotation for heatmap
rowAnnotation <- data.frame(ApModule = rev(rowColorOrder))
rownames(rowAnnotation) <- rev(rownames(eMatrix))
colAnnotation <- data.frame(BaModule = colColorOrder)
rownames(colAnnotation) <- colnames(eMatrix)
# define color to be assigned to each module name and arrange it in a list
rowColors <- rowColorOrder  
names(rowColors) <- rowColorOrder 
colColors <- colColorOrder 
names(colColors) <- colColorOrder
annotationColors <- list(ApModule = rowColors, BaModule = colColors)

## create a heatmap of eigengene correlation matrix
library("pheatmap")
library("RColorBrewer")
# establish colors to be used for heatmap
heatmapColors <- rev(colorRampPalette(brewer.pal(8, "Spectral"))(20))   # order is blue to red
# define legend breakpoints (for the 21 colors in heatmapColors)
breakList <- seq(-1, 1, by = 0.1)


png(paste0("eigengene_rlog_correlation_heatmap.png"), width = 800, height = 1000, 
    res=100)
pheatmap(eMatrix, color = heatmapColors, scale = "none",
         main = "Eigengene correlation heatmap",
         cluster_cols = F, cluster_rows = F,
         legend_breaks = breakList, 
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         annotation_row = rowAnnotation, annotation_col = colAnnotation, 
         annotation_names_row = T, annotation_names_col = T,
         annotation_colors = annotationColors, annotation_legend = F,
         breaks = breakList,
         angle_col = "45")
dev.off()




## correlate global correlation analysis module eigengenes with titer

# calculate module membership of each gene
MMAp = as.data.frame(cor(t(rldApFilt), MEsAp, method = "pearson"))
MMPvalueAp <- corPvalueStudent(as.matrix(MMAp), nSamples = 19)
MMBa = as.data.frame(cor(t(rldBaFilt), MEsBa, method = "pearson"))
MMPvalueBa <- corPvalueStudent(as.matrix(MMBa), nSamples = 19)

# format sampleTiter to correlate titer with module eigengenes
# sampleTiter with one AL4 and one TucC sample removed (titer = 62 and 40, respectively)
sampleTiter <- c(42, 42, 42, 82, 82, 82, 62, 62, 25, 25, 25, 15, 15, 15, 40, 40, 
                 35, 35, 35)
names(sampleTiter) <- colnames(rldApFilt)
# also include parentage for Buchnera
# 0 for Tucson, 1 for Austin
sampleParentage <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
names(sampleParentage) <- colnames(rldApFilt)
# combine titer and parentage into a single data.frame, with samples as rownames
datTraits <- data.frame(titer = sampleTiter, parentage = sampleParentage)

# correlate eigengenes with titer, calculate p-values of correlation
ApModuleTiterCor <- cor(MEsAp, sampleTiter, use = "p")
ApModuleTiterCorP <- corPvalueStudent(ApModuleTiterCor, nSamples = 19)
# order MEsVar to match ApModuleTiterCor
MEsVarApOrdered <- MEsVarAp[match( gsub("ME", "", rownames(ApModuleTiterCor)), 
                                   gsub("PVE", "", names(MEsVarAp)) )] 
ApModuleSummary <- data.frame(var.explained = MEsVarApOrdered,
                              cor.titer = ApModuleTiterCor[,1],
                              cor.titer.p.val = ApModuleTiterCorP[,1])
BaModuleTraitsCor <- cor(MEsBa, datTraits, use = "p")
BaModuleTraitsCorP <- corPvalueStudent(BaModuleTraitsCor, nSamples = 19)
# order MEsVar to match ApModuleTiterCor
MEsVapBaOrdered <- MEsVarBa[match( gsub("ME", "", rownames(BaModuleTraitsCor)), 
                                   gsub("PVE", "", names(MEsVarBa)) )]
BaModuleSummary <- data.frame(var.explained = MEsVapBaOrdered,
                              cor.titer = BaModuleTraitsCor[,1],
                              cor.titer.p.val = BaModuleTraitsCorP[,1],
                              cor.parentage = BaModuleTraitsCor[,2],
                              cor.parentage.p.val = BaModuleTraitsCorP[,2])
# add correlation coefficients with other modules and p-values
for (mod in 1:ncol(eMatrix)) {
  oldNames <- names(ApModuleSummary)
  ApModuleSummary <- data.frame(ApModuleSummary,
                                eMatrix[ , mod],
                                eMatrixPval[ , mod])
  names(ApModuleSummary) <- c(oldNames, 
                              paste0("Ba.", gsub("ME", "", colnames(eMatrix))[mod]), 
                              paste0("p.Ba.", gsub("ME", "", colnames(eMatrixPval))[mod]) )
}
rownames(ApModuleSummary) <- gsub("PVE", "Ap.", rownames(ApModuleSummary))
# and for Buchnera
for (mod in 1:nrow(eMatrix)) {
  oldNames <- names(BaModuleSummary)
  BaModuleSummary <- data.frame(BaModuleSummary,
                                eMatrix[mod, ],
                                eMatrixPval[mod, ])
  names(BaModuleSummary) <- c(oldNames, 
                              paste0("Ap.", gsub("ME", "", rownames(eMatrix))[mod]), 
                              paste0("p.Ap.", gsub("ME", "", rownames(eMatrixPval))[mod]) )
}
rownames(BaModuleSummary) <- gsub("PVE", "Ba.", rownames(BaModuleSummary))

# save module summaries as tables
write.csv(as.data.frame(ApModuleSummary), 
          file = paste0("ApModuleSummary_rld.csv"))
write.csv(as.data.frame(BaModuleSummary), 
          file = paste0("BaModuleSummary_rld.csv"))





## plot A. pisum module eigengenes against titer
library("ggplot2")
library("reshape")
# function to convert color name to hex code
col2hex <- function(color) {
  df <- col2rgb(color)
  hex <- rgb(df[1,], df[2,], df[3,], max = 255)
  return(hex)
}
# plot relationship of eigengenes to titer or both A. pisum and Buchnera
species <- c("Ap", "Ba")
for (i in 1:length(species)) {
  if (species[i] == "Ap") {
    MEs <- MEsAp
  }
  if (species[i] == "Ba") {
    MEs <- MEsBa
  }
  dat0 <- cbind(sampleTiter, MEs)
  dat <- melt(dat0, id.vars = "sampleTiter", value.name = "RLD", variable = "Module")
  dat$Module <- gsub("ME", "", dat$Module)
  dat <- dat[order(dat$Module), ]
  MEcol <- col2hex(unique(dat$Module))  
  # calculate mean rld counts for each module + titer level combination
  meandat <- aggregate(dat$value, by = list(sampleTiter = dat$sampleTiter,
                                            Module = dat$Module), FUN = mean)
  
  # plot of module trendlines for titer, from linear model of eigengenes
  line_plot0 <- ggplot(dat, aes(x = sampleTiter, y = value, color = Module)) +
    geom_smooth(method = lm, se = F, size =2, show.legend = F) +
    scale_color_manual(values = MEcol)
  line_plot <- line_plot0 + coord_cartesian(ylim = c(-0.4, 0.4)) +
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "grey90"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(7, "pt"),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color= "black", size = 20), 
          axis.title = element_blank() )
  pdf(paste0(species[i], "_EigengeneTrends_v_Titer.pdf"))
  plot(line_plot)
  dev.off()
  
  # plot of module eigengenes for titer, with each point representing a mean of 2-3 replicates
  line_plot0 <- ggplot(meandat, aes(x = sampleTiter, y = x, color = Module)) +
    geom_line(size = 2, show.legend = F) + 
    scale_color_manual(values = MEcol)
  line_plot <- line_plot0 + coord_cartesian(ylim = c(-0.4, 0.4)) +
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "grey90"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(7, "pt"),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color= "black", size = 20), 
          axis.title = element_blank() )
  pdf(paste0(species[i], "_Eigengenes_v_Titer.pdf"))
  plot(line_plot)
  dev.off()
  
  # plot legend alone
  # use a dataframe like dat but with NAs for x and y values
  emptydat <- data.frame(sampleTiter = rep(NA, nrow(dat)),
                         Module = dat$Module,
                         value = rep(NA, nrow(dat)) )
  legend_plot0 <- ggplot(emptydat, aes(x = sampleTiter, y = value, color = Module)) +
    geom_line(na.rm = T, size = 2) +
    scale_color_manual(values = MEcol)
  legend_plot <- legend_plot0 + theme(plot.title = element_blank(),
                                      panel.background = element_blank(),
                                      panel.grid.major =  element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.line = element_blank(),
                                      axis.text = element_blank(), 
                                      axis.title = element_blank(),
                                      axis.ticks = element_blank(),
                                      legend.title = element_blank(),
                                      legend.key = element_rect(fill = NA),
                                      legend.text = element_text(size = 20))
  pdf(paste0(species[i], "_legend.pdf"))
  plot(legend_plot)
  dev.off()
}



## correlate individual genes with titer, calculate p-values of correlation
ApTiterCor <- cor(t(rldApFilt), sampleTiter, use = "p")
ApTiterCorP <- corPvalueStudent(ApTiterCor, nSamples = 19)
ApTiterSignificance <- data.frame(correlation = ApTiterCor[,1],
                                  p.val = ApTiterCorP[,1] )
BaTraitsCor <- cor(t(rldBaFilt), datTraits, use = "p")
BaTraitsCorP <- corPvalueStudent(BaTraitsCor, nSamples = 19)
BaTraitsSignificance <- data.frame(titer.correlation = BaTraitsCor[,1],
                                  titer.p.val = BaTraitsCorP[,1],
                                  parentage.correlation = BaTraitsCor[,2],
                                  parentage.p.val = BaTraitsCorP[,2])

# load deseq2 results tables and format
resOrderedAp <- read.csv("deseq2_Ap_gene_LRT_resOrderedAp_noLSR1.csv",
                         quote = "\"", header = F, stringsAsFactors = F)
colnames(resOrderedAp) <- resOrderedAp[1, ]
resOrderedAp <- resOrderedAp[-1, ]
resOrderedBa <- read.csv("deseq2_Ba_gene_LRT_resOrderedBa_noLSR1.csv",
                         quote = "\"", header = F, stringsAsFactors = F)
colnames(resOrderedBa) <- resOrderedBa[1, ]
resOrderedBa <- resOrderedBa[-1, ]

# remove genes not included in the correlation analysis
resOrderedApFilt <- resOrderedAp[resOrderedAp$entrez %in% names(rowClustCutCol), ]
resOrderedBaFilt <- resOrderedBa[resOrderedBa$entrez %in% names(colClustCutCol), ]

# make sure gene order matches between resOrderedFilt, ClustCol, and TraitsCor
# are there any rows that don't match?
any(!(resOrderedApFilt$entrez == names(rowClustCutCol)))
any(!(resOrderedApFilt$entrez == row.names(ApTiterCor)))
any(!(resOrderedBaFilt$entrez == names(colClustCutCol)))
any(!(resOrderedBaFilt$entrez == row.names(BaTraitsCor)))
any(!(names(colClustCutCol) == row.names(BaTraitsCor)))
# order resOrderedBaFilt to match colClustCol
resOrderedBaFilt <- resOrderedBaFilt[match(names(colClustCutCol), resOrderedBaFilt$entrez), ]
any(!(resOrderedBaFilt$entrez == row.names(BaTraitsCor)))

# assign names to modules
modNamesAp <- substring(names(MEsAp), 3)
modNamesBa <- substring(names(MEsBa), 3)
# order modules by their significance for titer (or parentage for Buchnera)
modOrderAp <- order(-abs(cor(MEsAp, sampleTiter, use = "p")))
modOrderBa <- order(-abs(cor(MEsBa, datTraits[,2], use = "p")))

## construct an initial summary table for A. pisum
geneInfoAp0 = data.frame(resOrderedApFilt,
                         moduleColor = rowClustCutCol, 
                         GS.sampleTiter = ApTiterCor,
                         p.GS.sampleTiter = ApTiterCorP)
# add module membership information to the table in order of significance for titer
for (i in 1:ncol(MMAp)) {
  oldNames = names(geneInfoAp0)
  geneInfoAp0 = data.frame(geneInfoAp0,
                           MMAp[, modOrderAp[i]],
                           MMPvalueAp[, modOrderAp[i]])
  names(geneInfoAp0) = c(oldNames, paste("MM.", modNamesAp[modOrderAp[i]], sep=""),
                         paste("p.MM.", modNamesAp[modOrderAp[i]], sep=""))
}
# save the results
write.csv(geneInfoAp0, file = paste0("CorrelationAnalysis_Ap.csv"))


## make a similar table for Buchnera
geneInfoBa0 = data.frame(resOrderedBaFilt,
                         moduleColor = colClustCutCol, 
                         GS.sampleTiter = BaTraitsCor[,1],
                         p.GS.sampleTiter = BaTraitsCorP[,1],
                         GS.sampleParentage = BaTraitsCor[,2],
                         p.GS.sampleParentage = BaTraitsCorP[,2])

# add module membership information to the table in order of significance for titer
for (i in 1:ncol(MMBa)) {
  oldNames = names(geneInfoBa0)
  geneInfoBa0 = data.frame(geneInfoBa0,
                           MMBa[, modOrderBa[i]],
                           MMPvalueBa[, modOrderBa[i]])
  names(geneInfoBa0) = c(oldNames, paste("MM.", modNamesBa[modOrderBa[i]], sep=""),
                         paste("p.MM.", modNamesBa[modOrderBa[i]], sep=""))
}
# save the results
write.csv(geneInfoBa0, file = paste0("CorrelationAnalysis_Ba.csv"))

# rename some things
geneInfoAp <- geneInfoAp0
geneInfoBa <- geneInfoBa0

# save key objects to use later
save(MEsAp, MEsBa, MMAp, MMBa, eMatrix, eMatrixPval, 
     ApModuleSummary, BaModuleSummary, sampleTiter,
     geneInfoAp, geneInfoBa,
     file = paste0("EigengenesOfCorrelationAnalysis_rld_noLSR1.RData"))
