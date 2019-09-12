# R script for WGCNA Functional Enrichment Analysis for Buchnera genes
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
samp <- "Ba"

## defining the gene Universe (the gene background within which a subset of genes 
# can show enrichment of specific functions)

# construction of geneList with geneNumber replaced by ENTREZID and 
# including log2fold change
allGeneInfoBa <- read.csv(file = "WGCNAresBa_noLSR1.csv", 
                          stringsAsFactors = F, header = T)
allGeneInfoBa <- allGeneInfoBa[order(allGeneInfoBa$geneNumber), ]

# use iqr to identify genes whose expression profile varies little between samples, and remove them from the gene Universe
lnames <- load("DataInput_rld_noLSR1.RData")
lnames
datExprIqr <- apply(t(datExprBa), 1, IQR)
pdf(paste0("IQR_rld", samp, "_noLSR1.pdf"), width=12, height=9)
hist(datExprIqr[datExprIqr < 0.4], xlim = c(0,0.4), breaks = 40, labels=T)
dev.off()
iqrCutoff <- 0.10
selected <- datExprIqr > iqrCutoff
# what percentage of the total genes does this cutoff represent?
print(paste0("This represents ", ((sum(selected)/length(selected))*100), 
             "% of the total genes")) 
# filter genes by selected
datExprFilt <- datExprIqr[selected]
# match the geneNumber of datExprFilt to locus in allGeneInfoAp
IQRGeneInfoBa <- allGeneInfoBa[(allGeneInfoBa$geneNumber %in% names(datExprFilt)), ]

# create a vector containing log2FoldChange with entrez IDs for row names
geneList <- allGeneInfoBa$log2FoldChange
names(geneList) <- allGeneInfoBa$entrez
# KEGG enrichment uses "locus_tag" feature only (ENTREZ not supported)
geneList2 <- allGeneInfoBa$log2FoldChange
names(geneList2) <- allGeneInfoBa$locus_tag

# do the same for the IQR-filtered genes
IQRgeneList <- IQRGeneInfoBa$log2FoldChange
names(IQRgeneList) <- IQRGeneInfoBa$entrez
IQRgeneList2 <- IQRGeneInfoBa$log2FoldChange
names(IQRgeneList2) <- IQRGeneInfoBa$locus_tag


# choose interesting gene groups to do functional annotation enrichment analyses on
# I chose gene modules most strongly correlated with parentage, and also titer
intModules = c("brown", "yellow", "royalblue", "cyan", "tan", "turquoise", "darkgrey") 
intModulesList <- list()
IQRintModulesList <- list()
# pick out entrez numbers and log2foldchange for all genes in each module of interest
moduleColors <- allGeneInfoBa$moduleColor[allGeneInfoBa$entrez %in% names(geneList)]
IQRmoduleColors <- allGeneInfoBa$moduleColor[allGeneInfoBa$entrez %in% names(IQRgeneList)]
for (module in intModules) {
  # select module probes
  modGenes = (moduleColors==module)
  IQRmodGenes = (IQRmoduleColors==module)
  # get their geneID codes
  modLocustag = names(geneList2)[modGenes]
  modLocustag = modLocustag[!is.na(modLocustag)]
  IQRmodLocustag = names(IQRgeneList2)[IQRmodGenes]
  IQRmodLocustag = IQRmodLocustag[!is.na(IQRmodLocustag)]
  # save objects in a list
  intModulesList[[module]] <- modLocustag
  IQRintModulesList[[module]] <- IQRmodLocustag
}


## GO enrichment with clusterProfiler is not possible because no OrgDb for Buchnera

## KEGG enrichment with clusterProfiler is still possible because KEGG includes 
## Buchnera annotation use the IQR filtered gene list to remove genes with 
## low variation between samples

intModules = as.vector(names(intModulesList))
pval = 0.05
qval = 0.05
kkModulesList <- list()
for (i in 1:length(IQRintModulesList)) {
  # go enrichment for each module of interest
  geneNames <- IQRintModulesList[[i]]
  universe <- names(geneList2)
  kk <- enrichKEGG(gene = geneNames, universe = universe, 
                   organism = 'buc', pAdjustMethod = "fdr", pvalueCutoff = pval, 
                   qvalueCutoff = qval)
  objectName <- paste0(intModules[i], "KEGG")
  print(i)
  kkModulesList[[objectName]] <- kk
}
sapply(kkModulesList, nrow)

# only one module had enough genes for KEGG enrichment, so it's not worth making cnetplots
kkModulesList[["turquoiseKEGG"]]$ID


## packages for analyses requiring geneList object
library("DOSE")
library("ggplot2")
library("ggraph")
library("igraph")
library("enrichplot")
library("pathview")
library("KEGGREST")
library("png")



## KEGG pathview

# pull out KEGG pathway identifiers and do pathview on each one for each module of interest
# open and run pathviewTES_complete.R in workingDir to load necessary functions
source("pathviewTES_complete.R")

# identify pathways of interest
pathways <- c("buc01230", "buc00010", "buc00020", "buc00030", "buc00051", "buc00520", 
              "buc00190", "buc00230", "buc00240", "buc00250", "buc00260", "buc00270", 
              "buc00290", "buc00300", "buc00220", "buc00330", "buc00340", "buc00400", 
              "buc00670", "buc00790", "buc00760", "buc00550", "buc00970", "buc03010", 
              "buc03020", "buc03060", "buc04122", "buc03018", "buc03030", "buc02010", 
              "buc02060", "buc02040", "buc00410", "buc00480")

# start KEGG pathway
pathviewResultsList <- list()
folderName <- paste0("pathviewTES_Ba_gene_LRT_noLSR1")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName))
}
# since there are so few genes in each module, it makes little sense to do pathview 
# on each separate module
# unlist IQRintModulesList and do pathview on all genes from all modules of interest
IQRintModuleGenes <- unlist(IQRintModulesList) 
for (j in 1: length(pathways)) {
  k <- pathviewTES(gene.data = names(geneList), pathway.id = pathways[j], 
                   species = "buc", gene.idtype = "Entrez",
                   limit = list(gene= c(-1,1), cpd = 1))
  resultsName <- paste0(pathways[j])
  pathviewResultsList[[resultsName]] <- k
}
setwd(wkdir)

## still, it might be easier to look at Buchnera genes by hand, since there are so few
# which Ba genes are in each Buchnera module?
goiBa <- allGeneInfoBa[allGeneInfoBa$entrez %in% names(IQRgeneList), ]
# save table
write.csv(as.data.frame(goiBa), file = paste0("WGCNAgoiBa_rld_filtered.csv"), 
          row.names = F)

# rename and save enrichment results
WGCNAkeggListBa <- kkModulesList
WGCNApathviewListBa <- pathviewResultsList
WGCNAgoiBa <- goiBa

save(WGCNAkeggListBa, WGCNApathviewListBa, WGCNAgoiBa,
     file = paste0("FunctionalEnrichment_Ba_rld_noLSR1.RData"))
