# R script for WGCNA Functional Enrichment Analysis for A. pisum genes
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
samp <- "Ap"

## defining the gene Universe (the gene background within which a subset of genes
# can show enrichment of specific functions)

# construction of geneList with geneNumber replaced by ENTREZID and 
# including log2fold change
allGeneInfoAp <- read.csv(file = "WGCNAresAp_noLSR1.csv", 
                          stringsAsFactors = F, header = T)
allGeneInfoAp <- allGeneInfoAp[order(allGeneInfoAp$geneNumber), ]

# use iqr to identify genes whose expression profile varies little between samples, and remove them from the gene Universe
lnames <- load("DataInput_rld_noLSR1.RData")
lnames
datExprIqr <- apply(t(datExprAp), 1, IQR)
pdf(paste0("IQR_rld", samp, "_noLSR1.pdf"), width=12, height=9)
hist(datExprIqr[datExprIqr < 2], xlim = c(0,2), breaks = 30, labels=T)
dev.off()
iqrCutoff <- 0.25
selected <- datExprIqr > iqrCutoff
# what percentage of the total genes does this cutoff represent?
print(paste0("This represents ", ((sum(selected)/length(selected))*100), 
             "% of the total genes")) 
# filter genes by selected
datExprFilt <- datExprIqr[selected]
# match the geneNumber of datExprFilt to locus in allGeneInfoAp
IQRGeneInfoAp <- allGeneInfoAp[(allGeneInfoAp$entrez %in% names(datExprFilt)), ]

# create a vector containing log2FoldChange with entrez IDs for row names
geneList <- allGeneInfoAp$log2FoldChange
names(geneList) <- allGeneInfoAp$entrez

# do the same for the IQR-filtered genes
IQRgeneList <- IQRGeneInfoAp$log2FoldChange
names(IQRgeneList) <- IQRGeneInfoAp$entrez


# choose interesting gene groups to do functional annotation enrichment analyses on
# I chose gene modules most strongly correlated with titer
intModules = c("blue", "green") 
intModulesList <- list()
IQRintModulesList <- list()
# pick out entrez numbers and fold2change for all genes in each module of interest
moduleColors <- allGeneInfoAp$moduleColor[allGeneInfoAp$entrez %in% names(geneList)]
IQRmoduleColors <- allGeneInfoAp$moduleColor[allGeneInfoAp$entrez %in% names(IQRgeneList)]
for (module in intModules) {
  # select module probes
  modGenes = (moduleColors==module)
  IQRmodGenes = (IQRmoduleColors==module)
  # get their geneID codes
  modEntrez = names(geneList)[modGenes]
  modEntrez = modEntrez[!is.na(modEntrez)]
  IQRmodEntrez = names(IQRgeneList)[IQRmodGenes]
  IQRmodEntrez = IQRmodEntrez[!is.na(IQRmodEntrez)]
  # save objects in a list
  intModulesList[[module]] <- modEntrez
  IQRintModulesList[[module]] <- IQRmodEntrez
}


## GO enrichment with clusterProfiler
# use the IQR filtered gene list to remove genes with low variation between samples

ontologies <- c("MF", "BP", "CC")
intModules = as.vector(names(intModulesList))
pval = 0.05
qval = 0.05
# enrich by ENTREZID as opposed to any other keytype
# Set argument readable = F because conversion to gene symbolis not helpful for A. pisum,
# as there is little info availble under SYMBOL keytype.
egoModulesList <- list()
for (i in 1:length(IQRintModulesList)) {
  # perform enrichment for each ontology
  for (j in 1:length(ontologies)) {
    # go enrichment for each module of interest
    geneNames <- IQRintModulesList[[i]]
    universe <- names(geneList)
    ego <- enrichGO(gene = geneNames, universe = universe, 
                    OrgDb = apisum, ont = ontologies[j], pAdjustMethod = "fdr", 
                    pvalueCutoff = pval, qvalueCutoff = qval, keyType = "ENTREZID")
    objectName <- paste0(intModules[i], ontologies[j])
    egoModulesList[[objectName]] <- ego
  }
}
sapply(egoModulesList, nrow)

# simplify GO enrichment results by combining similar GO terms with simplify function
egoModulesSimp <- list()
Cutoff <- 0.6
for (i in 1:length(egoModulesList)) {
  ego <- egoModulesList[[i]]
  egoName <- paste0(names(egoModulesList)[i])
  egoSimp <- clusterProfiler::simplify(ego, cutoff = Cutoff, by = "p.adjust",
                                       select_fun = min)
  egoModulesSimp[[egoName]] <- egoSimp
}
sapply(egoModulesSimp, nrow)


## KEGG enrichment with clusterProfiler

intModules = as.vector(names(intModulesList))
pval = 0.05
qval = 0.05
# enrich with ENTREZID only (GENENAME not supported)
kkModulesList <- list()
for (i in 1:length(IQRintModulesList)) {
  # go enrichment for each module of interest
  geneNames <- IQRintModulesList[[i]]
  universe <- names(geneList)
  kk <- enrichKEGG(gene = geneNames, universe = universe, 
                   organism = 'api', pAdjustMethod = "fdr", pvalueCutoff = pval, 
                   qvalueCutoff = qval)
  objectName <- paste0(intModules[i], "KEGG")
  kkModulesList[[objectName]] <- kk
}
sapply(kkModulesList, nrow)



## packages for analyses requiring geneList object
library("DOSE")
library("ggplot2")
library("ggraph")
library("igraph")
library("enrichplot")
library("pathview")
library("KEGGREST")
library("png")

## cnetplot for GO terms

# generate cnetplots with ENTREZID geneLists
# open and run cnetplotTES_complete.R in workingDir to load necessary functions
source("cnetplotTES_complete.R")
# generate a folder to place the cnetplot output
folderName <- paste0("cnetplotTES_Ap_gene_LRT_noLSR1")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName))
}
# now use cnetplot
GOcnetplotResultsList <- list()
for (i in 1:length(egoModulesSimp)) {
  if ( nrow(egoModulesSimp[[i]]) > 0) {
    c <- cnetplotTES(egoModulesSimp[[i]], showCategory = nrow(egoModulesSimp[[i]]), 
                     foldChange = geneList,
                     limits = c(-1,1) )
    sublistName <- paste0(names(egoModulesSimp)[i])
    GOcnetplotResultsList[[sublistName]] <- c
    jpeg(paste0("cnetplotTES_ENTREZID_", names(egoModulesSimp)[i], 
                "_gene_LRT_noLSR1.jpg"), 
         width = 720, height = 720, pointsize = 20)
    print(c)
    dev.off()
  }
}

## cnetplot for KEGG terms, available only with ENTREZID
KEGGcnetplotResultsList <- list()
for (i in 1:length(kkModulesList)) {
  if ( nrow(kkModulesList[[i]]) > 0) {
    #moduleName <- gsub("[[:upper:]]", "", names(kkModulesList)[i] )
   # bookmark <- grep(moduleName, names(geneList))
    c <- cnetplotTES(kkModulesList[[i]], showCategory = nrow(kkModulesList[[i]]), 
                     foldChange = geneList,
                     limits = c(-1,1) )
    sublistName <- paste0(names(kkModulesList)[i])
    KEGGcnetplotResultsList[[sublistName]] <- c
    jpeg(paste0("cnetplotTES_KEGG_", names(kkModulesList)[i], "_gene_LRT_noLSR1.jpg"), 
         width = 720, height = 720, pointsize = 20)
    print(c)
    dev.off()
  }
}

## KEGG pathview
# pull out KEGG pathway identifiers and do pathview on each one for each module of interest
# open and run pathviewTES_complete.R in workingDir to load necessary functions
source("pathviewTES_complete.R")
pathviewResultsList <- list()
for (i in 1:length(kkModulesList)) {
  pathways <- kkModulesList[[i]]$ID
  moduleName <- gsub("[[:upper:]]", "", names(kkModulesList)[i] )
  bookmark <- IQRintModulesList[[grep(moduleName, names(IQRintModulesList))]]
  setwd(wkdir)
  folderName <- paste0("pathviewTES_", names(kkModulesList)[i], "_gene_LRT_noLSR1")
  if (file.exists(folderName)){
    setwd(file.path(wkdir, folderName))
  } else {
    dir.create(file.path(wkdir, folderName))
    setwd(file.path(wkdir, folderName))
  }
  for (j in 1: length(pathways)) {
    k <- pathviewTES(gene.data = geneList[names(geneList) %in% bookmark], 
                     pathway.id = pathways[j], species = "api", 
                     gene.idtype = "Entrez", limit = list(gene= c(-1,1), cpd = 1))
    resultsName <- paste0(intModules[i], ".", pathways[j])
    pathviewResultsList[[resultsName]] <- k
  }
}
setwd(wkdir)

# rename and save enrichment results
WGCNAgoListAp <- egoModulesSimp
WGCNAkeggListAp <- kkModulesList
WGCNAcnetplotGoListAp <- GOcnetplotResultsList
WGCNAcnetplotKeggListAp <- KEGGcnetplotResultsList
WGCNApathviewListAp <- pathviewResultsList


# save results
save(WGCNAgoListAp, WGCNAkeggListAp, WGCNAcnetplotGoListAp, 
     WGCNAcnetplotKeggListAp, WGCNApathviewListAp,
     file = paste0("FunctionalEnrichment_Ap_rld_noLSR1.RData"))
  
