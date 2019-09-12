# R script for functional enrichment of correlation analysis gene modules
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"
setwd(wkdir)
lnames <- load(file = "CorrelationAnalysis_rld_noLSR1.RData")
lnames
lnames <- load(file = "EigengenesOfCorrelationAnalysis_rld_noLSR1.RData")
lnames


## GO enrichment of aphid correlation modules

# put aphid correlation module names into a vector
intModules = unique(rowClustCutCol)
# prepare an empty list to place entrez IDs for each module in
intModulesList <- list()
# pick out entrez numbers for every gene in each module of interest
for (module in intModules) {
  # select module probes
  modGenes = (rowClustCutCol==module)
  # get their entrez ID
  modGenes = names(rowClustCutCol)[modGenes]
  # save objects in a list
  intModulesList[[module]] <- modGenes
}

# GO enrichment with clusterProfiler
# clusterProfiler uses resources from AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, "pisum")
apisum <- ah[["AH72781"]]
library("clusterProfiler")
ontologies <- c("MF", "BP", "CC")
intModules = as.vector(names(intModulesList))
pval = 0.05
qval = 0.05
# enrich by ENTREZID as opposed to any other keytype
# Set argument readable = F because conversion to gene symbols is not 
# helpful for A. pisum,
# as there is little info availble under SYMBOL keytype.
egoModulesList <- list()
for (i in 1:length(intModulesList)) {
  # perform enrichment for each ontology
  for (j in 1:length(ontologies)) {
    # go enrichment for each module of interest
    geneNames <- intModulesList[[i]]
    universe <- as.factor(unlist(intModulesList))
    ego <- enrichGO(gene = geneNames, universe = universe, 
                    OrgDb = apisum, ont = ontologies[j], pAdjustMethod = "fdr", 
                    pvalueCutoff = pval, qvalueCutoff = qval, keyType = "ENTREZID")
    objectName <- paste0(intModules[i], ontologies[j])
    egoModulesList[[objectName]] <- ego
  }
}
sapply(egoModulesList, nrow)
# simplify GO enrichment results by combining similar/related GO terms
egoModulesSimp <- list()
Cutoff <- 0.6
for (i in 1:length(egoModulesList)) {
  ego <- egoModulesList[[i]]
  egoName <- paste0(names(egoModulesList)[i])
  egoSimp <- clusterProfiler::simplify(ego, cutoff = Cutoff, by = "p.adjust", select_fun = min)
  egoModulesSimp[[egoName]] <- egoSimp
}
sapply(egoModulesSimp, nrow)

# KEGG enrichment with clusterProfiler
pval = 0.05
qval = 0.05
# enrich with ENTREZID only (GENENAME not supported)
kkModulesList <- list()
for (i in 1:length(intModulesList)) {
  # go enrichment for each module of interest
  geneNames <- intModulesList[[i]]
  universe <- as.factor(unlist(intModulesList))
  kk <- enrichKEGG(gene = geneNames, universe = universe, 
                   organism = 'api', pAdjustMethod = "fdr", pvalueCutoff = pval, 
                   qvalueCutoff = qval)
  objectName <- paste0(intModules[i], "KEGG")
  kkModulesList[[objectName]] <- kk
}
sapply(kkModulesList, nrow)
# can't use simplify function for KEGG enrichment

## cnetplots for GO and KEGG
library("DOSE")
library("ggplot2")
library("ggraph")
library("igraph")
library("enrichplot")
library("pathview")
library("png")

# cnetplot for GO terms
# set up foldchange values with gene identifiers as rownames
fcAp <- resApFilt$log2FoldChange
names(fcAp) <- rownames(resApFilt)

# open and run cnetplotTES_complete.R in wkdir to load necessary functions
scriptdir <- "/Users/Jason/Tom Smith/Moran/RNA-Seq/RNA-Seq_final_data/Scripts/WGCNA/Functions/"
setwd(scriptdir)
source("cnetplotTES_complete.R")
setwd(wkdir)
# establish a new directory to place the resulting files
folderName <- paste0("cnetplot_simp_rld_Ap")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName)) }
# start plotting
for (i in 1:length(egoModulesSimp)) {
  if ( nrow(egoModulesSimp[[i]]) > 0) {
    moduleName <- gsub("[[:upper:]]", "", names(egoModulesSimp)[i])
    moduleGenes <- intModulesList[[moduleName]]
    bookmark <- names(fcAp) %in% moduleGenes
    c <- cnetplotTES(egoModulesSimp[[i]], #showCategory = 10,
                     showCategory = nrow(egoModulesSimp[[i]]), 
                     foldChange = fcAp[bookmark], # foldChange between AL3 and TL2
                     limits = c(-1,1) ) # anything outside of limits will be appear at -1 or 1
    jpeg(paste0("cnetplotTES_cor_Ap_", names(egoModulesSimp)[i], 
                "simp_gene_LRT_noLSR1.jpg"), 
         width = 720, height = 720, pointsize = 20)
    print(c)
    dev.off()
  }
}

# Problem with scale and gene node colors. 
# I think the problem occurs if my data contains no neg numbers (is unidirectional). 
# There is a conditional statement somewhere in "cnetplotTES_complete.R" 
# that needs changing.

# cnetplot for KEGG terms
# open and run keggview.nativeTES.R in wkdir to load necessary functions
setwd(scriptdir)
source("keggview.nativeTES.R")
setwd(wkdir)
# establish a new directory to place the resulting files
folderName <- paste0("cnetplot_KEGG_rld_Ap")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName)) }
# start plotting
for (i in 1:length(kkModulesList)) {
  if ( nrow(kkModulesList[[i]]) > 0) {
    moduleName <- gsub("[[:upper:]]", "", names(kkModulesList)[i] )
    moduleGenes <- intModulesList[[moduleName]]
    bookmark <- names(fcAp) %in% moduleGenes
    c <- cnetplotTES(kkModulesList[[i]], #showCategory = 10,
                     showCategory = nrow(kkModulesList[[i]]), 
                     foldChange = fcAp[bookmark], # foldChange between AL3 and TL2
                     limits = c(-1,1) ) # anything outside of limits will be appear at -1 or 1
    jpeg(paste0("cnetplotTES_KEGG_cor_Ap_", names(kkModulesList)[i], 
                "_gene_LRT_noLSR1.jpg"), 
         width = 720, height = 720, pointsize = 20)
    print(c)
    dev.off()
  }
}
# Similar problem with scale bar colors



#### KEGG enrichment for Buchnera correlation modules

# put Buchnera correlation module names into a vector
intModulesBa = unique(colClustCutCol)
# prepare an empty list to place entrez IDs for each module in
intModulesListBa <- list()
# pick out entrez numbers for every gene in each module of interest
for (module in intModulesBa) {
  # select module probes
  modGenes = (colClustCutCol==module)
  # get their entrez ID from num2entrezBa; enrichment uses "BUXXX" numbers for Buchnera
  modGenes = names(colClustCutCol)[modGenes]
  modLocustag = geneInfoBa$locus_tag[match(modGenes, geneInfoBa$entrez)]
  # save objects in a list
  intModulesListBa[[module]] <- modLocustag
}

# KEGG enrichment with clusterProfiler
pval = 0.05
qval = 0.05
# enrich with locus_tag only for Buchnera ("BUXXX") 
kkModulesListBa <- list()
for (i in 1:length(intModulesListBa)) {
  # go enrichment for each module of interest
  geneNames <- intModulesListBa[[i]]
  universe <- as.factor(unlist(intModulesListBa))
  kk <- enrichKEGG(gene = geneNames, universe = universe, 
                   organism = 'buc', pAdjustMethod = "fdr", pvalueCutoff = pval, 
                   qvalueCutoff = qval)
  objectName <- paste0(intModulesBa[i], "KEGG")
  kkModulesListBa[[objectName]] <- kk
}
sapply(kkModulesListBa, nrow)

# cnetplots for KEGG
# set up foldchange values with gene identifiers as rownames
fcBa <- resBaFilt$log2FoldChange
names(fcBa) <- resBaFilt$locus_tag

# establish a new directory to place the resulting files
folderName <- paste0("cnetplot_KEGG_rld_Ba")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName)) }
# start plotting
for (i in 1:length(kkModulesListBa)) {
  if ( nrow(kkModulesListBa[[i]]) > 0) {
    moduleName <- gsub("[[:upper:]]", "", names(kkModulesListBa)[i] )
    moduleGenes <- intModulesListBa[[moduleName]]
    bookmark <- names(fcBa) %in% moduleGenes
    c <- cnetplotTES(kkModulesListBa[[i]], #showCategory = 10,
                     showCategory = nrow(kkModulesListBa[[i]]), 
                     foldChange = fcBa[bookmark],
                     limits = c(-1,1) )
    jpeg(paste0("cnetplotTES_KEGG_cor_Ba_", names(kkModulesListBa)[i], 
                "_gene_LRT_noLSR1.jpg"), 
         width = 720, height = 720, pointsize = 20)
    print(c)
    dev.off()
  }
}
# similar problem with scale bar colors


## KEGG pathview on Ap and Ba correlated genes, using KEGG pathways of interest

# load the necessary packages
library("DOSE")
library("ggplot2")
library("ggraph")
library("igraph")
library("enrichplot")
library("pathview")
library("KEGGREST")
library("png")

# open and run pathviewTES_complete.R in wkdir to load necessary functions
setwd(scriptdir)
source("pathviewTES_complete.R")
setwd(wkdir)

# identify KEGG pathways of interest for A. pisum
pathwaysAp <- c("api01230", "api00010", "api00020", "api00030", "api00051", "api00520", 
                "api00190", "api00230", "api00240", "api00250", "api00260", "api00270", 
                "api00290", "api00310", "api00220", "api00330", "api00340", "api00400", 
                "api00670", "api00790", "api00760", "api04130", "api00970", "api02010", 
                "api04310", "api03020", "api03010", "api03008", "api03030", "api04013", 
                "api04142", "api04150", "api04330", "api04341", "api04350", "api04391")

# set up a folder name to place the resulting image files
folderName <- paste0("pathviewTES_Ap_rld_noLSR1")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName))
}
# start plotting
pathviewResultsListAp <- list()
for (j in 1:length(pathwaysAp)) {
  k <- pathviewTES(gene.data = fcAp , pathway.id = pathwaysAp[j], species = "api", 
                   gene.idtype = "Entrez",
                   limit = list(gene= c(-1,1), cpd = 1))
  resultsName <- paste0(pathwaysAp[j])
  pathviewResultsListAp[[resultsName]][j] <- k
}
setwd(wkdir)

# for Buchnera, pathview maps genes using entrez and not locus-tag, so re-make fcBa object
names(fcBa) <- resBaFilt$entrez

# identify KEGG pathways of interest for Buchnera
pathwaysBa <- c("buc01230", "buc00010", "buc00020", "buc00030", "buc00051", "buc00520", 
                "buc00190", "buc00230", "buc00240", "buc00250", "buc00260", "buc00270", 
                "buc00290", "buc00300", "buc00220", "buc00330", "buc00340", "buc00400", 
                "buc00670", "buc00790", "buc00760", "buc00550", "buc00970", "buc03010", 
                "buc03020", "buc03060", "buc04122", "buc03018", "buc03030", "buc02010", 
                "buc02060", "buc02040", "buc00410", "buc00480")
# set up a folder name to place the resulting image files
folderName <- paste0("pathviewTES_Ba_rld_noLSR1")
if (file.exists(folderName)){
  setwd(file.path(wkdir, folderName))
} else {
  dir.create(file.path(wkdir, folderName))
  setwd(file.path(wkdir, folderName))
}
# start plotting
pathviewResultsListBa <- list()
for (j in 1:length(pathwaysBa)) {
  k <- pathviewTES(gene.data = fcBa, pathway.id = pathwaysBa[j], species = "buc", 
                   gene.idtype = "entrez",
                   limit = list(gene= c(-1,1), cpd = 1))
  resultsName <- paste0(pathwaysBa[j])
  pathviewResultsListBa[[resultsName]][j] <- k
}

# rename some things
intModulesListAp <- intModulesList
egoModulesSimpAp <- egoModulesSimp
kkModulesListAp <- kkModulesList

# save objects to be used later
setwd(datadir)
save(intModulesListAp, intModulesListBa, pathwaysAp, pathwaysBa,
     egoModulesSimpAp, kkModulesListAp, kkModulesListBa, pathviewResultsListAp, pathviewResultsListBa,
     file = paste0("FunctionalEnrichmentOfCorrelationAnalysis_rld_noLSR1.RData"))

