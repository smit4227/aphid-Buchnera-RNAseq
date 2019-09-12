# R script for generating publication quality PCAs from DESeq2 data
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"
datadir <- "/Users/User/Data"

### load DESeq2 data without LSR1 to obtain the gffInfo objects
lnames <- load("DESeq2ResultsAndGeneInfo_noLSR1.RData")
lnames

## we want to generate PCAs with and without LSR1 samples included, for comparison


####### perform DESeq2 with LSR1 included

# load DESeq2 and load data files
library("DESeq2")
setwd(datadir)

sampleFilesAp <- c("TES-1_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-2_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-3_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-4_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-5_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-6_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-7_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-8_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-9_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-10_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-11_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-12_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-13_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-14_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-15_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-16_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-17_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-18_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-19_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-20_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-21_cat_trim_filt_adapt_hisat2_Ap_gene.gff",
                   "TES-22_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-23_cat_trim_filt_adapt_hisat2_Ap_gene.gff", 
                   "TES-24_cat_trim_filt_adapt_hisat2_Ap_gene.gff")
sampleFilesBa <- c("TES-1_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-2_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-3_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-4_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-5_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-6_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-7_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-8_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-9_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-10_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-11_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-12_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-13_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-14_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-15_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-16_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-17_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-18_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-19_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-20_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-21_cat_trim_filt_adapt_bwa_Ba_gene.gff",
                   "TES-22_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-23_cat_trim_filt_adapt_bwa_Ba_gene.gff", 
                   "TES-24_cat_trim_filt_adapt_bwa_Ba_gene.gff")
sampleNames <- c("AL1-1", "AL1-2", "AL1-3", "AL3-1", "AL3-2", "AL3-3", 
                 "AL4-1", "AL4-2", "AL4-3", "TL1-1", "TL1-2", "TL1-3", 
                 "TL2-1", "TL2-2", "TL2-3", "TUCc-1", "TUCc-2", "TUCc-3", 
                 "AUSc-1", "AUSc-2", "AUSc-3", "LSR1c-1", "LSR1c-2", "LSR1c-3")
sampleTiter <- as.factor(c(42, 42, 42, 82, 82, 82, 62, 62, 62, 25, 25, 25, 
                           15, 15, 15, 40, 40, 40, 35, 35, 35, 50, 50, 50))
sampleColor <- as.factor(c("green", "green", "green", "red", "red", "red",
                           "red", "red", "red", "green", "green", "green",
                           "green", "green", "green", "green", "green", "green",
                           "green", "green", "green", "red", "red", "red"))
sampleParentage <- as.factor(c("A", "A", "A", "A", "A", "A", "A", "A", "A",
                               "T", "T", "T", "T", "T", "T", "T", "T", "T",
                               "A", "A", "A", "L", "L", "L"))
sampleTableAp <- data.frame(sampleName=sampleNames, fileName=sampleFilesAp, 
                            titer=sampleTiter, color=sampleColor, 
                            parentage=sampleParentage)
sampleTableBa <- data.frame(sampleName=sampleNames, fileName=sampleFilesBa, 
                            titer=sampleTiter, color=sampleColor, 
                            parentage=sampleParentage)

# use deseq2 to calculate differential expression using Likelihood Ratio Test
ddsHTSeqAp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableAp, directory = datadir, design = ~ titer)
ddsAp <- DESeq(ddsHTSeqAp, test="LRT", reduced= ~1)
nrow(ddsAp)
# remove rows with 0 counts
keep <- rowSums(counts(ddsAp)) >= 1
ddsAp <- ddsAp[keep,]
nrow(ddsAp)
# now for Buchnera
ddsHTSeqBa <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableBa, directory = datadir, design = ~ titer)
ddsBa <- DESeq(ddsHTSeqBa, test="LRT", reduced= ~1)
nrow(ddsBa)
keep <- rowSums(counts(ddsBa)) >= 1
ddsBa <- ddsBa[keep,]
nrow(ddsBa)

# summarize results in a table and order by adjusted p-value
resAp <- results(ddsAp)
resOrderedAp <- resAp[order(resAp$padj),]
resBa <- results(ddsBa)
resOrderedBa <- resBa[order(resBa$padj),]

# match the genes identified in RNA-Seq data of resOrdered to the info stored in gffInfo for all Ap genes 
probesAp <- row.names(resOrderedAp)
probes2gffInfoAp = match(probesAp, gffInfoAp$geneNumber)
resOrderedAp = data.frame(geneNumber = probesAp,
                          locus = gffInfoAp$locus[probes2gffInfoAp],
                          rnaNumber = gffInfoAp$rnaNumber[probes2gffInfoAp],
                          rnaType = gffInfoAp$rnaType[probes2gffInfoAp],
                          rnaName = gffInfoAp$rnaName[probes2gffInfoAp],
                          product = gffInfoAp$product[probes2gffInfoAp],
                          protName = gffInfoAp$protName[probes2gffInfoAp],
                          resOrderedAp)

# for A. pisum, replace gene numbers with entrez IDs using num2entrezAp object
# this is to eliminate duplicate genes that may have an affect on PCA distribution
# match NCBI locus identifiers with Entrez IDs from OrgDb
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, "pisum")
apisum <- ah[["AH72781"]]
keytypes(apisum)
apdb <- select(apisum, keys = keys(apisum), keytype = "ENTREZID", columns = c("REFSEQ"))
apdb2 <- select(apisum, keys = keys(apisum, keytype = "SYMBOL"), keytype = "SYMBOL", columns = c("ENTREZID"))
# create a data frame to put NCBI geneNumbers alongside entrez identifiers
num2entrezAp <- data.frame( geneNumber = rep(NA, nrow(resOrderedAp)),  
                            entrez = rep(NA, nrow(resOrderedAp)) )
# gather all of the entrez numbers from deseq2 data or OrgDb
for (i in 1:nrow(resOrderedAp) ) {
  geneNumber <- as.character(resOrderedAp$geneNumber[i])
  num2entrezAp$geneNumber[i] <- geneNumber
  locus <- as.character( resOrderedAp$locus[resOrderedAp$geneNumber == geneNumber] )
  rnaName <- as.character( resOrderedAp$rnaName[resOrderedAp$geneNumber == geneNumber] )
  if (grepl("LOC", locus) == T) {
    num2entrezAp$entrez[i] <- gsub("LOC", "", locus)
  }
  if (grepl("LOC", locus) == F) {
    if (!is.na(rnaName) == T) {
      probeRNA <- gsub("\\..*", "", rnaName)
      num2entrezAp$entrez[i] <- apdb$ENTREZID[grep(probeRNA, apdb$REFSEQ)]
    }
    if (is.na(rnaName) == T) {
      probeSymbol <- locus
      apdb2ID <- apdb2$ENTREZID[grep(probeSymbol, apdb2$SYMBOL)]
      if (length(apdb2ID) > 1) {
        num2entrezAp$entrez[i] <- apdb2ID[1]
      }  
      else num2entrezAp$entrez[i] <- apdb2ID
    } 
  }
} 
# identify duplicate entrez IDs within num2entrezAp
duplicates <- num2entrezAp[duplicated(num2entrezAp$entrez), ]
length(unique(duplicates$geneNumber))
length(unique(duplicates$entrez))
# only 49 of these entrez IDs are unique
duplicates <- duplicates[!duplicated(duplicates$entrez), ]

# add entrez IDs to resOrderedAp. Set these genes aside
resOrderedAp = data.frame(entrez = num2entrezAp$entrez,
                          resOrderedAp)

# for each duplicate entrez ID, select the geneNumber with the highest baseMean 
# value (higher expression)
dupInfo <- data.frame()
for (i in 1: length(duplicates$entrez)) { 
  dupEntrez <- as.character(duplicates$entrez[i])
  dupEntrezInfo <- resOrderedAp[grep(dupEntrez, resOrderedAp$entrez), ]
  # order by baseMean value, with 1 equal to highest value
  rowOrder <- order(dupEntrezInfo$baseMean, decreasing = T)
  # remove all rows containing the highest baseMean value
  entrezInfo <- dupEntrezInfo[(rowOrder != 1), ]
  dupInfo <- rbind(dupInfo, entrezInfo)
}
# remove these geneNumbers from resOrderedAp
resOrderedAp0 <- resOrderedAp[-(match(dupInfo$geneNumber, resOrderedAp$geneNumber)), ]
# are there equal numbers of unique entrez IDs and geneNumbers?
length(unique(resOrderedAp0$entrez)) == length(unique(resOrderedAp0$geneNumber))
# are there any entrez IDs or geneNumbers that are "NA"?
any(is.na(resOrderedAp0$entrez))
any(is.na(resOrderedAp0$geneNumber))

## save results in a .csv file
setwd(wkdir)
write.csv(as.data.frame(resOrderedAp0), 
          file = paste0("deseq2_Ap_gene_LRT_resOrderedAp_LSR1.csv"),
          row.names = F)

# remove duplicated samples from ddsAp
ddsAp <- ddsAp[match(as.character(resOrderedAp0$geneNumber), row.names(ddsAp)), ]
# replace geneNumbers with entrez IDs for ddsAp
row.names(ddsAp) <- resOrderedAp0$entrez[match(as.character(resOrderedAp0$geneNumber), row.names(ddsAp))]


### now for Buchnera

# match the genes identified in RNA-Seq data of resOrdered to the info stored in gffInfo for all Ap genes 
probesBa <- row.names(resOrderedBa)
probes2gffInfoBa = match(probesBa, gffInfoBa$geneNumber)
resOrderedBa = data.frame(geneNumber = probesBa,
                          geneSymbol = gffInfoBa$geneSymbol[probes2gffInfoBa],
                          locus = gffInfoBa$locus[probes2gffInfoBa],
                          locus_tag = gffInfoBa$locus_tag[probes2gffInfoBa],
                          rnaNumber = gffInfoBa$rnaNumber[probes2gffInfoBa],
                          rnaType = gffInfoBa$rnaType[probes2gffInfoBa],
                          protName = gffInfoBa$protName[probes2gffInfoBa],
                          product = gffInfoBa$product[probes2gffInfoBa],
                          resOrderedBa)
resOrderedBa0 <- data.frame(entrez = resOrderedBa$locus,
                            resOrderedBa)

## save results in a .csv file
write.csv(as.data.frame(resOrderedBa0), 
          file = paste0("deseq2_Ba_gene_LRT_resOrderedBa_LSR1.csv"),
          row.names = F)

# replace geneNumbers with entrez IDs for ddsBa
row.names(ddsBa) <- resOrderedBa0$entrez[match(row.names(ddsBa), resOrderedBa0$geneNumber)]



## transform the data

# log-transform data
rldAp <- rlog(ddsAp, blind = F)
rldBa <- rlog(ddsBa, blind = F)
# save results in a .csv file
write.csv(as.data.frame(assay(rldAp)), 
          file = paste0("RLD_deseq2_Ap_gene_LRT_rldAp_LSR1.csv"),
          row.names = T)
write.csv(as.data.frame(assay(rldBa)), 
          file = paste0("RLD_deseq2_Ba_gene_LRT_rldBa_LSR1.csv"),
          row.names = T)

# obtain counts
ddsAp0 <- estimateSizeFactors(ddsAp)
countsAp <- counts(ddsAp0, normalized=T)
ddsBa0 <- estimateSizeFactors(ddsBa)
countsBa <- counts(ddsBa0, normalized=T)

# rename objects to specify that this analysis was done with LSR1 samples included
resOrderedAp.LSR1 <- resOrderedAp0
resOrderedBa.LSR1 <- resOrderedBa0
ddsAp.LSR1 <- ddsAp
ddsBa.LSR1 <- ddsBa
rldAp.LSR1 <- rldAp
rldBa.LSR1 <- rldBa
countsAp.LSR1 <- countsAp
countsBa.LSR1 <- countsBa
sampleTiter.LSR1 <- sampleTiter

# save objects to be used in later analyses
save(resOrderedAp.LSR1, resOrderedBa.LSR1,
     ddsAp.LSR1, ddsBa.LSR1, 
     rldAp.LSR1, rldBa.LSR1, 
     countsAp.LSR1, countsBa.LSR1, 
     sampleTiter.LSR1, 
     file = paste0("DESeq2ResultsAndGeneInfo_LSR1.RData"))




### re-load DESeq2 rld data without LSR1 from .Rdata files
lnames <- load("DESeq2ResultsAndGeneInfo_noLSR1.RData")

# load ggplot2 tools necessary for PCA plotting
library("ggplot2")
library("ggfortify")

# plot PCAs
pcaAp.LSR1 <- plotPCA(rldAp.LSR1, intgroup = c("titer"), returnData = T)
pcaAp <- plotPCA(rldAp, intgroup = c("titer"), returnData = T)
pcaBa.LSR1 <- plotPCA(rldBa.LSR1, intgroup = c("titer"), returnData = T)
pcaBa <- plotPCA(rldBa, intgroup = c("titer"), returnData = T)

pcaList <- list(pcaAp.LSR1, pcaAp, pcaBa.LSR1, pcaBa)
pcaNames <- list("Ap_LSR1_", "Ap_noLSR1_", "Ba_LSR1_", "Ba_noLSR1")

setwd(wkdir)
for (i in 1:length(pcaList)) {
  pca <- pcaList[[i]]
  plotName <- pcaNames[i]
  # allow for different plot settings for Ap and Ba
  if (grepl("Ap", plotName)) {
    # label outlier point AL4-2
    outlier <- subset(pca, name == "AL4-2")
    # set unique cartesian coordinates for A. pisum plots, and label outlier
    pca_plot0 <- ggplot(pca, aes(PC1, PC2)) +
      coord_cartesian(ylim = c(-15, 20), xlim = c(-20, 25)) +
      geom_text(data = outlier, label = "AL4-2", size = 4, vjust = 3)
  }
  if (grepl("Ba", plotName)) {
    # label outlier point TUCc-3
    outlier <- subset(pca, name == "TUCc-3")
    # set unique cartesian coordinates for Buchnera plots, and label outlier
    pca_plot0 <- ggplot(pca, aes(PC1, PC2)) +
      coord_cartesian(ylim = c(-3, 3), xlim = c(-5, 12.5)) +
      geom_text(data = outlier, label = "Tucson-3", size = 4, vjust = 3)
  }
  # different shape order for whether LSR1 is present or not
  if (!grepl("no", plotName)) {
    pca_plot0 <- pca_plot0 +
      geom_point(aes(shape = titer), size = 4, stroke = 2, show.legend = F) +
      scale_shape_manual(values = c(15, 17, 1, 16, 2, 4, 6, 0))
  }
  if (grepl("no", plotName)) {
    pca_plot0 <- pca_plot0 +
      geom_point(aes(shape = titer), size = 4, stroke = 2, show.legend = F) +
      scale_shape_manual(values = c(15, 17, 1, 16, 2, 6, 0))
  }
  # plot common elements of each plot and save files
  pca_plot <- pca_plot0 + 
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(7, "pt"),
          axis.line = element_line(color = "black", size = 2),
          axis.text = element_text(color= "black", size = 15), 
          axis.title = element_blank() )
  pdf(paste0("PCA_deseq2LRT_", pcaNames[i], ".pdf"), 
      width = 4, height = 4)
  plot(pca_plot)
  dev.off()
  
  # plot legend alone
  pca_plot <- pca_plot0 +
    geom_point(aes(shape = titer), size = 4, stroke = 2, show.legend = T) +
    theme(plot.title = element_blank(),
          panel.background = element_blank(),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.key = element_rect(fill = NA) )
  pdf(paste0("PCA_deseq2LRT_legend_", ".pdf"), 
      width = 4, height = 4)
  plot(pca_plot)
  dev.off()
}


