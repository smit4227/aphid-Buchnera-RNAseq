# R script for calculation of differential gene expression using DESeq2
# written by Thomas E. Smith Jan, 2019 
# last updated Sep 2019

# set the stage for A. pisum samples
library(DESeq2)
wkdir <- "/Users/User/Work"
datadir <- "/Users/User/Data"
setwd(wkdir)

## first analyze A. pisum gene counts

# generate a sampleTable compatible with DESeq2
directory <- datadir
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
                   "TES-21_cat_trim_filt_adapt_hisat2_Ap_gene.gff")
sampleNames <- c("AL1-1", "AL1-2", "AL1-3", "AL3-1", "AL3-2", "AL3-3", 
                 "AL4-1", "AL4-2", "AL4-3", "TL1-1", "TL1-2", "TL1-3", 
                 "TL2-1", "TL2-2", "TL2-3", "TUCc-1", "TUCc-2", "TUCc-3", 
                 "AUSc-1", "AUSc-2", "AUSc-3")
sampleTiter <- as.factor(c(42, 42, 42, 82, 82, 82, 62, 62, 62, 25, 25, 25, 
                           15, 15, 15, 40, 40, 40, 35, 35, 35))
sampleColor <- as.factor(c("green", "green", "green", "red", "red", "red", 
                           "red", "red", "red", "green", "green", "green", 
                           "green", "green", "green", "green", "green", "green", 
                           "green", "green", "green"))
sampleParentage <- as.factor(c("A", "A", "A", "A", "A", "A", "A", "A", "A", 
                               "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                               "A", "A", "A"))
sampleTableAp <- data.frame(sampleName=sampleNames, fileName=sampleFilesAp, 
                            titer=sampleTiter, color=sampleColor, 
                            parentage=sampleParentage)

# use deseq2 to calculate differential expression using Likelihood Ratio Test
ddsHTSeqAp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableAp, 
                                         directory = directory, design = ~ titer)
ddsAp <- DESeq(ddsHTSeqAp, test="LRT", reduced= ~1)
nrow(ddsAp)
# remove rows with 0 counts
keep <- rowSums(counts(ddsAp)) >= 1
ddsAp <- ddsAp[keep,]
nrow(ddsAp)

## now the same process for Buchnera

# generate a sampleTable compatible with DESeq2
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
                   "TES-21_cat_trim_filt_adapt_bwa_Ba_gene.gff")
sampleTableBa <- data.frame(sampleName=sampleNames, fileName=sampleFilesBa, 
                            titer=sampleTiter)
ddsHTSeqBa <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableBa, 
                                         directory = directory, design = ~ titer)
ddsBa <- DESeq(ddsHTSeqBa, test="LRT", reduced= ~1)
nrow(ddsBa)
keep <- rowSums(counts(ddsBa)) >= 1
ddsBa <- ddsBa[keep,]
nrow(ddsBa)

# summarize results in a table and order by fdr adjusted p-value
resAp <- results(ddsAp)
resOrderedAp <- resAp[order(resAp$padj),]
resBa <- results(ddsBa)
resOrderedBa <- resBa[order(resBa$padj),]

## we need some information about each gene to include in our summary table
## we will get it from the gene annotation file itself

# load the A. pisum gene annotation file as a table
gffAp <- read.table("GCF_000142985.2_Acyr_2.0_genomic.gff", sep="\t", quote="", 
                    header=FALSE, fill=FALSE)
gffNames <- c("seqname", "source", "feature", "start", "end", "score", "strand", 
              "frame", "attribute")
names(gffAp) <- gffNames
table <- gffAp
library(rapportools)

# acquire information on "gene" features
# create an empty table to fill with gene information
geneInfoAp <- data.frame(geneNumber = rep(NA, nrow(table)), 
                         locus = rep(NA, nrow(table))) 
# search
for (i in 1:nrow(table)) 
{
  if (table$feature[i] == "gene"){
    s <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out geneNumber info and add to geneInfo table
    id <- grep("ID=", s[[1]])
    ids <- gsub("ID=", "", s[[1]][id])
    if (!is.empty(ids)) {
      geneInfoAp$geneNumber[i] <- ids
    }
    
    # pull out gene locus info and add to geneInfo table
    l <- grep("Name=", s[[1]])
    loc <- gsub("Name=", "", s[[1]][l])
    if (!is.empty(loc)) {
      geneInfoAp$locus[i] <- loc
    }
  }
}

# acquire information on "RNA" and "trasncript" features
rnaInfoAp <- data.frame(geneNumber = rep(NA, nrow(table)), 
                        rnaNumber = rep(NA, nrow(table)), rnaType = rep(NA, nrow(table)), 
                        rnaName = rep(NA, nrow(table)), product = rep(NA, nrow(table))) 
for (i in 1:nrow(table)) 
{
  if (any(grep("RNA", table$feature[i]))) {
    s0 <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent geneNumber info and add to rnaInfo table
    p <- grep("Parent=", s0[[1]])
    parent <- gsub("Parent=", "", s0[[1]][p])
    if (!is.empty(parent)) {
      rnaInfoAp$geneNumber[i] <- parent
    }
    
    # pull out parent rnaNumber info and add to rnaInfo table
    r <- grep("ID=", s0[[1]])
    rna <- gsub("ID=", "", s0[[1]][r])
    if (!is.empty(rna)) {
      rnaInfoAp$rnaNumber[i] <- rna
    }
    
    rnaInfoAp$rnaType[i] <- as.character(table$feature[i])
    
    # pull out transcript accession number info and add to rnaInfo table
    a <- grep("Name=", s0[[1]])
    acc <- gsub("Name=", "", s0[[1]][a])
    if (!is.empty(acc)) {
      rnaInfoAp$rnaName[i] <- acc
    }
    
    # pull out transcript product info and add to rnaInfo table
    p <- grep("product=", s0[[1]])
    product <- gsub("product=", "", s0[[1]][p])
    if (!is.empty(product)) {
      rnaInfoAp$product[i] <- product
    }
  }
  
  if (any(grep("transcript", table$feature[i]))) {
    s0 <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent geneNumber info and add to rnaInfo table
    p <- grep("Parent=", s0[[1]])
    parent <- gsub("Parent=", "", s0[[1]][p])
    if (!is.empty(parent)) {
      rnaInfoAp$geneNumber[i] <- parent
    }
    
    # pull out parent rnaNumber info and add to rnaInfo table
    r <- grep("ID=", s0[[1]])
    rna <- gsub("ID=", "", s0[[1]][r])
    if (!is.empty(rna)) {
      rnaInfoAp$rnaNumber[i] <- rna
    }
    
    rnaInfoAp$rnaType[i] <- as.character(table$feature[i])
    
    # pull out transcript accession number info and add to rnaInfo table
    n <- grep("Name=", s0[[1]])
    name <- gsub("Name=", "", s0[[1]][n])
    if (!is.empty(name)) {
      rnaInfoAp$rnaName[i] <- name
    }
    
    # pull out transcript product info and add to rnaInfo table
    p <- grep("product=", s0[[1]])
    product <- gsub("product=", "", s0[[1]][p])
    if (!is.empty(product)) {
      rnaInfoAp$product[i] <- product
    }
  }
}

# acquire information on "CDS" features
protInfoAp <- data.frame(rnaNumber = rep(NA, nrow(table)), protName = rep(NA, nrow(table))) 
for (i in 1:nrow(table)) 
{
  if (table$feature[i] == "CDS"){
    s <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent RNA info and add to protInfo table
    p <- grep("Parent=", s[[1]])
    parent <- gsub("Parent=", "", s[[1]][p])
    if (!is.empty(parent)) {
      protInfoAp$rnaNumber[i] <- parent
    }
    
    # pull out protein accession number info and add to protInfo table
    n <- grep("Name=", s[[1]])
    name <- gsub("Name=", "", s[[1]][n])
    if (!is.empty(name)) {
      protInfoAp$protName[i] <- name
    }
  }
}


# clean up tables and merge
geneInfo <- geneInfoAp[complete.cases(geneInfoAp$geneNumber), ]
rnaInfo <- rnaInfoAp[complete.cases(rnaInfoAp$geneNumber), ]
gffInfoAp0 <- merge(geneInfo, rnaInfo, by = "geneNumber", all = TRUE)
any(is.na(gffInfoAp0$geneNumber))
# add "CDS" information
protInfo <- protInfoAp[complete.cases(protInfoAp$rnaNumber), ]
gffInfoAp <- merge(gffInfoAp0, protInfo, by = "rnaNumber", all = TRUE)
# organize
gffInfoAp <- gffInfoAp[order(gffInfoAp$geneNumber),]

# genes in htseq-count files are annotated by geneNumber
# match these genes in resOrdered to the info collected in gffInfo for all Ap genes 
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


## finally, include entrez identifiers for ease of comparison with downstream analyses
## most "locus" features are entrez IDs, but some are not
## match NCBI locus identifiers with Entrez IDs obtained from OrgDb

# obtain Entrez IDs from OrgDb
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, "pisum")
apisum <- ah[["AH72781"]]
keytypes(apisum)
apdb <- select(apisum, keys = keys(apisum), keytype = "ENTREZID", 
               columns = c("REFSEQ"))
apdb2 <- select(apisum, keys = keys(apisum, keytype = "SYMBOL"), keytype = "SYMBOL", 
                columns = c("ENTREZID"))
# create a data frame to gather Entrez IDs from NCBI geneNumbers 
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
# any duplicates after conversion to entrez?
nrow(num2entrezAp[duplicated(num2entrezAp$geneNumber), ])
nrow(num2entrezAp[duplicated(num2entrezAp$entrez), ])
# 90 genes with duplicate entrez IDs but distinct geneNumbers
duplicates <- num2entrezAp[duplicated(num2entrezAp$entrez), ]
length(unique(duplicates$geneNumber))
length(unique(duplicates$entrez))
# only 46 of these entrez IDs are unique
duplicates <- duplicates[!duplicated(duplicates$entrez), ]

# add entrez IDs to resOrderedAp
resOrderedAp = data.frame(entrez = num2entrezAp$entrez,
                          resOrderedAp)

# for each duplicate entrez ID, select the geneNumber with the 
# highest baseMean value (higher expression)
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

# save results in a .csv file
write.csv(as.data.frame(resOrderedAp0), 
          file = paste0("deseq2_Ap_gene_LRT_resOrderedAp_noLSR1.csv"),
          row.names = F)

# remove duplicated samples from ddsAp
ddsAp <- ddsAp[match(as.character(resOrderedAp0$geneNumber), row.names(ddsAp)), ]
# replace geneNumbers with entrez IDs for ddsAp
row.names(ddsAp) <- resOrderedAp0$entrez[match(as.character(resOrderedAp0$geneNumber),
                                               row.names(ddsAp))]





## now for Buchnera

gffBa <- read.table("Buchnera_genomic.gff", sep="\t", quote="", header=FALSE, 
                    fill=FALSE)
gffNames <- c("seqname", "source", "feature", "start", "end", "score", "strand", 
              "frame", "attribute")
names(gffBa) <- gffNames
table <- gffBa
library(rapportools)

# "gene" features
geneInfoBa <- data.frame(geneNumber = rep(NA, nrow(table)), 
                         geneSymbol = rep(NA, nrow(table)),
                         locus = rep(NA, nrow(table)), 
                         locus_tag = rep(NA, nrow(table)))
for (i in 1:nrow(table)) 
{
  if (table$feature[i] == "gene"){
    s <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out geneNumber info and add to geneInfo table
    id <- grep("ID=", s[[1]])
    ids <- gsub("ID=", "", s[[1]][id])
    if (!is.empty(ids)) {
      geneInfoBa$geneNumber[i] <- ids
    }
    
    # pull out gene symbol info and add to geneInfo table
    sym <- grep("Name=", s[[1]])
    symbol <- gsub("Name=", "", s[[1]][sym])
    if (!is.empty(symbol)) {
      geneInfoBa$geneSymbol[i] <- symbol
    }
    
    # pull out gene locus info and add to geneInfo table
    l <- grep("Dbxref=", s[[1]])
    loc <- gsub("Dbxref=GeneID:", "", s[[1]][l])
    if (!is.empty(loc)) {
      geneInfoBa$locus[i] <- loc
    }
    
    # pull out gene locus tag info and add to geneInfo table
    l <- grep("locus_tag=", s[[1]])
    loc <- gsub("locus_tag=", "", s[[1]][l])
    if (!is.empty(loc)) {
      geneInfoBa$locus_tag[i] <- loc
    }
  }
}

# "RNA" or "transcript" features
rnaInfoBa <- data.frame(geneNumber = rep(NA, nrow(table)), 
                        rnaNumber = rep(NA, nrow(table)), rnaType = rep(NA, nrow(table)), 
                        product = rep(NA, nrow(table))) 
for (i in 1:nrow(table)) 
{
  if (any(grep("RNA", table$feature[i]))) {
    s0 <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent geneNumber info and add to rnaInfo table
    p <- grep("Parent=", s0[[1]])
    parent <- gsub("Parent=", "", s0[[1]][p])
    if (!is.empty(parent)) {
      rnaInfoBa$geneNumber[i] <- parent
    }
    
    # pull out parent rnaNumber info and add to rnaInfo table
    r <- grep("ID=", s0[[1]])
    rna <- gsub("ID=", "", s0[[1]][r])
    if (!is.empty(rna)) {
      rnaInfoBa$rnaNumber[i] <- rna
    }
    
    rnaInfoBa$rnaType[i] <- as.character(table$feature[i])
    
    # pull out transcript product info and add to rnaInfo table
    p <- grep("product=", s0[[1]])
    product <- gsub("product=", "", s0[[1]][p])
    if (!is.empty(product)) {
      rnaInfoBa$product[i] <- product
    }
  }
  
  if (any(grep("transcript", table$feature[i]))) {
    s0 <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent geneNumber info and add to rnaInfo table
    p <- grep("Parent=", s0[[1]])
    parent <- gsub("Parent=", "", s0[[1]][p])
    if (!is.empty(parent)) {
      rnaInfoBa$geneNumber[i] <- parent
    }
    
    # pull out parent rnaNumber info and add to rnaInfo table
    r <- grep("ID=", s0[[1]])
    rna <- gsub("ID=", "", s0[[1]][r])
    if (!is.empty(rna)) {
      rnaInfoBa$rnaNumber[i] <- rna
    }
    
    rnaInfoBa$rnaType[i] <- as.character(table$feature[i])
    
    # pull out transcript product info and add to rnaInfo table
    p <- grep("product=", s0[[1]])
    product <- gsub("product=", "", s0[[1]][p])
    if (!is.empty(product)) {
      rnaInfoBa$product[i] <- product
    }
  }
}

# "CDS" features
protInfoBa <- data.frame(geneNumber = rep(NA, nrow(table)), protName = rep(NA, nrow(table)),
                         product = rep(NA, nrow(table))) 
for (i in 1:nrow(table)) 
{
  if (table$feature[i] == "CDS"){
    s <- strsplit(as.character(table$attribute[i]), split = ";")
    
    # pull out parent RNA info and add to protInfo table
    p <- grep("Parent=", s[[1]])
    parent <- gsub("Parent=", "", s[[1]][p])
    if (!is.empty(parent)) {
      protInfoBa$geneNumber[i] <- parent
    }
    
    # pull out protein accession number info and add to protInfo table
    n <- grep("Name=", s[[1]])
    name <- gsub("Name=", "", s[[1]][n])
    if (!is.empty(name)) {
      protInfoBa$protName[i] <- name
    }
    
    # pull out protein accession number info and add to protInfo table
    p <- grep("product=", s[[1]])
    product <- gsub("product=", "", s[[1]][p])
    if (!is.empty(product)) {
      protInfoBa$product[i] <- product
    }
  }
}

# are there the same number of genes detected in geneInfo as in rnaInfo and protInfo? 
# If TRUE, then yes
length(unique(geneInfoBa$geneNumber)) == (length(unique(rnaInfoBa$geneNumber)) + length(unique(protInfoBa$geneNumber)))
# how many genes detected in geneInfo are unaccounted for in either rnaInfo or 
# protInfo?
length(unique(geneInfoBa$geneNumber)) - (length(unique(rnaInfoBa$geneNumber)) + 
                                           length(unique(protInfoBa$geneNumber)))
# count rows with NA and remove empty rows
length(geneInfoBa$geneNumber[is.na(geneInfoBa$geneNumber)])
geneInfoBa <- geneInfoBa[complete.cases(geneInfoBa$geneNumber), ]
length(geneInfoBa$geneNumber) == length(unique(geneInfoBa$geneNumber))
length(rnaInfoBa$geneNumber[is.na(rnaInfoBa$geneNumber)])
rnaInfoBa <- rnaInfoBa[complete.cases(rnaInfoBa$geneNumber), ]
length(rnaInfoBa$geneNumber) == length(unique(rnaInfoBa$geneNumber))
length(protInfoBa$geneNumber[is.na(protInfoBa$geneNumber)])
protInfoBa <- protInfoBa[complete.cases(protInfoBa$geneNumber), ]
length(protInfoBa$geneNumber) == length(unique(protInfoBa$geneNumber))
# pool genes detected in geneInfo, rnaInfo, and protInfo
detected <- c(geneInfoBa$geneNumber, rnaInfoBa$geneNumber, protInfoBa$geneNumber)
rows <- grep(1, table(detected))
gene <- detected[rows]
rows
gene
# merge tables
gffInfoBa0 <- merge(geneInfoBa, rnaInfoBa, by = "geneNumber", all = TRUE)
any(is.na(gffInfoBa0$geneNumber))
# add protInfo
gffInfoBa <- merge(gffInfoBa0, protInfoBa, by = "geneNumber", all = TRUE)
# organize and remove duplicate rows
gffInfoBa <- gffInfoBa[order(gffInfoBa$geneNumber),]
gffInfoBa <- gffInfoBa[grep(TRUE, !duplicated(gffInfoBa$geneNumber)),]
# combine the two "product" columns of gffInfoBa
for (i in 1:length(gffInfoBa$product.y)){
  if (is.na(gffInfoBa$product.y[i])){
    gffInfoBa$product.y[i] <- gffInfoBa$product.x[i]
  }
}
gffInfoBa$product.x <- NULL
names(gffInfoBa)[length(names(gffInfoBa))] <- "product"

# # match these genes in resOrdered to the info collected in gffInfo for all Ba genes 
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
# luckily for us, entrez IDs for Buchnera are the same as "locus" features
resOrderedBa0 <- data.frame(entrez = resOrderedBa$locus,
                            resOrderedBa)
# save results and gene annotations
write.csv(as.data.frame(resOrderedBa0), 
          file = paste0("deseq2_Ba_gene_LRT_resOrderedBa_noLSR1.csv"),
          row.names = F)
# replace geneNumbers with entrez IDs for ddsBa
row.names(ddsBa) <- resOrderedBa0$entrez[match(row.names(ddsBa), 
                                               resOrderedBa0$geneNumber)]



## data transformations and saving results

# log-transformed, normalized counts are required for all downstream analyses
rldAp <- rlog(ddsAp, blind = F)
rldBa <- rlog(ddsBa, blind = F)
# save log transformed data 
write.csv(as.data.frame(assay(rldAp)), 
          file = paste0("RLD_deseq2_Ap_gene_LRT_rldAp_noLSR1.csv"))
write.csv(as.data.frame(assay(rldBa)), 
          file = paste0("RLD_deseq2_Ba_gene_LRT_rldBa_noLSR1.csv"))


# save non-transformed normalized counts as .csv files 
ddsAp0 <- estimateSizeFactors(ddsAp)
countsAp <- counts(ddsAp0, normalized=T)
write.csv(as.data.frame(countsAp), 
          file = paste0("counts_deseq2_Ap_gene_LRT_countsAp_noLSR1.csv"))
ddsBa0 <- estimateSizeFactors(ddsBa)
countsBa <- counts(ddsBa0, normalized=T)
write.csv(as.data.frame(countsBa), 
          file = paste0("counts_deseq2_Ba_gene_LRT_countsBa_noLSR1.csv"))

# save gene annotation files
write.csv(as.data.frame(gffInfoAp), file = paste0("Ap_gffInfo_ncbi.csv"), 
          row.names = F)
write.csv(as.data.frame(gffInfoBa), file = paste0("Ba_gffInfo_ncbi.csv"), 
          row.names = F)

# save objects to be used in later analyses
save(ddsAp, ddsBa, rldAp, rldBa, countsAp,
     countsBa, gffInfoAp, gffInfoBa, sampleTiter, 
     num2entrezAp,
     file = paste0("DESeq2ResultsAndGeneInfo_noLSR1.RData"))

