# R script for generating Ap vs Ba correlation plots
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"

### load DESeq2 data without LSR1 to obtain the gffInfo objects
setwd(wkdir)
rldAp <- read.csv("RLD_deseq2_Ap_gene_LRT_rldAp_noLSR1.csv", header = T)
row.names(rldAp) <- rldAp[,1]
rldAp[,1] <- NULL
rldBa <- read.csv("RLD_deseq2_Ba_gene_LRT_rldBa_noLSR1.csv", header = T)
row.names(rldBa) <- rldBa[,1]
rldBa[,1] <- NULL

## plot Aphid versus Buchnera gene expression

# which genes are you interested in?
geneAp <- "100163677"
geneBa <- "1109512"

# extract rld data for each and remove outliers
geneRldAp <- as.vector(rldAp[rownames(rldAp) == geneAp, ])
geneRldAp <- geneRldAp[-(grep("AL4.2", names(geneRldAp))) ]
geneRldAp <- geneRldAp[-(grep("TUCc.3", names(geneRldAp))) ]
geneRldBa <- as.vector(rldBa[rownames(rldBa) == geneBa, ])
geneRldBa <- geneRldBa[-(grep("AL4.2", names(geneRldBa))) ]
geneRldBa <- geneRldBa[-(grep("TUCc.3", names(geneRldBa))) ]

geneRldAp0 <- as.numeric(geneRldAp)
geneRldBa0 <- as.numeric(geneRldBa)

# correlate gene expression
geneCor <- cor.test(geneRldAp0, geneRldBa0, method = "pearson")

cor(geneRldAp, geneRldBa, method = "pearson")

### scatter plot
library("ggplot2")

# place data in a data table
geneDat <- data.frame(Ap = geneRldAp0, Ba = geneRldBa0, 
                      genotype = gsub("\\.[[:digit:]]", "", names(geneRldAp)) )


q <- ggplot(geneDat, aes(x = Ap, y = Ba, shape = genotype)) +
  coord_cartesian(ylim = c(13.25, 14), xlim = c(8, 9.5)) +
  geom_point(size = 4, stroke = 2, show.legend = F) +
  scale_shape_manual(values = c(2, 0, 6, 1, 17, 15, 16)) +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 2, colour = "black"),
        axis.ticks.length = unit(7, "pt"),
        axis.line = element_line(color = "black", size = 2),
        axis.text = element_text(color= "black", size = 15), 
        axis.title = element_blank() )
pdf(paste0("ExampleGeneCorrelation.pdf"), 
    width = 4, height = 4)
plot(q)
dev.off()

