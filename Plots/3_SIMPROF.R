# R script to test of significance of clustering of correlated genes
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"

# load correlation analysis data
lnames <- load("CorrelationAnalysis_rld_noLSR1.RData")

library("clustsig")

## are the patterns we see in hierarchical clustering of host-symbiont gene 
## correlations real, or could we get the same pattern by chance?
## we will test this using a similarity profiles test

## these analyses take a long time. Typically, people do num.expected = 1000
## and num.simulated = 999, but corMatrix is too huge

## first for Buchnera
# calculate p-values for dendrogram splits

colSimprof <- simprof(corMatrix, num.expected = 100, num.simulated = 99, 
                      sample.orientation = "column",
                      method.cluster = "complete", method.distance = "euclidean",
                      alpha = 0.01, silent = FALSE, increment = 100)
colSimprofPlot <- simprof.plot(colSimprof)


## then for A. pisum
# calculate p-values for dendrogram splits
rowSimprof <- simprof(corMatrix, num.expected = 100, num.simulated = 99, 
                      sample.orientation = "row",
                      method.cluster = "complete", method.distance = "euclidean",
                      alpha = 0.01, silent = FALSE, increment = 100)
rowSimprofPlot <- simprof.plot(rowSimprof)


# save objects to be used later
save(colSimprof, colSimprofPlot,
     rowSimprof, rowSimprofPlot,
     file = paste0("SIMPROF.RData"))





