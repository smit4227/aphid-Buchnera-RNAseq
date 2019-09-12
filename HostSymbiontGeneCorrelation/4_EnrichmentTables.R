# R script for producing summary tables of functional enrichment analysis
# written by Thomas E. Smith
# last updated by TES Sep 2019

# setting the stage
wkdir <- "/Users/User/Work"
setwd(wkdir)

# load necessary packages
library("DOSE")

# load data
setwd(datadir)
lnames <- load("FunctionalEnrichmentOfCorrelationAnalysis_rld_noLSR1.RData")
lnames

### compile GO terms and KEGG terms into a single data frame

# gather up all GO terms in egoModulesSimpAp, avoiding lists with no results
GOdf <- data.frame()
GOnames <- names(egoModulesSimpAp)
for (i in 1:length(GOnames)) {
  res <- egoModulesSimpAp[[ GOnames[i] ]][, 1:9]
  if (nrow(res) > 0) {
    # add module name and ontology to a data frame
    mod <- gsub(pattern = "[[:upper:]]", replacement = "", GOnames[i])
    ont <- gsub(pattern = "[[:lower:]]", replacement = "", GOnames[i])
    ont <- gsub(pattern = "[[:digit:]]", replacement = "", ont)
    modont <- data.frame(module = rep(mod, nrow(res)), ontology = rep(ont, nrow(res)))
    # combine res with modont
    res <- cbind(modont, res)
    # bind results of this iteration with those of previous
    GOdf <- rbind(GOdf, res)
  }
}


# add to this df KEGG terms
KEGGdf <- data.frame()
KEGGnames <- names(kkModulesListAp)
for (i in 1:length(KEGGnames)) {
  res <- kkModulesListAp[[ KEGGnames[i] ]][, 1:9]
  if (nrow(res) > 0) {
    # add module name and ontology to a data frame
    mod <- gsub(pattern = "[[:upper:]]", replacement = "", KEGGnames[i])
    ont <- gsub(pattern = "[[:lower:]]", replacement = "", KEGGnames[i])
    ont <- gsub(pattern = "[[:digit:]]", replacement = "", ont)
    modont <- data.frame(module = rep(mod, nrow(res)), ontology = rep(ont, nrow(res)))
    # combine res with modont
    res <- cbind(modont, res)
    # bind results of this iteration with those of previous
    KEGGdf <- rbind(KEGGdf, res)
  }
}

# merge GO and KEGG data frames
df <- rbind(GOdf, KEGGdf)

# filter by significance (p.adjust < 0.01)
dffilt <- df[df$p.adjust < 0.01, ]

# rename it
ApEnrichmentResults <- dffilt

# export the table
setwd(wkdir)
write.csv(ApEnrichmentResults, row.names = F, file = "GOandKEGGtermsAp.csv")


### now for Buchnera KEGG terms

KEGGdf <- data.frame()
KEGGnames <- names(kkModulesListBa)
for (i in 1:length(KEGGnames)) {
  res <- kkModulesListBa[[ KEGGnames[i] ]][, 1:9]
  if (nrow(res) > 0) {
    # add module name and ontology to a data frame
    mod <- gsub(pattern = "[[:upper:]]", replacement = "", KEGGnames[i])
    ont <- gsub(pattern = "[[:lower:]]", replacement = "", KEGGnames[i])
    ont <- gsub(pattern = "[[:digit:]]", replacement = "", ont)
    modont <- data.frame(module = rep(mod, nrow(res)), ontology = rep(ont, nrow(res)))
    # combine res with modont
    res <- cbind(modont, res)
    # bind results of this iteration with those of previous
    KEGGdf <- rbind(KEGGdf, res)
  }
}

# filter by significance (p.adjust < 0.01)
dffilt <- KEGGdf[KEGGdf$p.adjust < 0.01, ]

# rename it
BaEnrichmentResults <- dffilt

# export the table
setwd(wkdir)
write.csv(BaEnrichmentResults, row.names = F, file = "GOandKEGGtermsBa.csv")



### compile KEGG pathview results into a data frame to be exported

# first for A. pisum
dfAp <- data.frame()
for (i in 1:length(pathwaysAp)) {
  pathway <- pathwaysAp[i]
  res <- pathviewResultsListAp[[pathway]][[i]]
  if (is.null(res) == FALSE) {
    pathwaydf <- data.frame(pathway = rep(pathway, nrow(res)))
    # combine res with pathwaydf
    res <- cbind(pathwaydf, res)
    # bind results of this iteration with those of previous
    dfAp <- rbind(dfAp, res)
  }
}

# now for Buchnera
dfBa <- data.frame()
for (i in 1:length(pathwaysBa)) {
  pathway <- pathwaysBa[i]
  res <- pathviewResultsListBa[[pathway]][[i]]
  if (is.null(res) == FALSE) {
    pathwaydf <- data.frame(pathway = rep(pathway, nrow(res)))
    # combine res with pathwaydf
    res <- cbind(pathwaydf, res)
    # bind results of this iteration with those of previous
    dfBa <- rbind(dfBa, res)
  }
}


## save data frames as .csv files
write.csv(dfAp, file = paste0("KEGGpathviewTableAp.csv"), row.names = F)
write.csv(dfBa, file = paste0("KEGGpathviewTableBa.csv"), row.names = F)


