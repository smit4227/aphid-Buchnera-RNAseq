fc_readable <- function (x, foldChange = NULL) 
{
  if (is.null(foldChange)) 
    return(NULL)
  if (x@readable) {                      # x@readable is not true for x = ego, so fc_readable returns foldchange as the final output
    gid <- names(foldChange)             # extract gene IDs from names of foldchange
    if (is(x, "gseaResult")) {           # is x a gseaResult? if so, match IDs between gid and x@geneList
      ii <- gid %in% names(x@geneList)
    }
    else {
      ii <- gid %in% x@gene              # if x something else? if so, match IDs between gid and x@gene. This is true for x = ego, but x@readable = F
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]    # gid filtered by those matched in x@gene equal to x@gene2Symbol at same row
    names(foldChange) <- gid             # replacing names of foldchange with gene symbol
  }
  return(foldChange)
}



## 2-25-19
## me playing with these functions

x <- ego
head(fc_readable(x, foldChange = foldchange))
head(foldchange)

gid <- names(foldchange)
ii <- gid %in% x@gene
gid[ii] <- x@gene2Symbol[gid[ii]]


# what it does: if foldChange is set to readable, this swaps out entrez ID for symbol
