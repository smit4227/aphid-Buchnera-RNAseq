extract_geneSets <- function (x, n) 
{
  n <- update_n(x, n)
  geneSets <- geneInCategory(x)  # library("DOSE")  # geneSets is a list of all gene IDS for every GO term that they represent
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]  # filter geneSets by those GO terms enriched in ego
  names(geneSets) <- y$Description # replace geneSets list names with GO term description (previously the GO ID)
  if (is.numeric(n)) {        # filter geneSets further by the top n enriched categories
    return(geneSets[1:n])
  }
  return(geneSets[n])    
}




## 2-25-19
## me playing with this function

x <- ego
n <- 5
names(extract_geneSets(x, n))

x <- ego
n <- 10
names(extract_geneSets(x, n))

x <- ego[1:10]
n <- 5
names(extract_geneSets(x, n))


# what it does: returns gene lists for top "n" enriched GO terms
