list2graph <- function (inputList) 
{
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed = FALSE)   # library("igraph"), creates an igraph graph from a data frame. graph contains edge list and edge/vertex attributes
  return(g)
}



list2df <- function (inputList) 
{
  ldf <- lapply(1:length(inputList), function(i) {           # make a list, separated by each GO term, in which each gene ID has it's GO term explicitly assigned to it
    data.frame(categoryID = rep(names(inputList[i]), length(inputList[[i]])), 
               Gene = inputList[[i]])
  })
  do.call("rbind", ldf)    # combine each list into a single data frame
}


## 2-25-19
## me playing with these functions

inputList <- geneSets
length(inputList)
inputList <- list2df(inputList)
dim(inputList)
g <- list2graph(inputList)
g



# what it does: converts GO terms and associated gene lists to an igraph graph, containing edge and vertices info


