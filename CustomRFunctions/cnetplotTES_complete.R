# created by Tom Smith on Feb 27, 2019
## modified from https://github.com/GuangchuangYu/enrichplot/blob/master/R/cnetplot.R

## cnetplotTES
##' @rdname cnetplotTES
##' @exportMethod cnetplotTES

setGeneric("cnetplotTES", function(x, showCategory = 5, foldChange = NULL, layout = "kk", ...)
standardGeneric("cnetplotTES") )

setMethod("cnetplotTES", signature(x = "enrichResult"),
          function(x, showCategory = 5,
                   foldChange = NULL, layout = "kk", ...) {
            cnetplotTES.enrichResult(x, showCategory = showCategory,
                                   foldChange = foldChange, layout = layout,
                                   ...)
          })


## fc_paletteTES
fc_paletteTES <- function (fc) 
{
  if (all(fc > 0, na.rm = TRUE)) {
    palette <- color_palette(c("blue", "red"))
  }
  else if (all(fc < 0, na.rm = TRUE)) {
    palette <- color_palette(c("green", "blue"))
  }
  else {
    # change the color scheme to something perfect!!!!
    palette <- color_palette(c("cyan", "blue", "grey60", "red", "yellow"))
  }
  return(palette)
}



## fc_readable, unmodified
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



# list2graph, unmodified
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


## extract_geneSets, unmodified
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


## update_n, unmodified
update_n <- function (x, showCategory) 
{
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  return(n)
}


##' @rdname cnetplotTES
##' @param colorEdge whether coloring edge by enriched terms
##' @param circular whether using circular layout
##' @param node_label whether display node label
##' @importFrom ggraph geom_edge_arc
##' @author Guangchuang Yu, modified by Thomas E Smith
## cnetplotTES.enrichResult
cnetplotTES.enrichResult <- function(x,
                                   showCategory = 5,
                                   foldChange   = NULL,
                                   layout = "kk",
                                   colorEdge = FALSE,
                                   circular = FALSE,
                                   node_label = TRUE,
                                   ...) {
  
  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  } else {
    geom_edge <- geom_edge_link   # library(ggraph)
  }
  
  geneSets <- extract_geneSets(x, showCategory)   # returns gene lists for the number of enriched GO terms as defined by showCategory
  
  g <- list2graph(geneSets)   # converts GO terms and associated gene lists to an igraph graph, containing edge and vertices info
  
  foldChange <- fc_readable(x, foldChange)    # if foldChange is set to readable, this swaps out entrez ID for gene symbol
  
  size <- sapply(geneSets, length)   # obtains length (number of genes) per GO term list
  V(g)$size <- min(size)/2          # set size of vertices of igraph graph to half the length of the smallest gene group
  
  n <- length(geneSets)         # number of geneLists/GO term groups
  V(g)$size[1:n] <- size        # set size of first n vertices to a size corresponding to how many genes are in the GO term gene groups
  # first n vertices represent those of the GO term nodes
  
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))  # set edge categories of igraph graph as the GO term group names
    edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)    # for geom_edge, defined earlier, set aesthetics to ??? ...
  } else {
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')   # ... or if no color = ~category, then set all to same color
  }
  
  if (!is.null(foldChange)) {        # if foldChange is not NULL ...
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]    # ... then filter foldchange, containing info for all genes in this geneList, by the vertice names of the GO enriched genes
    V(g)$color <- NA          # set up a column for color of vertices
    V(g)$color[(n+1):length(V(g))] <- fc    # set vertice color to foldchange of each GO-enriched gene
    palette <- fc_paletteTES(fc)              # select color palette, I changed from green-red to red-grey palette
    p <- ggraph(g, layout=layout, circular = circular) +
      edge_layer +
      geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
      scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494", ... )
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph(g, layout=layout, circular=circular) +
      edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size))
  }
  
  p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
    theme_void()
  
  if (node_label)
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
  
  return(p)
}
