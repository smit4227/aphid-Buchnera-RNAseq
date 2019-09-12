fc_palette <- function (fc) 
{
  if (all(fc > 0, na.rm = TRUE)) {
    palette <- color_palette(c("blue", "red"))
  }
  else if (all(fc < 0, na.rm = TRUE)) {
    palette <- color_palette(c("green", "blue"))
  }
  else {
    palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))
  }
  return(palette)
}




## 2-25-19
## me playing with these functions

x <- ego
fc <- foldchange[V(g)$name[(n+1):length(V(g))]]
length(fc_palette(fc))

library("RColorBrewer")
palette <- color_palette(c("darkgreen", "#0AFF34", "#B3B3B3", "#FF6347", "red"))   # original red-green palette
palette <- color_palette(brewer.pal(n=5, name = "RdGy"))    # my own palette, red to white to grey

palette <- color_palette(c("steelblue1", "steelblue4", "grey17", "indianred4", "indianred1"))  # my own palette, blue to black to red


