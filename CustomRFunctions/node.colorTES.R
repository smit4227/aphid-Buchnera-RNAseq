##node.colorTES
node.colorTES <- function (plot.data = NULL, discrete = FALSE, limit = 1, 
                           bins = 20, both.dirs = TRUE, 
                           low = "cyan", midlow = "blue", mid = "gray60", midhigh = "red", high = "yellow", 
                           na.col = "transparent", trans.fun = NULL) 
{
  if (is.null(plot.data)) 
    return(NULL)
  node.summary = plot.data[, -c(1:8)]
  if (length(dim(node.summary)) == 2) {
    node.summary = as.matrix(node.summary)
  }
  else names(node.summary) = rownames(plot.data)
  if (!is.null(trans.fun)) 
    node.summary = trans.fun(node.summary)
  if (both.dirs & length(limit) == 1) {
    limit = c(-abs(limit), abs(limit))
  }
  else if (length(limit) == 1) {
    limit = c(0, limit)
  }
  disc.cond1 = all(as.integer(limit) == limit)
  disc.cond2 = (limit[2] - limit[1])%%bins == 0
  if (discrete & disc.cond1 & disc.cond2) {
    node.summary[] = as.integer(node.summary)
    limit[2] = limit[2] + 1
    bins = bins + 1
  }
  else if (discrete) {
    message("Note: ", "limit or bins not proper, data not treated as discrete!")
  }
  node.summary[node.summary > limit[2]] = limit[2]
  node.summary[node.summary < limit[1]] = limit[1]
  if (both.dirs) {
    cols = colorpanel2(bins, low = low, midlow = midlow, mid = mid, midhigh = midhigh, high = high)
  }
  else cols = colorpanel2(bins, low = mid, mid = midhigh, high = high)
  na.col = colorpanel2(1, low = na.col, midlow = na.col, mid = na.col, midhigh = na.col, high = na.col)
  data.cuts = seq(from = limit[1], to = limit[2], length = bins + 
                    1)
  index.ts = cols.ts = node.summary
  index.ts[] = cut(node.summary, data.cuts, include.lowest = TRUE, 
                   right = F)
  cols.ts[] = cols[index.ts]
  cols.ts[is.na(cols.ts)] = na.col
  return(cols.ts)
}


## colorpanel2
colorpanel2 <- function (n, low, midlow, mid, midhigh, high) 
{
  if (missing(mid) || missing(high)) {
    low <- col2rgb(low)
    if (missing(high)) 
      high <- col2rgb(mid)
    else high <- col2rgb(high)
    red <- seq(low[1, 1], high[1, 1], length = n)/255
    green <- seq(low[3, 1], high[3, 1], length = n)/255
    blue <- seq(low[2, 1], high[2, 1], length = n)/255
  }
  else {
    isodd <- n%%2 == 1
    if (isodd) {
      n <- n + 1
    }
    low <- col2rgb(low)
    midlow <- col2rgb(midlow)
    mid <- col2rgb(mid)
    midhigh <- col2rgb(midhigh)
    high <- col2rgb(high)
    lower <- floor(n/4)
    midlower <- floor(n/4)
    midupper <- (n - lower - midlower)/2
    upper <- n - lower - midlower - midupper
    red <- c(seq(low[1, 1], midlow[1, 1], length = lower),
             seq(midlow[1, 1], mid[1, 1], length = midlower),
             seq(mid[1, 1], midhigh[1, 1], length = midupper),
             seq(midhigh[1, 1], high[1, 1], length = upper))/255
    green <- c(seq(low[3, 1], midlow[3, 1], length = lower),
               seq(midlow[3, 1], mid[3, 1], length = midlower),
               seq(mid[3, 1], midhigh[3, 1], length = midupper),
               seq(midhigh[3, 1], high[3, 1], length = upper))/255
    blue <- c(seq(low[2, 1], midlow[2, 1], length = lower),
              seq(midlow[2, 1], mid[2, 1], length = midlower),
              seq(mid[2, 1], midhigh[2, 1], length = midupper),
              seq(midhigh[2, 1], high[2, 1], length = upper))/255
    if (isodd) {
      red <- red[-(lower + 1)]
      green <- green[-(lower + 1)]
      blue <- blue[-(lower + 1)]
    }
  }
  rgb(red, blue, green)
}