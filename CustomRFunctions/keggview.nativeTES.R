## keggview.nativeTES
keggview.nativeTES <- function (plot.data.gene = NULL, plot.data.cpd = NULL, cols.ts.gene = NULL, 
                             cols.ts.cpd = NULL, node.data, pathway.name, out.suffix = "pathview", 
                             kegg.dir = ".", multi.state = TRUE, match.data = TRUE, same.layer = TRUE, 
                             res = 300, cex = 0.25, discrete = list(gene = FALSE, cpd = FALSE), 
                             limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
                             both.dirs = list(gene = T, cpd = T), 
                             low = list(gene = "cyan", cpd = "thistle"),
                             midlow = list(gene = "blue", cpd = "purple"),  # added
                             mid = list(gene = "grey60", cpd = "gray60"), 
                             midhigh = list(gene = "red", cpd = "green"),   # added
                             high = list(gene = "yellow", cpd = "lawngreen"), 
                             na.col = "transparent", 
                             new.signature = TRUE, plot.col.key = TRUE, key.align = "x", 
                             key.pos = "topright", ...) 
{
  img <- readPNG(paste(kegg.dir, "/", pathway.name, ".png",     #library("png") 
                       sep = ""))
  width <- ncol(img)
  height <- nrow(img)
  cols.ts.gene = cbind(cols.ts.gene)
  cols.ts.cpd = cbind(cols.ts.cpd)
  nc.gene = max(ncol(cols.ts.gene), 0)
  nc.cpd = max(ncol(cols.ts.cpd), 0)
  nplots = max(nc.gene, nc.cpd)
  pn.suffix = colnames(cols.ts.gene)
  if (length(pn.suffix) < nc.cpd) 
    pn.suffix = colnames(cols.ts.cpd)
  if (length(pn.suffix) < nplots) 
    pn.suffix = 1:nplots
  if (length(pn.suffix) == 1) {
    pn.suffix = out.suffix
  }
  else pn.suffix = paste(out.suffix, pn.suffix, sep = ".")
  na.col = colorpanel2(1, low = na.col, high = na.col)
  if ((match.data | !multi.state) & nc.gene != nc.cpd) {
    if (nc.gene > nc.cpd & !is.null(cols.ts.cpd)) {
      na.mat = matrix(na.col, ncol = nplots - nc.cpd, nrow = nrow(cols.ts.cpd))
      cols.ts.cpd = cbind(cols.ts.cpd, na.mat)
    }
    if (nc.gene < nc.cpd & !is.null(cols.ts.gene)) {
      na.mat = matrix(na.col, ncol = nplots - nc.gene, 
                      nrow = nrow(cols.ts.gene))
      cols.ts.gene = cbind(cols.ts.gene, na.mat)
    }
    nc.gene = nc.cpd = nplots
  }
  out.fmt = "Working in directory %s"
  wdir = getwd()
  out.msg = sprintf(out.fmt, wdir)
  message("Info: ", out.msg)
  out.fmt = "Writing image file %s"
  multi.state = multi.state & nplots > 1
  if (multi.state) {
    nplots = 1
    pn.suffix = paste(out.suffix, "multi", sep = ".")
    if (nc.gene > 0) 
      cols.gene.plot = cols.ts.gene
    if (nc.cpd > 0) 
      cols.cpd.plot = cols.ts.cpd
  }
  for (np in 1:nplots) {
    img.file = paste(pathway.name, pn.suffix[np], "png", 
                     sep = ".")
    out.msg = sprintf(out.fmt, img.file)
    message("Info: ", out.msg)
    png(img.file, width = width, height = height, res = res)
    op = par(mar = c(0, 0, 0, 0))
    plot(c(0, width), c(0, height), type = "n", xlab = "", 
         ylab = "", xaxs = "i", yaxs = "i")
    if (new.signature) 
      img[height - 4:25, 17:137, 1:3] = 1
    if (same.layer != T) 
      rasterImage(img, 0, 0, width, height, interpolate = F)
    if (!is.null(cols.ts.gene) & nc.gene >= np) {
      if (!multi.state) 
        cols.gene.plot = cols.ts.gene[, np]
      if (same.layer != T) {
        render.kegg.node(plot.data.gene, cols.gene.plot, 
                         img, same.layer = same.layer, type = "gene", 
                         cex = cex)
      }
      else {
        img = render.kegg.node(plot.data.gene, cols.gene.plot, 
                               img, same.layer = same.layer, type = "gene")
      }
    }
    if (!is.null(cols.ts.cpd) & nc.cpd >= np) {
      if (!multi.state) 
        cols.cpd.plot = cols.ts.cpd[, np]
      if (same.layer != T) {
        render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                         img, same.layer = same.layer, type = "compound", 
                         cex = cex)
      }
      else {
        img = render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                               img, same.layer = same.layer, type = "compound")
      }
    }
    if (same.layer == T) 
      rasterImage(img, 0, 0, width, height, interpolate = F)
    pv.pars = list()
    pv.pars$gsizes = c(width = width, height = height)
    pv.pars$nsizes = c(46, 17)
    pv.pars$op = op
    pv.pars$key.cex = 2 * 72/res
    pv.pars$key.lwd = 1.2 * 72/res
    pv.pars$sign.cex = cex
    off.sets = c(x = 0, y = 0)
    align = "n"
    ucol.gene = unique(as.vector(cols.ts.gene))
    na.col.gene = ucol.gene %in% c(na.col, NA)
    if (plot.col.key & !is.null(cols.ts.gene) & !all(na.col.gene)) {
      off.sets = col.keyTES(limit = limit$gene, bins = bins$gene, 
                         both.dirs = both.dirs$gene, discrete = discrete$gene, 
                         graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                         key.pos = key.pos, cex = pv.pars$key.cex, lwd = pv.pars$key.lwd, 
                         low = low$gene, 
                         midlow = midlow$gene,    # added
                         mid = mid$gene, 
                         midhigh = midhigh$gene,  # added
                         high = high$gene, 
                         align = "n")
      align = key.align
    }
    ucol.cpd = unique(as.vector(cols.ts.cpd))
    na.col.cpd = ucol.cpd %in% c(na.col, NA)
    if (plot.col.key & !is.null(cols.ts.cpd) & !all(na.col.cpd)) {
      off.sets = col.keyTES(limit = limit$cpd, bins = bins$cpd, 
                         both.dirs = both.dirs$cpd, discrete = discrete$cpd, 
                         graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                         key.pos = key.pos, off.sets = off.sets, cex = pv.pars$key.cex, 
                         lwd = pv.pars$key.lwd, 
                         low = low$gene, 
                         midlow = midlow$gene,    # added
                         mid = mid$gene, 
                         midhigh = midhigh$gene,  # added
                         high = high$gene, 
                         align = align)
    }
    if (new.signature) 
      pathview.stamp(x = 17, y = 20, on.kegg = T, cex = pv.pars$sign.cex)
    par(pv.pars$op)
    dev.off()
  }
  return(invisible(pv.pars))
}



## col.keyTES
col.keyTES <- function (discrete = FALSE, limit = 1, bins = 10, cols = NULL, 
          both.dirs = TRUE, 
          low = "cyan", midlow = "blue", mid = "gray60", midhigh = "red", high = "yellow",     # added mid-levels
          graph.size, node.size, size.by.graph = TRUE, key.pos = "topright", 
          off.sets = c(x = 0, y = 0), align = "n", cex = 1, lwd = 1) 
{
  if (both.dirs & length(limit) == 1) {
    limit = c(-abs(limit), abs(limit))
  }
  else if (length(limit) == 1) {
    limit = c(0, limit)
  }
  disc.cond1 = all(as.integer(limit) == limit)
  disc.cond2 = (limit[2] - limit[1])%%bins == 0
  discrete = discrete & disc.cond1 & disc.cond2
  if (discrete) {
    limit[2] = limit[2] + 1
    bins = bins + 1
  }
  width = graph.size[1]
  height = graph.size[2]
  if (size.by.graph == T) {
    xs = width/80
    ys = height/40
  }
  else if (!missing(node.sizes)) {
    xs = node.size[1] * 3/bins
    ys = node.size[2]
  }
  else {
    message("Note: ", "color key not plotted, node.size is needed\n when size.by.graph=FALSE!")
    return(off.sets)
  }
  if (align == "x") {
    off.sets["x"] = 2 * xs
    off.sets["y"] = off.sets["y"] + 3 * ys
  }
  if (align == "y") 
    off.sets = off.sets + c(x = 3 * xs, y = 0)
  if (align == "n") 
    off.sets = off.sets + c(x = 2 * xs, y = 2 * ys)
  if (length(grep("right", key.pos)) == 1) {
    off.sets["x"] = off.sets["x"] + bins * xs
    x = width - off.sets["x"]
  }
  else {
    x = off.sets["x"]
    off.sets["x"] = off.sets["x"] + bins * xs
  }
  if (length(grep("top", key.pos)) == 1) 
    y = height - off.sets["y"]
  else y = off.sets["y"]
  ckx = seq(x, x + bins * xs, length = bins + 1)
  cky = c(y, y + ys)
  if (is.null(cols)) {
    if (both.dirs) {
      cols = colorpanel2(bins, low = low, midlow = midlow, mid = mid, midhigh = midhigh, high = high)    # added mid-levels
    }
    else cols = colorpanel2(bins, low = mid, mid = midhigh, high = high)   # added mid-levels
  }
  data.cuts = seq(from = limit[1], to = limit[2], length = bins + 
                    1)
  image(x = ckx, y = cky, z = cbind(data.cuts[-1]), col = cols, 
        axes = FALSE, add = T)
  if (!discrete) {
    label = format(data.cuts[c(1, bins/2 + 1, bins + 1)], 
                   digits = 2)
    text(x = seq(x, x + bins * xs, length = length(label)), 
         y = rep(y - ys, length(label)), label = label, cex = cex)
  }
  else {
    label = paste(as.integer(data.cuts[c(1, bins)]))
    text(x = seq(x, x + bins * xs, length = length(label)) + 
           c(xs, -xs)/2, y = rep(y - ys, length(label)), label = label, 
         cex = cex)
  }
  cky = c(y - 0.25 * ys, y + ys)
  for (i in 1:(bins + 1)) lines(rep(ckx[i], 2), cky, lwd = lwd)
  return(off.sets)
}

## render.kegg.node, unmodified
render.kegg.node <- function (plot.data, cols.ts, img, same.layer = TRUE, type = c("gene", 
                                                                                   "compound")[1], text.col = "black", cex = 0.25) 
{
  width = ncol(img)
  height = nrow(img)
  nn = nrow(plot.data)
  pwids = plot.data$width
  if (!all(pwids == max(pwids))) {
    message("Info: ", "some node width is different from others, and hence adjusted!")
    wc = table(pwids)
    pwids = plot.data$width = as.numeric(names(wc)[which.max(wc)])
  }
  if (type == "gene") {
    if (same.layer != T) {
      rect.out = sliced.shapes(plot.data$x + 0.5, height - 
                                 plot.data$y, plot.data$width/2 - 0.5, plot.data$height/2 - 
                                 0.25, cols = cols.ts, draw.border = F, shape = "rectangle")
      text(plot.data$x + 0.5, height - plot.data$y, labels = as.character(plot.data$labels), 
           cex = cex, col = text.col)
      return(invisible(1))
    }
    else {
      img2 = img
      pidx = cbind(ceiling(plot.data$x - plot.data$width/2) + 
                     1, floor(plot.data$x + plot.data$width/2) + 1, 
                   ceiling(plot.data$y - plot.data$height/2) + 1, 
                   floor(plot.data$y + plot.data$height/2) + 1)
      cols.ts = cbind(cols.ts)
      ns = ncol(cols.ts)
      brk.x = sapply(plot.data$width/2, function(wi) seq(-wi, 
                                                         wi, length = ns + 1))
      for (k in 1:ns) {
        col.rgb = col2rgb(cols.ts[, k])/255
        pxr = t(apply(pidx[, 1:2], 1, function(x) x[1]:x[2])) - 
          plot.data$x - 1
        sel = pxr >= ceiling(brk.x[k, ]) & pxr <= floor(brk.x[k + 
                                                                1, ])
        for (i in 1:nn) {
          sel.px = (pidx[i, 1]:pidx[i, 2])[sel[i, ]]
          node.rgb = img[pidx[i, 3]:pidx[i, 4], sel.px, 
                         1:3]
          node.rgb.sum = apply(node.rgb, c(1, 2), sum)
          blk.ind = which(node.rgb.sum == 0 | node.rgb.sum == 
                            1, arr.ind = T)
          node.rgb = array(col.rgb[, i], dim(node.rgb)[3:1])
          node.rgb = aperm(node.rgb, 3:1)
          for (j in 1:3) node.rgb[cbind(blk.ind, j)] = 0
          img2[pidx[i, 3]:pidx[i, 4], sel.px, 1:3] = node.rgb
        }
      }
      return(img2)
    }
  }
  else if (type == "compound") {
    if (same.layer != T) {
      nc.cols = ncol(cbind(cols.ts))
      if (nc.cols > 2) {
        na.cols = rep("#FFFFFF", nrow(plot.data))
        cir.out = sliced.shapes(plot.data$x, height - 
                                  plot.data$y, plot.data$width[1], plot.data$width[1], 
                                cols = na.cols, draw.border = F, shape = "ellipse", 
                                lwd = 0.2)
      }
      cir.out = sliced.shapes(plot.data$x, height - plot.data$y, 
                              plot.data$width[1], plot.data$width[1], cols = cols.ts, 
                              shape = "ellipse", blwd = 0.2)
      return(invisible(1))
    }
    else {
      blk = c(0, 0, 0)
      img2 = img
      w = ncol(img)
      h = nrow(img)
      cidx = rep(1:w, each = h)
      ridx = rep(1:h, w)
      pidx = lapply(1:nn, function(i) {
        ii = which((cidx - plot.data$x[i])^2 + (ridx - 
                                                  plot.data$y[i])^2 < (plot.data$width[i])^2)
        imat = cbind(cbind(ridx, cidx)[rep(ii, each = 3), 
                                       ], 1:3)
        imat[, 1:2] = imat[, 1:2] + 1
        ib = which(abs((cidx - plot.data$x[i])^2 + (ridx - 
                                                      plot.data$y[i])^2 - (plot.data$width[i])^2) <= 
                     8)
        ibmat = cbind(cbind(ridx, cidx)[rep(ib, each = 3), 
                                        ], 1:3)
        ibmat[, 1:2] = ibmat[, 1:2] + 1
        return(list(fill = imat, border = ibmat))
      })
      cols.ts = cbind(cols.ts)
      ns = ncol(cols.ts)
      brk.x = sapply(plot.data$width, function(wi) seq(-wi, 
                                                       wi, length = ns + 1))
      for (i in 1:nn) {
        pxr = pidx[[i]]$fill[, 2] - 1 - plot.data$x[i]
        col.rgb = col2rgb(cols.ts[i, ])/255
        for (k in 1:ns) {
          sel = pxr >= brk.x[k, i] & pxr <= brk.x[k + 
                                                    1, i]
          img2[pidx[[i]]$fill[sel, ]] = col.rgb[, k]
        }
        img2[pidx[[i]]$border] = blk
      }
      return(img2)
    }
  }
  else stop("unrecognized node type!")
}

## pathview.stamp, unmodified
pathview.stamp <- function (x = NULL, y = NULL, position = "bottomright", graph.sizes, 
                            on.kegg = TRUE, cex = 1) 
{
  if (on.kegg) 
    labels = "Data on KEGG graph\nRendered by Pathview"
  else labels = "-Data with KEGG pathway-\n-Rendered  by  Pathview-"
  if (is.null(x) | is.null(y)) {
    x = graph.sizes[1] * 0.8
    y = graph.sizes[2]/40
    if (length(grep("left", position)) == 1) 
      x = graph.sizes[1]/40
    if (length(grep("top", position)) == 1) 
      y = graph.sizes[2] - y
  }
  text(x = x, y = y, labels = labels, adj = 0, cex = cex, font = 2)
}
