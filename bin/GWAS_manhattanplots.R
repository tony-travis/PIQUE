#' @(#)GWAS_manhattanplots.R  2018-06-21  A.Douglas and A.J.Travis
#'
#' Heavily modified from original code provided by Stephen Turner
#' http://GettingGeneticsDone.blogspot.com/
#'
#' R code for making Manhattan plots and QQ plots from plink format files.
#' must contain column headings CHR, BP, P
#'
#' AD added a couple of arguments to format the plots to allow stacking of separate
#' Manhattan plots (stack). Also added arguments to plot marker info (markers) and
#' QTL blocks(blocks)
#'
#' Marker info must be a dataframe with column headings 'BasePair', 'Chrom', 'colour' and 'line'.
#' BasePair = the basepair position of the marker WITHIN a chromosme
#' Chrom = the chromosome where the marker occurs
#' colour = the colour of the line to highlight the marker on the plot (must be a valid R colour - alphanumeric)
#' line = the style of the marker line (must be a valid R line style)
#'
#' QTL block info must be a dataframe with column headings 'chrom', 'start', 'stop', 'colour'
#' chrom = the chromosome number on which the QTL block occurs
#' start = the basepair position WITHIN a chromosome that denotes the start of the QTL block
#' stop = the basepair position WITHIN a charomosome that denotes the end of the QTL block
#' colour = the colour to fill the QTL block (must be a colour recognised by R - alphanumeric)     
#'
#' This is for testing purposes.
#' set.seed(42)
#' nchr=23
#' nsnps=1000
#' d=data.frame(
#'   SNP=sapply(1:(nchr*nsnps), function(x) paste('rs',x,sep='')),
#'   CHR=rep(1:nchr,each=nsnps),
#'   BP=rep(1:nsnps,nchr),
#'   P=runif(nchr*nsnps)
#' )
#' annotatesnps <- d$SNP[7550:7750]
#'

suppressMessages(library(scales))  # include check if installed
suppressMessages(library(Hmisc))  # include check if installed

#'
#' Manhattan plot using base graphics
#'
#' dataframe
#' stack = FALSE
#' markers = NULL
#' blocks = NULL
#' colors = c('gray10', 'gray50')
#' ymax = 'max'
#' xaxis.cex = 1
#' limitchromosomes = chr.limit
#' suggestiveline = -log10(1e-05)
#' genomewideline = -log10(5e-08)
#' logp = TRUE
#' highlight = NULL
#' annotate = NULL
#' gene.list = NULL
#' ...
#'
manhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
  "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
  highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = TRUE, bh = "BH", 
  stack = FALSE, markers = NULL, blocks = NULL, ymax = "max", xaxis.cex = 1, limitchromosomes = NULL, 
  gene.list = NULL, ...) {
  
  # Check for sensible dataset Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  
  # remove na's
  x <- na.omit(x)
  
  # limits chromosomes to plot. (23=x, 24=y, 25=par?, 26=mito?)
  if (any(as.integer(limitchromosomes))) 
    x <- x[x$CHR %in% limitchromosomes, ]
  
  #'
  #' Create a new data.frame with columns called CHR, BP, and P.
  #' d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA)
  #' with millions of SNPs, create dataframe at once rather than dynamically
  #' allocated(see line 72-73, and remove line 87 and line 91 )
  #'
  d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], pos = NA, index = NA)
  if (!is.null(x[[snp]])) {
    
    # If the input data frame has a SNP column, add it to the new data frame
    d$SNP <- x[[snp]]
  }
  if (!is.null(x[[bh]])) {
    
    # If the input data frame has a BH column, add it to the new data frame
    d$BH <- x[[bh]]
  }
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    
    # use -log10(p) value
    d$logp <- -log10(d$P)
  } else {
    
    # use raw P value
    d$logp <- d$P
  }
  
  # Fixes the bug where one chromosome is missing by adding a sequential index
  # column
  d$index <- rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR, length))
  
  # Create a vector of alternating colors for chromosomes
  col <- rep_len(col, max(d$index))
  
  # set maximum value on y axis (-log10p scale)
  if (ymax == "max") {
    ymax <- ceiling(max(d$logp))
    if (ymax < 8) 
      ymax <- 8
  } else ymax <- as.numeric(ymax)
  
  # gap between chromosomes
  gap <- ceiling(max(d$BP)) * 0.04
  
  # This section sets up positions and ticks. Ticks should be placed in the middle
  # of a chromosome. The a new pos column is added that keeps a running sum of the
  # positions of each successive chromsome. For example: chr bp pos 1 1 1 1 2 2 2 1
  # 3 2 2 4 3 1 5
  nchr <- length(unique(d$CHR))
  if (nchr == 1) {
    lastbase <- tail(d[d$CHR == limitchromosomes, ]$BP, 1)/1e+06
    d$pos <- d$BP/1e+06
    xlabel <- paste("Chromosome", unique(d$CHR), "position (Mb)")
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      lastbase <- lastbase + head(d[d$CHR == i, ]$BP, 1)
      d[d$CHR == i, ]$pos <- d[d$CHR == i, ]$BP + lastbase
      nextbase <- lastbase + tail(d[d$CHR == i, ]$BP, 1)
      ticks <- c(ticks, (lastbase + nextbase)/2)
      lastbase <- nextbase + gap
    }
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }
  
  #' Initialize plot
  #' xmin <- floor(lastbase * -0.03)
  #' xmax <- ceiling(lastbase * 1.03)
  xmin <- 0
  xmax <- lastbase
  
  # default options
  xlim <- c(xmin, xmax)
  ylim <- c(0, ymax)
  xlab <- xlabel
  ylab <- expression(-log[10](italic(p)))
  #' xaxt = 'n'
  #' bty = 'n'
  #' xaxs = 'i'
  #' yaxs = 'i'
  #' las = 1
  #' pch = 20
  #'
  
  # create plot
  par(mgp = c(2, 0.5, 0), cex = 1, xaxs = "i", yaxs = "i", pch = 20, pty = "m", 
    bty = "n")
  if (nchr == 1) {
    
    # single chromosome
    if (stack == TRUE) {
      
      # suppress x-axis when stacking manahttan plots
      with(d, plot(pos, logp, cex = 0.5, xlim = xlim, ylim = ylim, ylab = "", 
        xaxt = "n", yaxt = "n", ...))
    } else {
      
      # x-axis: chromosomal position, y-axis: -log10(p)
      with(d, plot(pos, logp, cex = 0.5, xlim = xlim, ylim = ylim, ylab = ylab, 
        xlab = xlab, yaxt = "n", ...))
      minor.tick(nx = 5, ny = 1)
    }
  } else {
    
    # multiple chromosomes: first plot with no x-axis (xaxt='n')
    if (stack) {
      
      # suppress x-axis when stacking manahttan plots
      with(d, plot(pos, logp, cex = 0.5, xlim = xlim, ylim = ylim, ylab = "", 
        xaxt = "n", yaxt = "n", type = "n", ...))
    } else {
      
      # x-axis: chromosome number instead of position
      with(d, plot(pos, logp, cex = 0.5, xlim = xlim, ylim = ylim, ylab = "", 
        xlab = xlabel, xaxt = "n", yaxt = "n", ...))
      axis(1, at = ticks, labels = labs)
    }
  }
  
  # y-axis with custom tick marks
  major <- seq(0, ymax, 4)
  axis(2, at = major)
  minor <- seq(0, ymax, 1)
  axis(2, at = minor, tck = -0.04, labels = rep("", length(minor)))
  
  # SNPs
  icol <- 1
  for (i in unique(d$CHR)) {
    with(d[d$CHR == i, ], points(pos, logp, cex = 0.5, col = col[icol], ...))
    icol <- icol + 1
  }
  
  # marker lines
  if (!is.null(markers)) {
    markers <- markers
    markers$marker.pos <- NA
    base <- 0
    nchr <- length(unique(markers$Chrom))
    if (nchr == 1) {
      markers$marker.pos <- markers$BasePair
    } else {
      for (i in markers$Chrom) {
        if (i == 1) {
          markers[markers$Chrom == i, ]$marker.pos <- markers[markers$Chrom == 
          i, ]$BasePair
        } else {
          base <- tail(subset(d, CHR == i - 1)$pos, 1)
          markers[markers$Chrom == i, ]$marker.pos <- markers[markers$Chrom == 
          i, ]$BasePair + base
        }
      }
    }
    for (m in 1:length(markers$BasePair)) {
      abline(v = markers$marker.pos[m], lwd = 1.2, col = rgb(red = t(col2rgb(markers$colour[m])), 
        alpha = 90, maxColorValue = 255), lty = markers$line[m])
    }
  }
  
  # highlight blocks
  if (!is.null(blocks)) {
    qtl.block <- blocks
    
    # put block code here
    qtl.block$qtl.start.pos <- NA
    qtl.block$qtl.stop.pos <- NA
    base <- 0
    nchr <- length(unique(qtl.block$chrom))
    if (nchr == 1) {
      qtl.block$qtl.start.pos <- qtl.block$start
      qtl.block$qtl.stop.pos <- qtl.block$stop
    } else {
      for (i in qtl.block$chrom) {
        if (i == 1) {
          qtl.block[qtl.block$chrom == i, ]$qtl.start.pos <- qtl.block[qtl.block$chrom == 
          i, ]$start
          qtl.block[qtl.block$chrom == i, ]$qtl.stop.pos <- qtl.block[qtl.block$chrom == 
          i, ]$stop
          
        } else {
          base <- tail(subset(d, CHR == i - 1)$pos, 1)
          qtl.block[qtl.block$chrom == i, ]$qtl.start.pos <- qtl.block[qtl.block$chrom == 
          i, ]$start + base
          qtl.block[qtl.block$chrom == i, ]$qtl.stop.pos <- qtl.block[qtl.block$chrom == 
          i, ]$stop + base
        }
      }
    }
    out.blocks <<- qtl.block  #test
    for (b in 1:length(qtl.block$chrom)) {
      rect(qtl.block$qtl.start.pos[b], -0.5, qtl.block$qtl.stop.pos[b], ymax + 
        0.5, col = rgb(red = t(col2rgb(qtl.block$colour[b])), alpha = 80, 
        maxColorValue = 255), border = NA)
    }
  }
  
  # create a new data frame with rows from the original data frame where SNP is in
  # annotate character vector.  then plot those points over the original graph, but
  # with a larger point size and a different color.
  if (!is.null(annotatePval)) {
    d.annotatePval <- d[which(d$SNP %in% annotate.Pval), ]
    with(d.annotatePval, points(pos, logp, cex = 0.5, col = "red", cex = 2.5, 
      ...))
  }
  
  # Highlight significant SNPs
  if (!is.null(highlight)) {
    
    # extract significant SNPs
    topHits <- subset(d, BH <= highlight)
    with(topHits, points(pos, -log10(P), cex = 0.5, col = "red"), ...)
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue", lty = 2)
  if (genomewideline) 
    abline(h = genomewideline, col = "red", lty = 2)
  
  # add line in position of gene
  if (!is.null(gene.list)) {
    gene.list <- gene.list[gene.list$seqid %in% limitchromosomes, ]
    with(gene.list, {
      abline(v = start/1000, col = "red", lty = 2)
      text(start/1000, ymax, attributes, adj = 1, srt = 90)
    })
  }
}

# Base graphics qq plot
qq <- function(d, ymax = "max", limitchromosomes = 1:23, highlight = 0.05, showBH = FALSE, 
  xaxis.cex = 1, stack = FALSE, ...) {
  
  # sort by observed probabilities
  d <- d[order(d$P), ]
  
  # limit chromosomes to plot. (23=x, 24=y, 25=par?, 26=mito?)
  if (any(as.integer(limitchromosomes))) 
    d <- d[d$CHR %in% limitchromosomes, ]
  
  # remove na's
  d <- na.omit(d)
  
  # raw observed probabilities
  o1 <- -log10(d$P)
  
  # Benjamini-Hochberg corrected observed probabilities
  o2 <- -log10(d$BH)
  
  # expected probabilities
  e1 <- -log10(ppoints(length(o1)))
  
  # maximum value on y axis (-log10p scale).
  if (ymax == "max") {
    ymax <- ceiling(max(o1))
    if (ymax < 8) 
      ymax <- 8
  } else ymax <- as.numeric(ymax)
  
  # create plot
  par(mgp = c(2, 0.5, 0), cex = 1, xaxs = "i", yaxs = "i", pch = 20, pty = "s", 
    bty = "l")
  if (stack) {
    
    # suppress x-axis when stacking QQ plots
    plot(e1, o1, cex = 0.5, xlim = c(0, ymax), ylim = c(0, ymax), xlab = "", 
      ylab = "", asp = 1, xaxt = "n", yaxt = "n", ...)
  } else {
    
    # x-axis: Expected -log10(P) values, y-axis Observed -log10(P) values
    plot(e1, o1, cex = 0.5, xlim = c(0, ymax), ylim = c(0, ymax), xlab = expression(Expected - 
      log[10](italic(p))), ylab = expression(Observed - log[10](italic(p))), 
      asp = 1, xaxt = "n", yaxt = "n", ...)
    
    # x-axis annotation
    minor.tick(nx = 2, ny = 1)
    axis(1, at = seq(0, ymax, 4))
  }
  
  # y-axis annotation
  minor.tick(nx = 1, ny = 2)
  axis(2, at = seq(0, ymax, 4))
  
  # highlight significant SNPs
  top <- data.frame(e1, o1, o2)
  top <- subset(top, o2 >= -log10(highlight))
  points(top$e, top$o1, cex = 0.5, col = "red")
  
  # Benjamini-Hochberg corrected observed probabilities
  if (showBH) {
    e2 <- e1
    points(e2, o2, cex = 0.5, col = "blue")
    abline(h = -log10(highlight), col = "blue", lty = 2)
  }
  
  # 1:1 line
  abline(0, 1, col = "blue", lty = 2)
}
