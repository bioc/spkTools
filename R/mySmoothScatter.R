mySmoothScatter <- function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", 
    brewer.pal(9, "Blues"))), nrpoints = 100, transformation = function(x) x^0.25, 
    xlab = NULL, ylab = NULL, postPlotHook = box, pch = ".", 
    cex = 1, xlim, ylim, col = "black", xaxs = par("xaxs"), yaxs = par("yaxs"), 
    ...) 
{
    if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) != 
        1)) 
        stop("'nrpoints' should be numeric scalar with value >= 0.")
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel)
    xlab <- if (is.null(xlab)) 
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab)) 
        xy$ylab
    else ylab
    x <- cbind(xy$x, xy$y)[!(is.na(xy$x) | is.na(xy$y)), ]
    if (!missing(xlim) & !is.null(xlim)) {
        stopifnot(is.numeric(xlim), length(xlim) == 2, !any(is.na(xlim)))
        x <- x[(x[, 1] >= xlim[1]) & (x[, 1] <= xlim[2]), ]
    }
    else {
        xlim <- range(x[, 1], na.rm = TRUE)
    }
    if (!missing(ylim) & !is.null(ylim)) {
        stopifnot(is.numeric(ylim), length(ylim) == 2, !any(is.na(ylim)))
        x <- x[(x[, 2] >= ylim[1]) & (x[, 2] <= ylim[2]), ]
    }
    else {
        ylim <- range(x[, 2], na.rm = TRUE)
    }
    map <- .smoothScatterCalcDensity(x, nbin, bandwidth)
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    dens <- array(transformation(dens), dim = dim(dens))
    image(xm, ym, z = dens, col = colramp(256), xlab = xlab, 
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs, 
        ...)
    if (!is.null(postPlotHook)) 
        postPlotHook()
    if (nrpoints != 0) {
        stopifnot(length(xm) == nrow(dens), length(ym) == ncol(dens))
        ixm <- round((x[, 1] - xm[1])/(xm[length(xm)] - xm[1]) * 
            (length(xm) - 1))
        iym <- round((x[, 2] - ym[1])/(ym[length(ym)] - ym[1]) * 
            (length(ym) - 1))
        idens <- dens[1 + iym * length(xm) + ixm]
        nrpoints <- min(nrow(x), ceiling(nrpoints))
        sel <- order(idens, decreasing = FALSE)[1:nrpoints]
        points(x[sel, 1:2], pch = pch, cex = cex, col = col)
    }
}
