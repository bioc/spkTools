spkViolinplot <- function(x, ylim, orientation = "vertical", bw = "nrd0",
                           names = NULL, pars = NULL, i1=NULL){
    N <- length(i1)
    groups <- if (is.list(x))
        x
    if (0 == (n <- length(groups)))
        stop("invalid first argument")
    if (length(class(groups)))
        groups <- unclass(groups)
    if (!missing(names))
        attr(groups, "names") <- names[i1]
    else {
        if (is.null(attr(groups, "names")))
            attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    xvals <- matrix(0, nr = 512, nc = n)
    yvals <- matrix(0, nr = 512, nc = n)
    center <- which(i1)
    for (i in 1:n) {
        tmp.dens <- density(groups[[i]], bw = bw)
        xvals[, i] <- tmp.dens$x
        yvals.needtoscale <- tmp.dens$y
        yvals.scaled <- 7/16 * yvals.needtoscale/max(yvals.needtoscale)
        yvals[, i] <- yvals.scaled
    }
    if (orientation == "vertical") {
        xrange <- c(1/2, N + 1/2)
        yrange <- range(xvals)
    }
    else {
        xrange <- range(xvals)
        yrange <- c(min(yvals), max(yvals))
    }
    plot.new()
    if(is.null(ylim)) plot.window(xlim = xrange, ylim = yrange)
    if(!is.null(ylim)) plot.window(xlim = xrange, ylim = ylim)
    for (i in 1:n) vlnplt(xvals[, i], yvals[, i], center[i],
        bordercolor = rainbow(i), bgcolor = rainbow(n - i), orientation = orientation)
    axis(1, at = 1:N, labels = names)
    axis(2)
}
