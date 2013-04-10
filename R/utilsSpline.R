

#' Spline Knot Selector
#'
#' @description takes a grid of points and places spline knots such that a) # of basis funs = # of data sites and b) the placement of knots 
knot.select <- function(degree,grid,plotit=FALSE){
# returns a knotvector for a grid of data sites and a spline degree such that number of data sites = number of basis functions
    n <- length(grid)
    p <- n+degree+1     # number of nodes required for solving Ax=b exactly
    knots <- rep(0,p)
    knots <- replace(knots,1:(degree+1),rep(grid[1],degree+1))
    knots <- replace(knots,(p-degree):p,rep(tail(grid,1),degree+1))
# this puts multiplicity of first and last knot in order
# if there are anough gridpoints, compute intermediate ones
    if (n<(degree+1)) stop("to few grid points for clamped curve")
    if (n>(degree+1)){
        for (j in 2:(n-degree)) knots[j+degree] <- mean(grid[j:(j+degree-1)])
    }
	if (plotit) {
	  plot(x=grid,y=rep(1,n))
	  points(x=knots,y=rep(1,p),pch=3)
	}
    return(knots)
}



#' get spline representation
splineKnotsMpec <- function(support) {
  fit   = smooth.spline(support,support)
  knots = fit$fit$knot * fit$fit$range + fit$fit$min
  return(list(knots = knots, N = fit$fit$nk))
}


#' select knots for splineDesign at quantiles of data
#'
#' user selects number of desired basis functions and degree of spline, function returns knot vector with knot multiplicity \emph{degree} at both ends. 
#' @param degree positive integer for spline degree
#' @param x numeric vector of data sites at which spline will be evaluated
#' @param num.basis integer of desired number of basis functions
#' @param plotit logical of whether plot result
#' @param stretch numeric value indicating if and by which percentage you want to stretch the knots over the data support. Useful when approximating values at or near bounds of data.
#' @detail the crucial relationship is num.basis = \code{length(knots)} - degree - 1, which we use to find \code{length(knots)} = num.basis + degree + 1. The minimum number of basis functions to obtain a valid knot vector with correct multiplicity is \code{min(num.basis) = deg + 1}. 
#' @return numeric vector of spline knots of class \emph{knotVec} with attribute \emph{num.basis}
#' @examples 
#' knot.select2(degree=3,x=1:10,num.basis=6,stretch=0.01,plotit=TRUE)
knot.select2 <- function(degree,x,num.basis,stretch=NULL,plotit=FALSE){
    n <- length(x)
	x <- sort(x)
	if (!is.null(stretch)){
		r <- range(x)
		low <- r[1] - stretch * diff(r)
		hi  <- r[2] + stretch * diff(r)
	} else {
		low <- x[1]
		hi <- x[n]
	}
	kdown <- rep(low,times=degree+1)
	kup   <- rep(hi,times=degree+1)
	if (num.basis < degree + 1) {
		num.basis <- degree + 1
		warning(c("you chose too few basis functions. I selected the required minimum of ",num.basis))
	}
	# any remaining points for interior knots?
	z <- num.basis - (degree+1)
	iknots <- NULL
	if (z > 0) {
		quants <- seq.int(from=0,to=1,length.out=z+2)[-c(1,z+2)]	# quantiles of data sites do consider (excluding first and last, since those are included in kdown and kup)
		iknots <- quantile(x,quants)
	}
	knots <- c(kdown, iknots, kup)
	names(knots) <- NULL
	stopifnot(length(knots) == num.basis + degree + 1)
	if (plotit){
		plot(x=x,y=rep(1,n),yaxt='n',ylab="",ylim=c(0.8,1.4),xlim=range(knots),xlab="data index",main=sprintf("Spline Knots and Data\nsetup implies %s basis functions",num.basis),sub=sprintf("note knot multiplicity of %s at first and last data point",degree+1))
		points(x=knots,y=rep(1.2,length(knots)),pch=3)
		legend("bottomright",legend=c("data","knots"),pch=c(1,3))
	}
	class(knots) <- 'knotVec'
	attr(knots,'num.basis') <- num.basis
	return(knots)
}





  




























function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
    Boundary.knots = range(x)) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1) 
        stop("'degree' must be integer >= 1")
    if (!missing(df) && missing(knots)) {
        nIknots <- df - ord + (1 - intercept)
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used  ", ord - 
                (1 - intercept))
        }
        knots <- if (nIknots > 0) {
            knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                2)[-c(1, nIknots + 2)]
            stats::quantile(x[!outside], knots)
        }
    }



