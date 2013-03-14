# is setGeneric similar to a "virtual" class in c++?
##########################
# knot selector
##########################

knot.select <- function(degree,grid){
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
    return(knots)
}



#' get spline representation
splineKnotsMpec <- function(support) {
  fit   = smooth.spline(support,support)
  knots = fit$fit$knot * fit$fit$range + fit$fit$min
  return(list(knots = knots, N = fit$fit$nk))
}
