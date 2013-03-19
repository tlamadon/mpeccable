require(splines)
require(Matrix)

#' One dimensional spline
#'
#' @param xsupp grid to be used for the spline
#' @export
#' @family frep
#' @example examples/example-SplineFittingShapeRestriction.r
F_Spline <- function(xsupp) {

  #xknots = splineKnotsMpec(xsupp)$knots
  #Nx     = splineKnotsMpec(xsupp)$N

  xknots = knot.select(3,xsupp)
  Nx     = length(xsupp)

  ff <- function(xin,gin,deriv=0) {

    # take care of the function representation
    # for sure gin is a fdiff, otherwise, makes not sense
    if (class(gin) == "FDiff") {      
      M = splineDesign(xknots,xin,derivs = rep(deriv,length(xin)),ord=4,sparse=TRUE)        
      F = M %*% gin@F
      if (gin@coloring) M = (M!=0)*1;

      vars = list(v1 = Nx)
      names(vars) <- names(gin@vars[[1]])
      R = new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars,coloring=gin@coloring)

    } else {
      stop('function representation is not a parameter, this seems odd')
    }

  return(R)
  }

  class(ff) = 'frep'
  attr(ff,'ng') = Nx
  return(ff)
}