#' implements a 1 dimensional spline representation of a function.
#' it inherits from FRep which you should refer to to get more
#' details on how to use those functions
#'
#' @export
#' @family frep
F_Spline1D <-  function(support,cc,order=4) {

  fit   = smooth.spline(support,support)
  knots = fit$fit$knot * fit$fit$range + fit$fit$min

  ff <- function(xin,gin) {

    # first we get the levels and compute the design matrices
    vars = list()
    if (class(xin)=='FDiff') {
      xout = xin@F
      vars = xin@vars
    } else {
      xout = xin
    }

    # checking if functional representation is a parameter
    # ----------------------------------------------------
    if (class(gin)=='FDiff') {
      gval = gin@F
      if (length(intersect(names(vars), names(gin@vars)))>0) 
        stop('parmeters and collocation have same variables -- not yet implemented');
      vars = c(vars,gin@vars)
    } else {
      gval = gin
    }
    D  = splineDesign(knots,xout,order,outer.ok=TRUE)
    D1 = splineDesign(knots,xout,order,outer.ok=TRUE,deriv=rep(1,length(xout))) 

    # compute the vector of levels
    # ----------------------------
    F = c(D %*% gval)
    J=NULL

    # checking if x is a parameter (type fdiff)
    # -----------------------------------------
    if (class(xin)=='FDiff') {
      xin@J = Matrix(diag(F),sparse=TRUE) %*% xin@J
      J = expandJacDomain(xin,vars)@J
    }

    # checking if g is of type FDiff
    if (class(gin)=='FDiff') {
      gin@J = Matrix(D,sparse=TRUE)
      J2 = expandJacDomain(gin,vars)@J
      if (is.null(J)) J=J2 else J=J+J2;
    }

    # check if we return a FDiff at all
    # ---------------------------------
    if (is.null(J)) {
      return(F)
    } else {
      return(new("FDiff",F=c(F),J=J,vars=vars))
    }
  }

  class(ff) = 'frep'
  attr(ff,'ng') = length(fit$fit$coef)

  return(ff)
}
