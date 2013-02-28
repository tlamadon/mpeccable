# initializes a 1D spline represenation function
#' @export
#' @family frep
F_Spline2D <-  function(support,cc,order=4)  {
  # you can control the knots in each dimensions
  # but you get the product of each dimension so maybe use (5,5,30)

  # preparing the knots for the spline
  fit   = smooth.spline(support,support)
  rknots = fit$fit$knot * fit$fit$range + fit$fit$min

  # prepcompute the design matrices for x and y
  # cause that does not change
  Dxy = cbind( 
          rep(1,length(cc$x)),
          cc$x,
          cc$x^2)
  # then I create a sparse Matrix that will 
  # multiply the spline design
  Dr   = Matrix(splineDesign(rknots,cc$v,order,outer.ok=TRUE),sparse=TRUE)
  nc = length(rknots) - 4
  D = Dxy[,1]%*% t(rep(1,nc))
  for (n in 2:ncol(Dxy)) {
    D = cBind(D, Dr * Dxy[,n] %*% t(rep(1,nc)))    
  }

  # creating the function that evaluates
  # ------------------------------------
  ff <- function(rin,gin,xin=NULL,deriv=0) {

    # preparation 
    # -----------

    # dealing with x and y which are exogenously
    # given but potentially at a different collocation
    if (is.null(xin)) {
      xin = cc$x
    } else {
       Dxy = cbind( 
          rep(1,length(xin)),
          xin,
          xin^2)
    }

    # processing the rho argument which might be a control
    vars = list()
    if (class(rin)=='FDiff') {
      rout = rin@F
      vars = rin@vars
    } else {
      rout = rin
    }

    # processing the represenation parameters
    if (class(gin)=='FDiff') {
      gval = gin@F
      if (length(intersect(names(vars), names(gin@vars)))>0) 
        stop('parmeters and collocation have common variables -- not yet implemented');
      vars = c(vars,gin@vars)
    } else {
      gval = gin
    }
    Dr1  = Matrix(splineDesign(rknots, rout,order, outer.ok=TRUE, deriv=rep(deriv+1,length(rout))),sparse=TRUE)
    Dr   = Matrix(splineDesign(rknots, rout,order, outer.ok=TRUE, deriv=rep(deriv,  length(rout))),sparse=TRUE)

    # overall design matrix
    # we need to cbind replication of Dr for each column of xy
    nc   = length(rknots) - 4
    D    =  Dr  * Dxy[,1]%*% t(rep(1,nc))
    D1   =  Dr1 * Dxy[,1]%*% t(rep(1,nc))
    for (n in 2:ncol(Dxy)) {
      D  = cBind(D,  Dr  * Dxy[,n] %*% t(rep(1,nc)))    
      D1 = cBind(D1, Dr1 * Dxy[,n] %*% t(rep(1,nc)))    
    }

    # compute the value
    # -----------------
    F = as.numeric(D %*% gval)

    # compute the jacobian
    # --------------------
    J=NULL
    # checking if x is of type FDiff
    if (class(rin)=='FDiff') {
      rin@J = Matrix(diag(F),sparse=TRUE) %*% rin@J
      J = expandJacDomain(rin,vars)@J
    }

    # checking if g is of type FDiff
    if (class(gin)=='FDiff') {
      gin@F = F
      gin@J = Matrix(D,sparse=TRUE)
      J2 = expandJacDomain(gin,vars)@J
      if (is.null(J)) J=J2 else J=J+J2;
    }

    # check if we return a FDiff
    if (is.null(J)) {
      return(F)
    } else {
      return(new("FDiff",F=c(F),J=J,vars=vars))
    }
  }

  class(ff) = 'frep'
  attr(ff,'ng') = 3 * length(fit$fit$coef)

  return(ff)
}
