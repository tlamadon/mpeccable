require(splines)
require(Matrix)

#' Spline function with 1 continuous dimension and 1 discrete dimension
#'
#' Creates a functional representation for a 1 dimensional splines indexed by a discrete variable.  
#' This is typically used in the case where you have a function of some continuous variable but also 
#' indexed by time, or by some other discrete variable like a location.
#' This can be combined with control variables and with shape restrictions like monotonicity of 
#' or convexity.
#'
#' @param xsupp a vector of support point for the the continuous dimension
#' @param ivals a vector of discrete values for the support of the discrete dimension
#' @return frep an object of class frep that can be sued to evaluate functions and get levels and Jacobians
#' @export
#' @family frep
#' @example examples/example-MultiSplineFitting.r
F_SplineInt1D <- function(xsupp,ivals) {

  # get the slpines knots and parameter length
  xknots = knot.select(3,xsupp)
  Nx     = length(xsupp)

  ff <- function(ain,zin,gin,deriv=0) {

  # take care of the function representation
  # for sure gin is a fdiff, otherwise, makes not sense
  if (class(gin) == "FDiff") {
    M = matrix(0,nrow=length(zin) , ncol = length(ivals) * Nx,sparse=TRUE)
    F = array(0,length(zin))
    
    for (i in ivals) { # quite inneficient!
      I = which(zin==i)
      J = ((i-1)*Nx+1) : ((i-1)*Nx+Nx)

      # check for coloring
      if (is.fdiff(ain) & gin@coloring) {
        D = array(1,c(length(I),length(J)))
      } else {
        D = splineDesign(xknots,ain[I],outer.ok = TRUE,derivs = rep(deriv,length(ain[I])),ord=4,sparse=TRUE)        
      }
      M[I,J] = D
      F[I]   = F[I] + c(D %*% gin@F[J])
    }
    vars = list(v1 = length(ivals) * Nx)
    names(vars) <- names(gin@vars[[1]])

    if (gin@coloring) M = (M!=0)*1;
    R = new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars,coloring=gin@coloring)

  } else {
    stop('function representation is not a parameter, this seems odd')
  }

  # check if we have an exogneous or endogenous variable
  if (class(ain)=="FDiff") {
    M = matrix(0,nrow=length(xsupp)*length(ivals) , ncol = length(ain@F))
    
    for (i in ivals) { # quite inneficient!
      I  = which(zin==i)
      J  = 1:Nx # we need all functional parameters
      D = splineDesign(xknots,ain[I],deriv=rep(1+deriv,length(ain[I])),outer.ok = TRUE,ord=4,sparse=TRUE)
      if (ain@coloring) {
        M[I,I] = Diagonal(length(I))  
      } else {
        M[I,I] = Diagonal(length(I),c(D %*% gin[J]))   
      }
    }
    R = appendJac(R,Matrix(M,sparse=T),ain@vars)
  }
  return(R)
}

class(ff) = 'frep'
attr(ff,'ng') = length(ivals) * Nx
return(ff)

}