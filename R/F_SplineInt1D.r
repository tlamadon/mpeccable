require(splines)
require(Matrix)

#' implements a 1 dimensional spline representation of a function.
#' it inherits from FRep which you should refer to to get more
#' details on how to use those functions
#'
#' @export
#' @family frep
#' @example examples/F_SplineInt1D_example.r
F_SplineInt1D <- function(xsupp,ivals) {
  xknots = splineKnotsMpec(xsupp)$knots
  Nx     = splineKnotsMpec(xsupp)$N

 ff <- function(ain,zin,gin,deriv=0) {

  # take care of the function representation
  # for sure gin is a fdiff, otherwise, makes not sense
  if (class(gin) == "FDiff") {
    M = matrix(0,nrow=length(zin) , ncol = length(ivals) * Nx)
    F = array(0,length(zin))
    
    for (i in ivals) { # quite inneficient!
      I = which(zin==i)
      J = ((i-1)*Nx+1) : ((i-1)*Nx+Nx)

      # check for coloring
      if (is.fdiff(ain) & gin@coloring) {
        D = array(1,c(length(I),length(J)))
      } else {
        D = splineDesign(xknots,ain[I],outer.ok = TRUE,derivs = rep(deriv,length(ain[I])))        
      }
      M[I,J] = D
      F[I]   = F[I] + c(D %*% gin@F[J])
    }
    vars = list(v1 = length(ivals) * (length(xknots)-4))
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
      D = splineDesign(xknots,ain[I],deriv=rep(1+deriv,length(ain[I])),outer.ok = TRUE)
      if (ain@coloring) {
        M[I,I] = diag(length(I))  
      } else {
        M[I,I] = diag(c(D %*% gin[J]))   
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