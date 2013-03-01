# scratch file


# This is an example of a carthesian spline function
# It's defined by:
# - one continous variable
# - a discrete index
#

# if I give N support points

SplineInt1D <- function(xsupp,ivals) {
  xknots = splineKnots(xsupp)$knots
  Nx     = splineKnots(xsupp)$N

 return(function(ain,gin) {

  # take care of the function representation
  # for sure gin is a fdiff, otherwise, makes not sense
  if (class(gin) == "FDiff") {
    M = matrix(0,nrow=nrow(cc) , ncol = length(ivals) * Nx)
    F = array(0,nrow(cc))
    
    for (i in ivals) { # quite inneficient!
      I = which(cc$z==i)
      J = ((i-1)*Nx+1) : ((i-1)*Nx+Nx)
      D = splineDesign(xknots,ain[I])
      M[I,J] = D
      F[I]   = F[I] + c(D %*% gin@F[J])
    }
    vars = list(v1 = length(ivals) * (length(xknots)-4))
    names(vars) <- names(ain@vars[[1]])

    R = new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars)

  } else {
    stop('function representation is not a parameter, this seems odd')
  }

  # check if we have an exogneous or endogenous variable
  if (class(ain)=="FDiff") {
    M = matrix(0,nrow=nrow(cc) , ncol = length(ain@F))
    
    for (i in ivals) { # quite inneficient!
      I = which(cc$z==i)
      J = 1:Nx # we need all functional parameters
      D = splineDesign(xknots,ain[I],deriv=rep(1,length(ain[I])))
      M[I,I] = diag(c(D %*% gin[J]))
    }

    R = appendJac(R,Matrix(M,sparse=T),gin@vars)
  }

  return(R)
})

}

test.it <- function() {

  # define the support of the function
  xsupp  = seq(0,1,l=10)
  ivals  = 1:3

  # define a simple collocation
  N=10
  cc = expand.grid(a=xsupp,z=1:3)

  # input arguments
  g. = FDiff(rep(1,length(ivals) * Nx),'g')
  a_ = FDiff(cc$a,'a')

  V = SplineInt1D(xsupp,ivals)
}

# one thing is missing, is how to do t-1

