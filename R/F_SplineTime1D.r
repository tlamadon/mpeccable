
require(splines)
require(Matrix)

#' Spline function with 1 continuous dimension and 1 discrete dimension defined on finite Time
#'
#' Creates a functional representation for a 1 dimensional splines indexed by a discrete variable and time.
#' time is discrete and assumed to end at period T. This functional is useful for situations where the 
#' form of the structural function V is known in T, i.e. does not have to be approximated.
#'
#' @param xsupp a vector of support point for the the continuous dimension
#' @param ivals a vector of discrete values for the support of the discrete dimension
#' @param degree integer of desired spline degree
#' @param nbasis.fun integer for desired number of basis functions
#' @return frep an object of class frep that can be sued to evaluate functions and get levels and Jacobians
#' @export
#' @seealso \code{\link{F_SplineInt1D.r}} 
#' @family frep
#' @example examples/example-MultiSplineFitting.r
F_SplineTime1D <- function(xsupp,ivals,degree,nbasis.funs) {

  # get the slpines knots and parameter length
  xknots <- knot.select2(degree=degree,x=xsupp,num.basis=nbasis.funs,stretch=0.01)
  Nx     <- length(xsupp)            # number of data points in each discrete bin
  Nm     <- attr(xknots,'num.basis') # number of spline coefficients in each discrete bin

  ff <- function(ain,zin,gin,deriv=0) {


	if (class(ain)=="FDiff"){
	 if (ain@coloring) browser()
	}
	# spline parameter gin
	# --------------------------

	  # compute value of approximation to V(ain,gin) and return in F
	  # compute partial derivative of approximation to V(ain,gin) w.r.t gin and return in M
	
	  if (class(gin) == "FDiff") {

		M = Matrix(0,nrow=length(zin) , ncol = length(ivals) * Nm, sparse=TRUE)	# M not necessarily square
		F = array(0,length(zin))
		
		# here it would be good to be more flexbile
		# with timing, there is always a last period
		# e.g. i don't want to know the euler equation in the last period, it's not defined.
		#     key <- data.table(expand.grid(ia=1:Na,iz=1:2,it=1:5),key="it")
		# key[,index := 1:nrow(key)]
		for (i in unique(zin)){
		  I = which(zin==i)	# row indices of M for discrete bin i
		  J = ((i-1)*Nm+1) : ((i-1)*Nm+Nm)	# col indices
		
		  # check for final period
		  if (i == maxtime){
			F[I] = log( ain[I] )	# suppose final function is log
			#         M[I,J] = diag(length(I))	there is zero impact of the coefficients in the last period, so don't change.
		  } else {
			  # check for coloring
			  if (is.fdiff(ain) & gin@coloring) {
				D = array(1,c(length(I),length(J)))
			  } else {
				# return the value of the basis function evaluated at ain
				D = splineDesign(xknots,ain[I],derivs = rep(deriv,length(ain[I])),outer.ok=TRUE,ord=(degree+1),sparse=TRUE)        
			  }
			  M[I,J] = D
			  F[I]   = F[I] + as.numeric(D %*% gin@F[J])
			}

		  }
		vars = list(v1 = length(ivals) * Nm)
		names(vars) <- names(gin@vars[[1]])

		if (gin@coloring) M = (M!=0)*1;
		R = new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars,coloring=gin@coloring)

	} else {
	  stop('function representation is not a parameter, this seems odd')
	}

	# endogenous choice variable ain
	# ------------------------------

	# check if we have an exogneous or endogenous variable
	if (class(ain)=="FDiff") {
	  M = Matrix(0, nrow=length(zin) , ncol = length(ain@F),sparse=T)	# M is square
	  
	  for (i in unique(zin)){
		I = which(zin==i)
		J  = 1:Nm
		# check for final period
		if (i == maxtime){
		  M[I,I] = Diagonal(length(I), 1/ ain[I] )	# suppose final function is log
		} else {
		  D = splineDesign(xknots,ain[I],deriv=rep(1+deriv,length(ain[I])),outer.ok = TRUE,ord=degree+1,sparse=TRUE)
		  if (ain@coloring) {
			M[I,I] = Diagonal(length(I))  
		  } else {
			M[I,I] = Diagonal(length(I), as.numeric(D %*% gin[J]))   
		  }
		}
	  }
	  R = appendJac(R,Matrix(M,sparse=T),ain@vars)
	}
	return(R)
  }
  class(ff) = 'frep'
  attr(ff,'ng') = length(ivals) * Nm
  return(ff)
}
