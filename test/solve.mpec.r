

# scratch file 2


SplineInt1D <- function(xsupp,ivals) {
  xknots = splineKnots(xsupp)$knots
  Nx     = splineKnots(xsupp)$N

 return(function(ain,zin,gin) {

  # take care of the function representation
  # for sure gin is a fdiff, otherwise, makes not sense
  if (class(gin) == "FDiff") {
    M = matrix(0,nrow=nrow(cc) , ncol = length(ivals) * Nx)
    F = array(0,nrow(cc))
    
    for (i in ivals) { # quite inneficient!
      I = which(zin==i)
      J = ((i-1)*Nx+1) : ((i-1)*Nx+Nx)

      # check for coloring
      if (is.fdiff(ain) & getOption('mpeccable.coloring')) {
        D = array(1,c(length(I),length(J)))
      } else {
        D = splineDesign(xknots,ain[I])        
      }
      M[I,J] = D
      F[I]   = F[I] + c(D %*% gin@F[J])
    }
    vars = list(v1 = length(ivals) * (length(xknots)-4))
    names(vars) <- names(gin@vars[[1]])

    R = new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars)

  } else {
    stop('function representation is not a parameter, this seems odd')
  }

  # check if we have an exogneous or endogenous variable
  if (class(ain)=="FDiff") {
    M = matrix(0,nrow=nrow(cc) , ncol = length(ain@F))
    
    for (i in ivals) { # quite inneficient!
      I  = which(zin==i)
      J  = 1:Nx # we need all functional parameters
      D = splineDesign(xknots,ain[I],deriv=rep(1,length(ain[I])))
      if (getOption('mpeccable.coloring')) {
        M[I,I] = diag(length(I))  
      } else {
        M[I,I] = diag(c(D %*% gin[J]))   
      }
    }

    R = appendJac(R,Matrix(M,sparse=T),ain@vars)
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

  options(mpeccable.coloring=TRUE)
  R = V(cc$a,cc$z,g.)
  R = V(a_,  cc$z,g.)

  #options(mpeccable.coloring=FALSE)
  # create a low of motion matrix, where 1->2 , 2->3, 3->1
  z1 = (cc$z %% 3)+ 1 

  # create a simple Euler Equation
  R2 = a_ + V(cc$a,cc$z,g.) + 0.9 * V(a_,z1,g.)
  image(R2@J)

}

# one thing is missing, is how to do t-1





# develop a class SolveMpec that wraps the whole thing and sends it off to a solver

setClass(Class="SolveMpec", representation(Func = "FDiff", x0="numeric", ub="numeric", lb="numeric",sol="list"))

setMethod("print","SolveMpec",
		  function(x,...){
			  cat("*** Class SolveMpec, print method *** \n")
			  cat("* initial value = "); print(x@x0)
			  cat("* ub = "); print(x@ub)
			  cat("* lb = "); print(x@lb)
			  cat("* Func = "); print(x@Func)
			  cat("* Solution = "); print(x@sol)
			  cat("*** End print (SolveMpec) *** \n")
		  }
		  )

setMethod(f="show",signature="SolveMpec",
		  definition=function(object){
			  cat("*** Class SolveMpec, show method *** \n")
			  cat("\n")
			  cat("* optimization parameters\n")
			  cat("* -----------------------\n")
			  cat("* initial value = "); print(head(object@x0,20))
			  cat("* ub = "); print(head(object@ub,20))
			  cat("* lb = "); print(head(object@lb,20))
			  cat("\n")
			  cat("* Elements of Functional Representation:\n")
			  cat("* --------------------------------------\n")
			  cat("* function Values F = \n"); print(head(object@Func@F,20))
			  cat("* choice variables  = \n"); 
			  print(object@Func@vars)
			  if (length(object@Func@J)>0){
				  cat("* Jacobian[max 30 by 30] = "); 
				  rows <- nrow(object@Func@J)
				  cols <- ncol(object@Func@J)
				  print(object@Func@J[1:min(min(rows,cols),30),1:min(min(rows,cols),30)])
			  } else {
				  cat("* empty jacobian slot\n")
			  }
			  cat("*** End show (SolveMpec) *** \n")
		  }
		  )


# try print and show
M <- new("SolveMpec",Func=R,x0=rnorm(100),ub=rep(3,10),lb=runif(19))
M

# initialize an empty object and show it
M2 <- new("SolveMpec")
M2



# the core function: solve()
setMethod(f="solve",signature="SolveMpec",
		  definition=function(obj){
			  ...
			  res <- ipoptr(x0=obj@x0,
							lb=obj@lb,
							ub=obj@ub,
							eval_f=objective,
							eval_grad_f=grad,
							eval_g=constfun,
							eval_jac_g=jacfun,
							eval_jac_g=constfun,
							eval_jac_g_structure=jac.struc)
			  # you get jac.struc outside of here by
			  # jac.struc <- make.sparse( R@J )		and this will remain fixed throughout the process.

			  # process result

			  # store results in respective slots

			  # return solver status code
		  }
		  )




# objective
objective <- function(x) {
	# find a way to extract errors from x. 
	# right now there is a data.table with a key "variable":
	return( sum(  x[ h[c("err.1","err.2")][,idx] ]^2   ))
}

# gradient of objective
grad <- function(x) {
	return( 2 * x[ h[c("err.1","err.2")][,idx] ] )
}

constfun <- function(x,FDiffobj){
# input arguments
	
	# which elemnents of x are errors, a_ and g. respectively?
	newg <- x[ "gamma" ]
	newa <- x[ "endog" ]
	errors <- x[ "errors" ]

  g. = FDiff(newg,'g')
  a_ = FDiff(newa,'a')

  V = SplineInt1D(xsupp,ivals)	# this can be stored on disk/memory

  R = V(cc$a,cc$z,g.)	# cc is the grid. stored outside.
  R = V(a_,  cc$z,g.)

  rval <- R@F - errors
  return(rval)
}

# you see that the exact same computations are required to get constraint and Jacobian
# i.e. you need to to V(...) twice
# there is an easy way to return both constr and jac in one function, we should definitely do that.
jacfun <- function(x,FDiffobj){

	newg <- x[ "gamma" ]
	newa <- x[ "endog" ]
  g. = FDiff(newg,'g')
  a_ = FDiff(newa,'a')

  V = SplineInt1D(xsupp,ivals)	# this can be stored on disk/memory

  R = V(cc$a,cc$z,g.)	# cc is the grid. stored outside.
  R = V(a_,  cc$z,g.)

  rval <- t(R@J)@x
  return( rval )
}



# test for turning around dgCMatrix into dgRMatrix
  A <- matrix(c(0,0,1,2,0,3,4,5,0,6,7,0,8,9,0,10),4,4,byrow=T)
  Ac = as(A,"dgCMatrix")
  Ar = as(A,"dgRMatrix")
  tAc = t(Ac)
  Ar@x
  tAc@x




