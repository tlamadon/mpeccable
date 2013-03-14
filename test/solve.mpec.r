

# from user point of view
#' solves a set of constraints (equality and inequality)
#' the argument is a function cFunc that takes a parameter and returns the 
#' evaluated list of constraints, the sparse Jacobian and the constraints type
#' the second argument is an intial parameter value
#' cFunc should return the coloring if called with coloring=TRUE
#' @export
#' @example examples/mpeccable.csolve.r
mpeccable.csolve <- function( cFunc, x0, vars ) {

  # Compute Jacobian Structure
  # --------------------------

  # get the Jacobian from constraints
  theta = mpec.vars.collate(x0,vars,coloring=TRUE)
  res = cFunc(ps)

  # create the error for each equality constraint, not for inequality
  I = res$C == 'MSE'

  








}




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

# ipopt wrapper
ipopt.const <- function(x) {

  errors = ggetFromX

  # remove errors and call the user function
  R = constfun(x)

  # add errors to it
  rval <- R@F - errors

  # append to Jacobian

  return(   )
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




