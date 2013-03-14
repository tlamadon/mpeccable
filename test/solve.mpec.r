

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

  # @todo bounds on parameters
  # @todo deal with inequality constraints

  # get the Jacobian from constraints
  theta = mpec.vars.collate(x0,vars,coloring=TRUE)
  res = cFunc(ps)

  # create the error for each equality constraint, not for inequality
  C.MSE = res[['C.MSE']]

  # I create the error and append it
  optim.err = FDiff(rep(0,dim(C.MSE)[1]),'optim.err') 
  C.MSE = C.MSE + optim.err
  constraint_lb = 0
  constraint_ub = 0

  # append optim.err to list of vars
  vars[['optim.err']] = dim(optim.err)[1]

  # exctract the sparse structure of the jacobian
  eval_jac_g_structure <- make.sparse( C.MSE@J!=0) # needs to 

  # Create the objective function and the constraint function for ipopt
  # -------------------------------------------------------------------

  #  ------  objective function ------
  eval_f <- function(x) { 
    # extract the squared errors and sum them  
    ps  = mpec.vars.collate(x,vars)
    return(sum( (ps[['optim.err']]@F)^2))
  }
  #  ------  objective function gradient ------
  eval_grad_f <- function(x) { 
    # extract the errors and return them
    ps  = mpec.vars.collate(x,vars)
    return(2*ps[['optim.err']]@F)
  }
  #  ------  constraint function gradient ------
  eval_g <- function( x ) {    
    ps        = mpec.vars.collate(x,vars)    # extract the parameters
    optim.err = ps[['optim.err']]             # extract the errors
    res       = cFunc(ps)                     # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err    # add the errors
    return(C.MSE@F)
  }
  #  ------  constraint function gradient ------
  eval_jac_g <- function( x ) {    
    ps        = mpec.vars.collate(x,vars)    # extract the parameters
    optim.err = ps[['optim.err']]             # extract the errors
    res       = cFunc(ps)                     # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err    # add the errors

    # extract sparse structure of the JAC
    return(t(C.MSE@J)@x)
  }

  x0.augmented = c(x0,optim.err@F)

  # Call the optimizer
  # ------------------
  res = ipoptr( 
     x0=x0, 
     eval_f=eval_f, 
     eval_grad_f=eval_grad_f, 
     eval_g=eval_g, 
     eval_jac_g=eval_jac_g,
     eval_jac_g_structure=eval_jac_g_structure,
     constraint_lb=constraint_lb,
     constraint_ub=constraint_ub,
     ipoptr_environment = environment())

  return(res)
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






