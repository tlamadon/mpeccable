#' solves a set of constraints (equality and inequality)
#' the argument is a function cFunc that takes a parameter and returns the 
#' evaluated list of constraints, the sparse Jacobian and the constraints type
#' the second argument is an intial parameter value
#' cFunc should return the coloring if called with coloring=TRUE
#' @export
#' @example examples/mpeccable.csolve.ex.r
mpeccable.csolve <- function( cFunc, x0, vars ) {

  # Compute Jacobian Structure
  # --------------------------

  # @todo deal with inequality constraints

  # get the Jacobian from constraints
  theta = mpec.vars.collate(x0,vars,coloring=TRUE)
  res = cFunc(theta)

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
  eval_f <- function(x,private) { 
    # extract the squared errors and sum them  
    ps  = mpec.vars.collate(x,private$vars)
    return(sum( (ps[['optim.err']]@F)^2))
  }
  #  ------  objective function gradient ------
  eval_grad_f <- function(x,private) { 
    # extract the errors and return them
    ps     = mpec.vars.collate(x,private$vars)
    ranges = computeVarRanges(private$vars)
    r      = mpec.vars.init(private$vars,0)
    r[ranges[['optim.rr']]] = ps[['optim.err']]@F
    return(r)
  }
  #  ------  constraint function gradient ------
  eval_g <- function( x ,private) {    
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]             # extract the errors
    res       = private$cFunc(ps)                     # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err    # add the errors
    return(C.MSE@F)
  }
  #  ------  constraint function gradient ------
  eval_jac_g <- function( x ,private) {   
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]             # extract the errors
    res       = private$cFunc(ps)                     # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err    # add the errors

    # extract sparse structure of the JAC
    res = getSparseValues(C.MSE@J,private$eval_jac_g_structure)
    return(res)
  }

  x0.augmented = c(x0,optim.err@F)
  lb = rep(-Inf,length(x0.augmented))
  ub = rep(Inf,length(x0.augmented))
  
  private = list(vars=vars,eval_jac_g_structure=eval_jac_g_structure,cFunc=cFunc)

  length(eval_jac_g(x0.augmented,private))
  length(unlist(eval_jac_g_structure))

  # Call the optimizer
  # ------------------
  res = ipoptr( 
     x0=x0.augmented, 
     lb=lb, 
     ub=ub,
     eval_f=eval_f, 
     eval_grad_f=eval_grad_f, 
     eval_g=eval_g, 
     eval_jac_g=eval_jac_g,
     eval_jac_g_structure=eval_jac_g_structure,
     private=private)

  return(res)
}







