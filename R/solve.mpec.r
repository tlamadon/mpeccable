#' solves a set of constraints (equality and inequality)
#' the argument is a function cFunc that takes a parameter and returns the 
#' evaluated list of constraints, the sparse Jacobian and the constraints type
#' the second argument is an intial parameter value
#' cFunc should return the coloring if called with coloring=TRUE
#' @export
#' @example examples/example-MpecSolve.r
mpec.solve <- function( objFunc=function(mpec) {mpec}, cFunc, x0, mpec ,opts ) {

  # Compute Jacobian Structure
  # --------------------------

  # intialize variable with x0
  mpec.color = mpec.setVarsFromVector(mpec,x0)
  mpec.color$coloring=TRUE

  # evaluate the constraints
  mpec.color = cFunc(mpec.color)
  x0 = mpec.getVarsAsVector(mpec.color)

  # get jacobian of constraints
  eval_jac_g_structure <- ipoptr.sparse( mpec.color$CST.INEQ.FDiff@J ) 
  constraint_lb = mpec.color$CST.INEQ.lb
  constraint_ub = mpec.color$CST.INEQ.ub

  # set the parameter constratints
  bounds = mpec.getBoundsAsVector(mpec.color)
  lb = bounds$lb
  ub = bounds$ub
  mpec.color= mpec.reset(mpec.color)

  # Create the objective function and the constraint function for ipopt
  # -------------------------------------------------------------------
  private = list(mpec=mpec.color, cFunc=cFunc, objFunc=objFunc,eval_jac_g_structure=eval_jac_g_structure,evals=list(),last_x=0)

  last_x    = -1
  last_mpec = NA
  hit = 0
  miss =0

  #  ------  objective function ------
  eval_f <- function(x,private) { 
    if (all(x==last_x)) {
      mpec = last_mpec
      hit <<- hit +1
    } else {
      mpec = mpec.setVarsFromVector(private$mpec,x)
      mpec = private$objFunc(mpec)
      last_x <<- x
      last_mpec <<- mpec
      miss <<- miss +1
    }
    res  = as.numeric(mpec.getObjective(mpec)@F)
    if (any(!is.finite(res))) stop("error in eval_grad_f");
    return(res) 
  }
  #  ------  objective function gradient ------
  eval_grad_f <- function(x,private) { 
    if (all(x==last_x)) {
      mpec = last_mpec
      hit <<- hit +1
    } else {
      mpec = mpec.setVarsFromVector(private$mpec,x)
      mpec = private$objFunc(mpec)
      last_x <<- x
      last_mpec <<- mpec
      miss <<- miss +1
    }
    mpec  = mpec.setVarsFromVector(private$mpec,x)
    mpec  = private$objFunc(mpec)
    res   = as.numeric(mpec.getObjective(mpec)@J)
    if (any(!is.finite(res))) stop("error in eval_grad_f");
    return(res)
  }
  #  ------  constraint function gradient ------
  eval_g <- function( x ,private) {  
    if (all(x==last_x)) {
      mpec = last_mpec
      hit <<- hit +1
    } else {
      mpec = mpec.setVarsFromVector(private$mpec,x)
      mpec = private$objFunc(mpec)
      last_x <<- x
      last_mpec <<- mpec
      miss <<- miss +1
    }
    mpec  = mpec.setVarsFromVector(private$mpec,x)
    mpec  = private$cFunc(mpec)                    # evaluate the constraints
    res   = mpec$CST.INEQ.FDiff@F
    if (any(!is.finite(res))) stop("error in eval_g");
    return(res)
  }
  #  ------  constraint function gradient ------
  eval_jac_g <- function( x ,private) { 
     if (all(x==last_x)) {
      mpec = last_mpec
      hit <<- hit +1
    } else {
      mpec = mpec.setVarsFromVector(private$mpec,x)
      mpec = private$objFunc(mpec)
      last_x <<- x
      last_mpec <<- mpec
      miss <<- miss +1
    }
    mpec  = mpec.setVarsFromVector(private$mpec,x)
    mpec  = private$cFunc(mpec)                    # evaluate the constraints
 
    # extract sparse structure of the JAC
    res = getSparseValues(mpec$CST.INEQ.FDiff@J,private$eval_jac_g_structure)
    if (any(!is.finite(res))) stop("error in jac");
    return(res)
  }

  # Call the optimizer
  # ------------------
  res = ipoptr( 
     x0=x0, 
     lb=lb, 
     ub=ub,
     eval_f=eval_f, 
     eval_grad_f=eval_grad_f, 
     eval_g=eval_g, 
     eval_jac_g=eval_jac_g,
     eval_jac_g_structure=eval_jac_g_structure,
     constraint_lb=constraint_lb,
     constraint_ub=constraint_ub,
     opts=opts,
     private=private,
     ipoptr_environment=environment())

  theta.opt = mpec.vectorToList(mpec.color,res$solution)
  res$solution = theta.opt

  return(res)
}



