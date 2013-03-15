#' solves a set of constraints (equality and inequality)
#' the argument is a function cFunc that takes a parameter and returns the 
#' evaluated list of constraints, the sparse Jacobian and the constraints type
#' the second argument is an intial parameter value
#' cFunc should return the coloring if called with coloring=TRUE
#' @export
#' @example examples/mpeccable.csolve.ex.r
mpeccable.csolve <- function( objFunc=function(params) {list(OBJ=0)}, cFunc, x0, vars ,opts ) {

  # Compute Jacobian Structure
  # --------------------------

  # @todo deal with inequality constraints

  # get the Jacobian from constraints
  theta = mpec.vars.collate(x0,vars,coloring=TRUE)
  res = cFunc(theta)
  reso= objFunc(theta)

  # create the error for each equality constraint, not for inequality
  C.MSE = res[['C.MSE']]
  OBJ   = reso[['OBJ']]

  # I create the error and append it
  optim.err = FDiff(rep(0,dim(C.MSE)[1]),'optim.err') 
  C.MSE = C.MSE + optim.err
  OBJ   = OBJ + sum(optim.err)

  # append optim.err to list of vars
  vars[['optim.err']] = dim(optim.err)[1]

  # exctract the sparse structure of the jacobian
  eval_jac_g_structure <- make.sparse( C.MSE@J!=0) # needs to 
  constraint_lb = rep(0,length(eval_jac_g_structure))
  constraint_ub = rep(0,length(eval_jac_g_structure))

  # Create the objective function and the constraint function for ipopt
  # -------------------------------------------------------------------

  #  ------  objective function ------
  eval_f <- function(x,private) { 
    ps   = mpec.vars.collate(x,private$vars)
    res  = private$objFunc(ps)
    optim.err = ps[['optim.err']]                    # extract the errors
    OBJ  = res[['OBJ']] + sum(optim.err^2)    # add the errors
    return(OBJ@F)
  }
  #  ------  objective function gradient ------
  eval_grad_f <- function(x,private) { 
    ps   = mpec.vars.collate(x,private$vars)
    res  = private$objFunc(ps)
    optim.err = ps[['optim.err']]                    # extract the errors
    OBJ  = res[['OBJ']] + sum(optim.err^2)    # add the errors
    return(as.numeric(OBJ@J))
  }
  #  ------  constraint function gradient ------
  eval_g <- function( x ,private) {    
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]                    # extract the errors
    res       = private$cFunc(ps)                    # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err           # add the errors
    return(C.MSE@F)
  }
  #  ------  constraint function gradient ------
  eval_jac_g <- function( x ,private) {   
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]                    # extract the errors
    res       = private$cFunc(ps)                    # evaluate the constraints
    C.MSE     = res[['C.MSE']] + optim.err           # add the errors

    # extract sparse structure of the JAC
    res = getSparseValues(C.MSE@J,private$eval_jac_g_structure)
    return(res)
  }

  x0.augmented = c(x0,optim.err@F)
  lb = rep(-Inf,length(x0.augmented))
  ub = rep(Inf ,length(x0.augmented))
  # adding the contraint that the errors are positive
  ranges = computeVarRanges(vars)
  lb[ranges[['optim.err']]] = 0

  private = list(vars=vars,eval_jac_g_structure=eval_jac_g_structure,cFunc=cFunc,objFunc=objFunc)

  #length(eval_jac_g(x0.augmented,private))
  #length(unlist(eval_jac_g_structure))
  #length(eval_grad_f(x0.augmented,private))

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
     constraint_lb=constraint_lb,
     constraint_ub=constraint_ub,
     opts=opts,
     private=private)

  theta.opt = mpec.vars.collate(res$solution,private$vars)
  res$solution = theta.opt

  return(res)
}


#' solves a set of constraints (equality and inequality)
#' the argument is a function cFunc that takes a parameter and returns the 
#' evaluated list of constraints, the sparse Jacobian and the constraints type
#' the second argument is an intial parameter value
#' cFunc should return the coloring if called with coloring=TRUE
#' @export
mpeccable.asolve <- function( objFunc=function(params) {list(OBJ=0)}, cFunc, x0, vars ,opts ) {

  # Compute Jacobian Structure
  # --------------------------

  # @todo deal with inequality constraints

  # get the Jacobian from constraints
  theta = mpec.vars.collate(x0,vars,coloring=TRUE)
  res = cFunc(theta)
  reso= objFunc(theta)

  # create the error for each equality constraint, not for inequality
  C.ABS = res[['C.ABS']]
 
  # I create the error and append it
  C.ABS.p     =  C.ABS 
  C.ABS.n     = -C.ABS 
  C.ABS       = rbind2(C.ABS.p,C.ABS.n)
  optim.err   = FDiff(rep(0,dim(C.ABS)[1]),'optim.err') 
  C.ABS       = C.ABS - optim.err
  constraint_lb = rep(-Inf   ,dim(C.ABS)[1])
  constraint_ub = rep(0      ,dim(C.ABS)[1])
 
  # append optim.err to list of vars
  vars[['optim.err']] = dim(optim.err)[1]

  # exctract the sparse structure of the jacobian
  eval_jac_g_structure <- make.sparse( C.ABS@J!=0) # needs to 

  # Create the objective function and the constraint function for ipopt
  # -------------------------------------------------------------------

  #  ------  objective function ------
  eval_f <- function(x,private) { 
    ps   = mpec.vars.collate(x,private$vars)
    res  = private$objFunc(ps)
    optim.err = ps[['optim.err']]           # extract the errors
    OBJ  = res[['OBJ']] + sum(optim.err)    # add the errors
    return(OBJ@F)
  }
  #  ------  objective function gradient ------
  eval_grad_f <- function(x,private) { 
    ps   = mpec.vars.collate(x,private$vars)
    res  = private$objFunc(ps)
    optim.err = ps[['optim.err']]           # extract the errors
    OBJ  = res[['OBJ']] + sum(optim.err)    # add the errors
    vv <- private$vars
    class(vv)<-"list"
    OBJ = expandJacDomain(OBJ,vv) # make sure the support of variable is correct
    return(as.numeric(OBJ@J))
  }
  #  ------  constraint function gradient ------
  eval_g <- function( x ,private) {    
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]                    # extract the errors
    res       = private$cFunc(ps)                    # evaluate the constraints
    C.ABS     = rbind2(res[['C.ABS']],-res[['C.ABS']]) - optim.err           # add the errors
    return(C.ABS@F)
  }
  #  ------  constraint function gradient ------
  eval_jac_g <- function( x ,private) {   
    ps        = mpec.vars.collate(x,private$vars)    # extract the parameters
    optim.err = ps[['optim.err']]                    # extract the errors
    res       = private$cFunc(ps)                    # evaluate the constraints
    C.ABS     = rbind2(res[['C.ABS']],-res[['C.ABS']]) - optim.err # add the errors

    # extract sparse structure of the JAC
    res = getSparseValues(C.ABS@J,private$eval_jac_g_structure)
    return(res)
  }

  x0.augmented = c(x0,optim.err@F)
  lb = rep(-Inf,length(x0.augmented))
  ub = rep(Inf,length(x0.augmented))
  
  private = list(vars=vars,eval_jac_g_structure=eval_jac_g_structure,cFunc=cFunc,objFunc=objFunc)

  #length(eval_jac_g(x0.augmented,private))
  #length(unlist(eval_jac_g_structure))
  #length(eval_grad_f(x0.augmented,private))

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
     constraint_lb=constraint_lb,
     constraint_ub=constraint_ub,
     opts=opts,
     private=private)

  theta.opt = mpec.vars.collate(res$solution,private$vars)
  res$solution = theta.opt

  return(res)
}





