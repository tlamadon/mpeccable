#' initializes the mpec problem
#' an mpec problem is formed of an objective function (FDIFF), a set of equality 
#' constraints and a set of inequality constraints, as well as a list of variables
#' with bounds
#' @export
#' @family mpec
mpec.new <- function() {

  mp = list()

  mp$OBJ.FDiff         = 0 # objective function
  mp$CST.EQ            = 0 # constraints
  mp$CST.INEQ.FDiff    = NULL # constraints
  mp$OBJ.ABS.vars      = list()    # constraints
  mp$CST.INEQ.ub = NULL # constraints
  mp$CST.INEQ.lb = NULL # constraints
  mp$vars        = list() # variable to optimize
  mp$vars.ub     = list() # variable to optimize
  mp$vars.lb     = list() # variable to optimize
  mp$vars.level  = list() # variable to optimize
  mp$coloring    = FALSE

  class(mp) = 'mpec'
  return(mp)
}

#' creates an set of constraints in absolute value
#' to be added to an mpec problem. In practice this
#' creates 2 set of constraints and adds 1 epsilons always positive to both
#' that appear in the constraint and in the obective function
#' @export
#' @family mpec
mpec.addAbsConstraint <- function(mpec,R,name) {

  error.name = paste('err.abs.',name,sep='')

  # check if we have a value for it already
  if (error.name %in% names(mpec$vars)) {
    e.   = mpec.getVar(mpec,error.name)
  } else {
    e.   = FDiff(rep(0,dim(R)[1]),paste('err.abs.',name,sep=''))
    mpec = mpec.addVar(mpec,e.,lb=0)
  }

  # combine the constraints and append it to be negative
  R    = rbind2( R - e.  , -R - e. )
  mpec = mpec.addInequalityConstraint(mpec,R,name,ub=0,lb=-Inf)

  # append the error to the obejctive
  mpec$OBJ.ABS.vars  = c(mpec$OBJ.ABS.vars,error.name)

  return(mpec)
}

#' add an inequatllity constraint to the mpec object
#' @export
#' @family mpec
mpec.addInequalityConstraint <- function(mpec,R,name,lb=-Inf,ub=Inf) {

  if (length(R@F)>length(lb)) lb = rep(lb,length(R@F));
  if (length(R@F)>length(ub)) ub = rep(ub,length(R@F));

  # append inequality constraint
  if (is.null(mpec$CST.INEQ.FDiff)) {
    mpec$CST.INEQ.FDiff = R
    mpec$CST.INEQ.ub = ub
    mpec$CST.INEQ.lb = lb
  } else {
    mpec$CST.INEQ.FDiff= rbind2( mpec$CST.INEQ.FDiff , R)
    mpec$CST.INEQ.ub = c(mpec$CST.INEQ.ub,ub)
    mpec$CST.INEQ.lb = c(mpec$CST.INEQ.lb,lb)
  }

  # append the variables
  mpec$vars = mergevars(mpec$vars,R@vars)

  # return the mpec object
  return(mpec)
}

#' Creates the object that stores information about the different
#' variables used in the optimzation problem. It takes an fdiff argument
#' and possible lower and upper bounds information
#' @param var  a FDiff object
#' @param lb the lower bound, either one value or a vector of same length as FDiff
#' @param ub the upper bound, either one value or a vector of same length as FDiff
#' @export
#' @family mpec.vars
mpec.addVar <- function(mpec,F,lb=-Inf,ub=Inf) {
  # makes sure size is correct
  if (length(F@F)>length(lb)) lb = rep(lb,length(F@F))
  if (length(F@F)>length(ub)) ub = rep(ub,length(F@F))

  # append the variables
  mpec$vars = mergevars(mpec$vars,F@vars)

  # add the levels
  nn <- names(F@vars)[1]
  mpec$vars.level[[  nn   ]] = F@F

  # add the bounds
  mpec$vars.ub[[nn]] = ub
  mpec$vars.lb[[nn]] = lb

  return(mpec)
}

#' given a parameter value list, a variable description
#' and a variable name, returns an initialized FDiff object
#' with the correct level 
#' @export
mpec.getVar <- function(mpec,varname) {
  if (! varname %in% names(mpec$vars)) stop(paste(varname,' is not in list of variables'));
  return(FDiff(mpec$vars.level[[varname]],varname,coloring=mpec$coloring))
}

#' @export
mpec.getVarRange <- function(mpec,varname) {
  if (! varname %in% names(mpec$vars)) stop(paste(varname,' is not in list of variables'));
  ranges = computeVarRanges(mpec$vars)
  return(ranges[[varname]])
}

#' creates a list of FDiff objects from a vector of values. This is
#' useful when processing the parameters given by the optimizer before
#' sending it to the function that computes the constraint
#' @export
mpec.setVarsFromVector <- function(mpec,x,coloring=FALSE) {
  ranges = computeVarRanges(mpec$vars)
  for (vs in names(mpec$vars)) {
    mpec$vars.level[[vs]] = x[ ranges[[vs]] ]
  }
  return(mpec)
}
#' @export
mpec.getVarsAsVector <- function(mpec) {
  res = c()
  for (n in names(mpec$vars)) {
    res = c(res,mpec$vars.level[[n]])
  }
  return(res)
}

#' @export
mpec.getBoundsAsVector <- function(mpec) {
  lb = c()
  ub = c()
  for (n in names(mpec$vars)) {
    lb = c(lb,mpec$vars.lb[[n]])
    ub = c(ub,mpec$vars.ub[[n]])
  }
  return(list(ub=ub,lb=lb))
}

#' @export
mpec.getObjective <- function(mpec) {
  # make sure Jac support is correct
  if (class(mpec$OBJ.FDiff)=='FDiff' || length(mpec$OBJ.ABS.vars)>0  ) {
    OBJ = mpec$OBJ
    for (v in mpec$OBJ.ABS.vars) {
      OBJ = OBJ + sum(mpec.getVar(mpec,v))
    }
    OBJ = expandJacDomain( OBJ ,mpec$vars)
  } else {
    return(0)
  }
  return(OBJ)
}

#' @export
mpec.reset <- function(mpec) {
  mpec$OBJ         = 0 # objective function
  mpec$CST.EQ      = 0 # constraints
  mpec$CST.INEQ.FDiff    = NULL # constraints
  mpec$vars.level  = list() # variable to optimize
  mpec$coloring    = FALSE
  return(mpec)
}

#' @export
mpec.vectorToList <- function(mpec,x) {
  mpec = mpec.setVarsFromVector(mpec,x)
  res = list()
  for (vs in names(mpec$vars)) {
    res[[vs]] = mpec.getVar(mpec,vs)
  }
  return(res)
}

