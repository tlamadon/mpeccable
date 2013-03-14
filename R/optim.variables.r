
#' creates the variable structure that will describe the
#' optimization problem
#' @param vars is a list FDiff objects
#' @export
mpec.vars.desc <- function(vars) {

  # combines the variables in a list, the order will be important
  ll = list()
  for (v in vars) {
    ll = c(ll,v@vars)
  }

  class(ll) = "mpec.vars"
  return(ll)
}

#' given a parameter value list, a variable description
#' and a variable name, returns an initialized FDiff object
#' with the correct level 
#' @export
mpec.vars.extract <- function(theta,vars,varname,coloring=FALSE) {
  if (! varname %in% names(vars)) stop(paste(varname,' is not in list of variables'));
  ranges = computeVarRanges(vars)
  return(FDiff(theta[ranges[[varname]]],varname,coloring=coloring))
}

#' create the list of params, correctly initialized
#' @export
mpec.vars.collate <- function(x0,vars,coloring=FALSE) {

  ranges = computeVarRanges(vars)
  ll = list()
  for (vs in names(ranges)) {
    ll = c(ll, mpec.vars.extract(x0,vars,vs,coloring))
  }

  names(ll) = names(vars)
  return(ll)
}