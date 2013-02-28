#' this makes sure both Jacobians cover the same
#' variables. If they do not, a common Jocabian 
#' support is returned with appropriate values
expandJacDomain <- function(e1, vars) {

  # check if anything is to be done
  if (setequal(names(e1@vars), names(vars))) return(e1);

  # we merge the variable lists
  range1 = computeVarRanges(e1@vars)
  rangeall = computeVarRanges(vars)

  J1 = Matrix(0,length(e1@F),sum(unlist(vars)) , sparse=T)
  for (v in names(vars)) {
    if (v %in% names(range1)) {
      J1[,rangeall[[v]]] = e1@J[,range1[[v]]]
    }
  }

  e1@J = J1
  e1@vars = vars
  return(e1)
}

#' computes the indices associated with each names 
#' variables for the Jacobian columns
#' @export
computeVarRanges <- function(vars) {
  r = list()
  s=1
  for (v in names(vars)) {
    r[[v]] = s:(s + vars[[v]]-1)
    s = s + vars[[v]] 
  }
  return(r)
}

#' merge two list but don't keep duplicates
mergevars <- function(v1,v2) {
  # we add to v1 the variable from v2 not in v1
  nn = setdiff(names(v2),names(v1))
  if (length(nn)>0) {
    v1[nn] = v2[nn]
  }
  v1
}