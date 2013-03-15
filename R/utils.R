#' this makes sure both Jacobians cover the same
#' variables. If they do not, a common Jocabian 
#' support is returned with appropriate values
#' @export
expandJacDomain <- function(e1, vars) {

  # check if anything is to be done
  # we need both the same set of variables
  # and the same order
  if (setequal(names(e1@vars), names(vars))) {
    if(all(names(e1@vars)==names(vars))) return(e1);
  }

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

Var <- function(name,size) {
  l = list()
  l[[name]]=size
  return(l)
}

is.fdiff <- function(x) (class(x)=='FDiff')

# this is a bit dirty, sorry....
applyColoring <- function(F) {
  if (F@coloring) {
    F@J = (F@J!=0)*1;
  }
  return(F)
}

#' variables for the Jacobian columns
#' @export
getSparseValues <- function(M,S) 
{
  r = rep(0,length(unlist(S)))
  i = 1
  for (row in 1:length(S)) {
    for (col in S[[row]]) {
      r[i] = M[row,col];
      i=i+1
    }
  }

  return(r)
}

#' transform a sparse matrix in column compressed format to the 
#' format required by ipoptr
#' @export
ipoptr.sparse <- function(A) {

  # first we get the total number of rows
  nr = nrow(A)

  nvals = diff(A@p) # this is the number of values in each column
  j = 1
  for (col in 1:length(nvals)) {
    for (i in 1:nvals[col]) {
      if ( length(S) >=(A@i[j] + 1)) {
        S[[ A@i[j] + 1 ]] = c(  S[[ A@i[j] + 1 ]]   , col )
      } else {
        S[[ A@i[j] + 1 ]] = col     
      }
      j = j+1
    }
  }

  return(S)
}


