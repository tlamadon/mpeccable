require(splines)
require(Matrix)

#' implements a linear form
#' @export
#' @family frep
#' @example examples/F_Linear.ex.r
F_Linear <- function(X) {

 XX = as.matrix(X)

 ff <- function(gin) {

  # take care of the function representation
  # for sure gin is a fdiff, otherwise, makes not sense
  if (class(gin) == "FDiff") {
    F = XX %*% gin@F
    J = XX
    R = new("FDiff",F=c(F),J=Matrix(J,sparse=T),vars=gin@vars,coloring=gin@coloring)
  } else {
    stop('function representation is not a parameter, this seems odd')
  }

  return(R)
}

class(ff) = 'frep'
attr(ff,'ng') = ncol(X)
return(ff)

}