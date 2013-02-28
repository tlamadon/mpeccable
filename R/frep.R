#' get initial value for the functional representation 
#' @export
param0 <- function(frep,name,x0=0) {
  return(FDiff(rep(x0,attr(frep,'ng')),name  ))
}
