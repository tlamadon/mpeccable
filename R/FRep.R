#' get initial value for the functional representation 
#' @export
param0 <- function(frep,name,x0=0) {
  return(FDiff(rep(x0,attr(frep,'ng')),name  ))
}


#' define a parameter constant over the entire collocation
#' @export
fixedParam <- function(xin,name,x0=0) {
  return(FDiff(rep(x0,length(xin)),name  ))
}
