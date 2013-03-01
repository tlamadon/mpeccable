#' fdiff
#'
#' FDiff is a pacakge to facilitate the solving of functional equation
#' with control variables. The method is built around the general finite
#' element approach. 
#'
#' @import Matrix
#' @docType package
#' @name fdiff
#' @rd-name fdiff


# rutn the following command to reload dev package
# roxygenise('~/git/fdiff');install('~/git/fdiff')

# creating a class for mutlivariate differentiation
# the structure is only a vector of values and a Jocabian

# Main class
# ==========
# Defining the main dif class that will sotre the 
# values and the hessian

#' @export
setClass("FDiff", 
    representation(F = "numeric", J = "sparseMatrix", vars="list", coloring='logical'))

#'@export
setMethod("+", c("FDiff","numeric"), function(e1,e2) { e1@F = e1@F+e2; e1}) 
#'@export
setMethod("-", c("FDiff","numeric"), function(e1,e2) { e1@F = e1@F-e2;e1}) 
#'@export
setMethod("/", c("FDiff","numeric"), function(e1,e2) { e1@F = e1@F/e2;e1@J = e1@J/e2 ; e1}) 
#'@export
setMethod("*", c("FDiff","numeric"), function(e1,e2) { e1@F = e1@F*e2;e1@J = e1@J*e2 ; e1}) 
#'@export
setMethod("+", c("numeric","FDiff"), function(e1,e2) { e2@F = e2@F+e1;e2}) 
#'@export
setMethod("-", c("numeric","FDiff"), function(e1,e2) { e2@F = e2@F-e1;e2}) 
#'@export
setMethod("/", c("numeric","FDiff"), function(e1,e2) { e2@F = e2@F/e1;e2@J = e2@J/e1 ; e2}) 
#'@export
setMethod("*", c("numeric","FDiff"), function(e1,e2) { e2@F = e2@F*e1;e2@J = e2@J*e1 ; e2}) 
#'@export
setMethod("%*%", c("matrix","FDiff"), function(x,y) { y@F = as.numeric(x %*% y@F); y@J = Matrix(x,sparse=TRUE)%*%y@J ; y}) 

#'@export
setMethod("+", c("FDiff","FDiff"), 
  function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    e1@F = e1@F + e2@F
    e1@J = e1@J + e2@J
    e1
  })
#'@export
setMethod("-", c("FDiff","FDiff"), 
  function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    e1@F = e1@F - e2@F
    e1@J = e1@J - e2@J
    e1
  })
#'@export
setMethod("*", c("FDiff","FDiff"), 
  function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    e1@J = Matrix(diag(e2@F),sparse=T) %*% e1@J + Matrix(diag(e1@F),sparse=T) %*% e2@J
    e1@F = e1@F * e2@F
    e1
  })
#'@export
setMethod("/", c("FDiff","FDiff"), 
  function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    e1@J = Matrix(diag( (e2@F)^2 ),sparse=T) %*% ( - Matrix(diag(e2@F),sparse=T) %*% e1@J + Matrix(diag(e1@F),sparse=T) %*% e2@J) 
    e1@F = e1@F / e2@F
    e1
  })
#'@export
setMethod("^", c("FDiff","numeric"), 
  function (e1, e2) {
    e1@J = e2 * Matrix(diag((e1@F)^(e2-1)),sparse=T) %*% e1@J 
    e1@F = (e1@F)^e2
    e1
  })
# initializes a variable it's NxN, F is x and J is diag(1)
#'@export
setMethod("initialize", "FDiff", function(.Object,F,J,vars) { 
  .Object@F <- F
  .Object@J <- J
  .Object@vars <- vars
  .Object@coloring=FALSE
  .Object
  })

# initializes a variable it's NxN, F is x and J is diag(1)
#'@export
setMethod("print", "FDiff", function(x) { 
  cat('FDiff object (', length(x@F) , 'x' , ncol(x@J), ') dense:', length(x@J@x)/prod(dim(x@J))  ,'% \n')
  })

# initializes a variable it's NxN, F is x and J is diag(1)
#' @export
setGeneric("pMax", function(e1,e2) {
  standardGeneric("pMax")
})

# initializes a variable it's NxN, F is x and J is diag(1)
#' @export
setGeneric("appendJac", function(x,J,vs) {
    standardGeneric("appendJac")
})


#' like pmax, take the max between two values. It works on vector and
#' computes the Jacobian correctly
#' @docType methods
#' @rdname fdiff-methods
setMethod("pMax", c("FDiff","numeric"), function(e1,e2) { 
  I = e1@F < e2
  if(!any(I)) return(e1);

  e1@J[I,] = 0
  if (length(e2)>1) {
    e1@F[I]  = e2[I]
  } else {
    e1@F[I]  = e2
  }
  e1
})

#' append to Jac
#' @docType methods
#' @rdname fdiff-methods
#' @export
setMethod("appendJac", c("FDiff","dsCMatrix","list"), function(x,J,vs) { 
  x@J = cBind(x@J,J)
  x@vars = c(x@vars,vs)
  x
})

#' like rBind, combines the levels and jacobian, however
#' it makes the variables of the Jacobian consistent 
#' @docType methods
#' @rdname fdiff-methods
#' @export
setMethod("rbind2", c("FDiff","FDiff"), function(x,y) { 
  vars = mergevars(x@vars,y@vars)
  x   = expandJacDomain(x,vars)
  y   = expandJacDomain(y,vars)
  x@F = c(x@F,y@F)
  x@J = rBind(x@J,y@J)
  x@vars = vars
  x
})

#' allows accessing the levels directly, this is convenient
#' within the code
#' @export
setMethod(
  f= "[",
  signature="FDiff",
  definition=function(x,i){
    return(x@F[i])
  }
)

#' creates a variable to be tracked by the computation of the Jacobian
#' @param x a vector of current values
#' @param name a anme that uniquly identify this variable in the overall
#' @export
FDiff <- function(x, name) {
  N = length(x)
  vars = list()
  vars[name] = N
  new("FDiff",x,Matrix(diag(rep(1,N)),sparse=TRUE),vars)
}




setMethod("Ops", c("FDiff","FDiff"), function(e1,e2) {})

# Simple function representation
# ==============================

# we want to represent a function in a finite subspace
# and evaluate it at the collocation. 
# a function representation is an R function that takes
# a collocation and a parameter vector and returns 
# the correspong FDiff object.

# let's first look at a one dimensional spline






 

