# TODO: diag( x ) returns an error
# if x is of length 1, and it is not equal to an integer.
# If it is an integer, then a diagonal matrix with ones of size integer is returned.
# This is not what we want, so we have to specify ncol, nrow in all diag statements.

# rutn the following command to reload dev package
# roxygenise('~/git/fdiff');install('~/git/fdiff')

# creating a class for mutlivariate differentiation
# the structure is only a vector of values and a Jocabian

# Main class
# ==========
# Defining the main dif class that will sotre the 
# values and the hessian
require(Matrix)

#' @export
setClass("FDiff", 
    representation(
        F                = "numeric", 
        J                = "sparseMatrix",
        vars             = "list", 
        coloring         = "logical" ) )

# Constructor.
# initializes a variable it's NxN, F is x and J is diag(1)
#'@export
setMethod("initialize", "FDiff", function(.Object,F,J,vars,name,coloring=FALSE) { 
    .Object@F <- F
    .Object@J <- J
    .Object@vars <- vars
    .Object@coloring=coloring
    .Object
})

# Print FDiff object.
#'@export
setMethod("print", "FDiff", function(x) { 
    if (length(x@F)!=0){
        cat('FDiff object (', length(x@F) , 'x' , ncol(x@J), ') dense:', length(x@J@x)/prod(dim(x@J))  ,'% \n')
    } else {
        cat("empty FDiff object\n")
    }
})

##
#
# Unary operators.
#
##

# Unary operator(+)
#'@export
setMethod("+", c("FDiff","missing"), function(e1,e2) {
    e1@F = e1@F
    e1@J = e1@J
    return( e1 )
})

# Unary operator(-)
#'@export
setMethod("-", c("FDiff","missing"), function(e1,e2) {
    e1@F = -e1@F
    e1@J = -e1@J
    return( e1 )
})

##
#
# Binary operators("FDiff","numeric")
#
##

# By default, the operator is not defined.
setMethod("Ops", c("FDiff","numeric"), function(e1,e2) {
    stop('Operator on FDiff and numeric not yet implemented.')
})

#'@export
setMethod("+", c("FDiff","numeric"), function(e1,e2) {
    if ( length(e1@F) == length(e2) ) {
        e1@F = e1@F + e2
    } else if ( length(e1@F) == 1 & length(e2) > 1 ) {
        # Duplicate single row e1@J to get length(e2) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), length(e2) ), ncol=ncol(e1@J), nrow=length(e2), byrow=TRUE, sparse=TRUE )
        e1@F = e1@F + e2
        e1@J = Jac1
    } else if ( length(e1@F) > 1 & length(e2) == 1 ) {
        e1@F = e1@F + e2
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e1) ) 
})

#'@export
setMethod("-", c("FDiff","numeric"), function(e1,e2) {
    if ( length(e1@F) == length(e2) ) {
        e1@F = e1@F - e2
    } else if ( length(e1@F) == 1 & length(e2) > 1 ) {
        # Duplicate single row e1@J to get length(e2) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), length(e2) ), ncol=ncol(e1@J), nrow=length(e2), byrow=TRUE, sparse=TRUE )
        e1@F = e1@F - e2
        e1@J = Jac1
    } else if ( length(e1@F) > 1 & length(e2) == 1 ) {
        e1@F = e1@F - e2
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e1) )
})

#'@export
setMethod("/", c("FDiff","numeric"), function(e1,e2) {
    if ( length(e1@F) == length(e2) ) {
        e1@F = e1@F / e2
        e1@J = e1@J / e2
    } else if ( length(e1@F) == 1 & length(e2) > 1 ) {
        # Duplicate single row e1@J to get length(e2) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), length(e2) ), ncol=ncol(e1@J), nrow=length(e2), byrow=TRUE, sparse=TRUE )
        e1@F = e1@F / e2
        e1@J = Jac1 / e2
    } else if ( length(e1@F) > 1 & length(e2) == 1 ) {
        e1@F = e1@F / e2
        e1@J = e1@J / e2
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e1) )
})

#'@export
setMethod("*", c("FDiff","numeric"), function(e1,e2) {
    if ( length(e1@F) == length(e2) ) {
        e1@F = e1@F * e2
        e1@J = e1@J * e2
    } else if ( length(e1@F) == 1 & length(e2) > 1 ) {
        # Duplicate single row e1@J to get length(e2) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), length(e2) ), ncol=ncol(e1@J), nrow=length(e2), byrow=TRUE, sparse=TRUE )
        e1@F = e1@F * e2
        e1@J = Jac1 * e2
    } else if ( length(e1@F) > 1 & length(e2) == 1 ) {
        e1@F = e1@F * e2
        e1@J = e1@J * e2
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e1) )
})

#'@export
setMethod("^", c("FDiff","numeric"), function (e1, e2) {
    if ( length(e1@F) == length(e2) ) {
        e1@J = e2 * Diagonal( x = (e1@F)^(e2-1) ) %*% e1@J 
        e1@F = (e1@F)^e2
    } else if ( length(e1@F) == 1 & length(e2) > 1 ) {
        # Duplicate single row e1@J to get length(e2) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), length(e2) ), ncol=ncol(e1@J), nrow=length(e2), byrow=TRUE, sparse=TRUE )
        e1@J = e2 * (e1@F)^(e2-1) * Jac1
        e1@F = (e1@F)^e2
    } else if ( length(e1@F) > 1 & length(e2) == 1 ) {
        e1@J = e2 * Diagonal( x = (e1@F)^(e2-1) ) %*% e1@J 
        e1@F = (e1@F)^e2
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }
    
    return( applyColoring(e1) )
})

##
#
# Binary operators("numeric","FDiff")
#
##

# By default, the operator is not defined.
setMethod("Ops", c("numeric","FDiff"), function(e1,e2) {
    stop('Operator on numeric and FDiff not yet implemented.')
})

#'@export
setMethod("+", c("numeric","FDiff"), function(e1,e2) {
    if ( length(e1) == length(e2@F) ) {
        e2@F = e2@F + e1
    } else if ( length(e1) == 1 & length(e2@F) > 1 ) {
        e2@F = e2@F + e1
    } else if ( length(e1) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), length(e1) ), ncol=ncol(e2@J), nrow=length(e1), byrow=TRUE, sparse=TRUE )
        e2@J = Jac2
        e2@F = e2@F + e1
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e2) )
})

#'@export
setMethod("-", c("numeric","FDiff"), function(e1,e2) {
    if ( length(e1) == length(e2@F) ) {
        e2@J = -e2@J
        e2@F = e1 - e2@F
    } else if ( length(e1) == 1 & length(e2@F) > 1 ) {
        e2@J = -e2@J
        e2@F = e1 - e2@F
    } else if ( length(e1) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), length(e1) ), ncol=ncol(e2@J), nrow=length(e1), byrow=TRUE, sparse=TRUE )
        e2@J = -Jac2
        e2@F = e1 - e2@F
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e2) )
})

#'@export
setMethod("/", c("numeric","FDiff"), function(e1,e2) {
    # Division of scalar by a sparse matrix (elementwise) by definition results in a 
    # dense matrix, because all zero elements are now Inf. In order to comply with the
    # definition of the Jacobian in FDiff, we explicitly convert the matrix to a sparse
    # matrix.
    if ( length(e1) == length(e2@F) ) {
        e2@J = Diagonal( x = -e1/((e2@F)^2) ) %*% e2@J
        e2@F = e1 / e2@F
    } else if ( length(e1) == 1 & length(e2@F) > 1 ) {
        e2@J = Diagonal( x = -e1/((e2@F)^2) ) %*% e2@J
        e2@F = e1 / e2@F
    } else if ( length(e1) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), length(e1) ), ncol=ncol(e2@J), nrow=length(e1), byrow=TRUE, sparse=TRUE )
        e2@J = Diagonal( x = -e1/((e2@F)^2) ) %*% Jac2
        e2@F = e1 / e2@F
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e2) )
})

#'@export
setMethod("*", c("numeric","FDiff"), function(e1,e2) {
    if ( length(e1) == length(e2@F) ) {
        e2@J = e1 * e2@J
        e2@F = e1 * e2@F
    } else if ( length(e1) == 1 & length(e2@F) > 1 ) {
        e2@J = e1 * e2@J
        e2@F = e1 * e2@F
    } else if ( length(e1) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), length(e1) ), ncol=ncol(e2@J), nrow=length(e1), byrow=TRUE, sparse=TRUE )
        e2@J = e1 * Jac2
        e2@F = e1 * e2@F
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e2) )
})

##
#
# Binary operators("FDiff","FDiff")
#
##

# By default, the operator is not defined.
setMethod("Ops", c("FDiff","FDiff"), function(e1,e2) {
    stop('Operator on FDiff and FDiff not yet implemented.')
})

#'@export
setMethod("+", c("FDiff","FDiff"), function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    if ( length(e1@F) == 1 ) {
        # Duplicate single row e1@J to get nrow(e2@F) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), nrow(e2@J) ), ncol=ncol(e1@J), nrow=nrow(e2@J), byrow=TRUE, sparse=TRUE )
    } else {
        Jac1 = e1@J
    }
    if ( length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1@F) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), nrow(e1@J) ), ncol=ncol(e2@J), nrow=nrow(e1@J), byrow=TRUE, sparse=TRUE )
    } else {
        Jac2 = e2@J
    }
    e1@F = e1@F + e2@F
    e1@J = Jac1 + Jac2
    return( applyColoring(e1) )
})

#'@export
setMethod("-", c("FDiff","FDiff"), function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    if ( length(e1@F) == 1 ) {
        # Duplicate single row e1@J to get nrow(e2@F) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), nrow(e2@J) ), ncol=ncol(e1@J), nrow=nrow(e2@J), byrow=TRUE, sparse=TRUE )
    } else {
        Jac1 = e1@J
    }
    if ( length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1@F) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), nrow(e1@J) ), ncol=ncol(e2@J), nrow=nrow(e1@J), byrow=TRUE, sparse=TRUE )
    } else {
        Jac2 = e2@J
    }
	# we check for coloring because Jac1-Jac2 = 1 - 1 = 0 and we want 1.
	# @FIXME this should be done everywhere
	if (e1@coloring){
		e1@F = e1@F - e2@F
		e1@J = Jac1 | Jac2
	} else {
		e1@F = e1@F - e2@F
		e1@J = Jac1 - Jac2
	}
    return( applyColoring(e1) )
})

#'@export
setMethod("*", c("FDiff","FDiff"), function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    
    # Perform different calculation depending on whether e1 or e2 is scalar.
    if ( length(e1@F) == length(e2@F) ) {
        e1@J = Diagonal( x = e2@F ) %*% e1@J + 
               Diagonal( x = e1@F ) %*% e2@J
        e1@F = e1@F * e2@F
    } else if ( length(e1@F) == 1 & length(e2@F) > 1 ) {
        # Duplicate single row e1@J to get nrow(e2@F) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), nrow(e2@J) ), ncol=ncol(e1@J), nrow=nrow(e2@J), byrow=TRUE, sparse=TRUE )
        e1@J = Diagonal( x = e2@F ) %*% Jac1 + e1@F * e2@J
        e1@F = e1@F * e2@F
    } else if ( length(e1@F) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1@F) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), nrow(e1@J) ), ncol=ncol(e2@J), nrow=nrow(e1@J), byrow=TRUE, sparse=TRUE )
        e1@J = e2@F * e1@J + Diagonal( x = e1@F ) %*% Jac2
        e1@F = e1@F * e2@F
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }

    return( applyColoring(e1) )
})

#'@export
setMethod("/", c("FDiff","FDiff"), function (e1, e2) {
    vars = mergevars(e1@vars,e2@vars)
    e1   = expandJacDomain(e1,vars)
    e2   = expandJacDomain(e2,vars)
    
    # Perform different calculation depending on whether e1 or e2 is scalar.
    if ( length(e1@F) == length(e2@F) ) {
        e1@J = Diagonal( x = (1/e2@F)^2 ) %*% ( 
                      Diagonal( x = e2@F ) %*% e1@J -
                      Diagonal( x = e1@F ) %*% e2@J )
        e1@F = e1@F / e2@F
    } else if ( length(e1@F) == 1 & length(e2@F) > 1 ) {
        # Duplicate single row e1@J to get nrow(e2@F) rows.
        Jac1 = Matrix( rep( as.vector(e1@J), nrow(e2@J) ), ncol=ncol(e1@J), nrow=nrow(e2@J), byrow=TRUE, sparse=TRUE )
        e1@J = Diagonal( x = (1/e2@F)^2 ) %*% ( 
                      Diagonal( x = e2@F ) %*% Jac1 -
                      e1@F * e2@J )
        e1@F = e1@F / e2@F
    } else if ( length(e1@F) > 1 & length(e2@F) == 1 ) {
        # Duplicate single row e2@J to get length(e1@F) rows.
        Jac2 = Matrix( rep( as.vector(e2@J), nrow(e1@J) ), ncol=ncol(e2@J), nrow=nrow(e1@J), byrow=TRUE, sparse=TRUE )
        e1@J = (1/e2@F)^2 * ( e2@F * e1@J -
                    Diagonal( x = e1@F ) %*% Jac2 )
        e1@F = e1@F / e2@F
    } else {
        stop('Dimensions of ', deparse(substitute(e1)), ' and ', deparse(substitute(e2)), ' do not agree.', sep='')
    }
    return( applyColoring(e1) )
})

##
#
# Operators on matrix and FDiff.
#
##

#'@export
setMethod("%*%", c("matrix","FDiff"), function(x,y) {
    y@F = as.numeric(x %*% y@F)
    y@J = Matrix(x, sparse=TRUE) %*% y@J
    return( applyColoring(y) )
})

##
#
# Functions working separate elements of FDiff.
#
##

#'@export
# TODO: throw some error maybe if we have non-positive numbers?
setMethod("log", "FDiff", function(x) {
    # Order of defining J and F matters, as J is defined in terms of the original F (i.e. before taking the log).
    x@J = Diagonal( x = 1/(x@F) ) %*% x@J
    x@F = log(x@F)
    return( applyColoring(x) ) 
})

#'@export
# TODO: throw some error maybe if we have non-positive numbers?
setMethod("abs", "FDiff", function(x) {
    # Order of defining J and F matters, as J is defined in terms of the original F (i.e. before taking the log).
    x@J = sign(x@J)
    x@F = abs(x@F)
    return( applyColoring(x) ) 
})

##
#
# Functions working on all elements of FDiff at the same time.
#
##

#'@export
# TODO: is there a better way to find out the number of variables in a FDiff, rather
# than using length( x@F )?
setMethod("sum", "FDiff", function(x) {
    x@J = Matrix( 1, nrow=1, ncol=length(x@F), sparse=TRUE ) %*% x@J
    x@F = sum( x@F )
    return( x )
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
setMethod("appendJac", c("FDiff","sparseMatrix","list"), function(x,J,vs) { 
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
#' @method [ FDiff
#' @name extract 
setMethod(
  f= "[",
  signature="FDiff",
  definition=function(x,i){
    return(x@F[i])
  }
)

#' allows accessing the levels directly, this is convenient
#' within the code
#' @name extract 
setMethod(
  f= "dim",
  signature="FDiff",
  definition=function(x){
    return(dim(x@J))
  }
)

#' creates a variable to be tracked by the computation of the Jacobian
#' @param x a vector of current values
#' @param name a anme that uniquly identify this variable in the overall
#' @export
FDiff <- function(x, name, vars=NULL, coloring=FALSE) {
    x.length = length(x)
    if ( is.null(vars) ) {
        vars       <- list()
        vars[name] <- x.length
        F          <- x
        J          <- Matrix(diag(rep(1,x.length), nrow=x.length, ncol=x.length), sparse=TRUE)
    } else {
        # Run some checks on input.
        if ( is.null( vars[[name]] ) ) {
            stop('vars should contain at least name.')
        }
        if ( vars[[name]] != x.length ) {
            stop('length of name in vars is not the same as length of x.')
        }
    
        # Determine start and end index of name.
        total.cols <- sum( unlist(vars) )
        name.idx   <- (1:length(names(vars)))[ names( vars ) == name ]
        col.start  <- sum( unlist(vars[1:name.idx]) ) - vars[[name]] + 1
        
        # Set F and J.
        F              <- x
        J              <- Matrix( 0, nrow=x.length, ncol=total.cols, sparse=TRUE )
        J[ 1:x.length, col.start:(col.start + x.length - 1) ] <- Diagonal(x=rep(1,x.length))
    }
    
    # Return new FDiff object.
    return( new("FDiff",
        F    = F,
        J    = J,
        vars = vars,
        coloring=coloring ) )
}

#' function that extracts 

