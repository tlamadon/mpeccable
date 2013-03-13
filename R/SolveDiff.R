# Simple function representation
# ==============================

# we want to represent a function in a finite subspace
# and evaluate it at the collocation. 
# a function representation is an R function that takes
# a collocation and a parameter vector and returns 
# the correspong FDiff object.

# let's first look at a one dimensional spline


#' SolveDiff object
#' takes an FDiff object and runs the opimization on ipoptr
#' @param F a functional representation of the problem
#' @param x0 a starting value for the optimization
#' @param ub a vector of upper bounds
#` @param lb a vector of lower bounds
setClass("SolveDiff", 
    representation(a="numeric") )

# methods: 
# setConstraintBounds
# setChoicevarBounds
# setStopping
# getConstraintBounds
# getChoicevarBounds
# getStopping

