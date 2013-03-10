# Copyright (C) 2011-2013 Jelmer Ypma. All Rights Reserved.
# This code is published under the L-GPL.
#
# File:   helper-finite-diff.R
# Author: Jelmer Ypma
# Date:   24 July 2011
#
# Approximate the gradient of a function using finite differences.
#
# Input: 
#		func       : Calculate finite difference approximation for the gradient of this function.
#		.x         : Evaluate at this point.
#		indices    : Calculate finite difference approcximation only for .x-indices in this vector (optional).
#		stepSize   : Relative size of the step between x and x+h (optional).
#       method     : Method to be used to calculate finite difference. Values can be forward (default), backward or central.
#       returnList : Boolean. If FALSE (default) only the gradient is returned. If TRUE, a list with "objective" and "gradient" is returned.
#       fx         : Value of function evaluated at .x. This value can be supplied if it is known, to save one function evaluation.
#       ...        : arguments that are passed to the user-defined function (func).
#
# Output: Matrix with finite difference approximations, or list with the value of the function ("objective") and its gradient ("gradient").
#
# Example:
#    finite.diff( func=function(x) { return( log( x ) ) }, .x=c(1,2,3), method="central" )
#
# Notes:
#    http://en.wikipedia.org/wiki/Numerical_differentiation
finite.diff <- function( 
    func, 
    .x, 
    indices    = 1:length(.x), 
    stepSize   = sqrt( .Machine$double.eps ), 
    method     = "forward", 
    returnList = FALSE,
    fx         = NULL,
    ... )
{
    
    # Choose stepSize depending on whether |.x| is bigger or smaller than 1.
    stepSizeVec <- pmax( abs(.x), 1 ) * stepSize 
    
    # Calculate function value at x if it has not been supplied.
    if ( is.null( fx ) ) {
        fx <- func( .x, ... )
    }
    
    # Locally define function to approximate gradient.
    if ( method == "forward" ) {
        approx.gradf.index <- function(i, .x, func, fx, stepSizeVec, ...) {
            x_prime <- .x
            x_prime[i] <- .x[i] + stepSizeVec[i]
            stepSizeVec[i] <- x_prime[i] - .x[i]
            fx_prime <- func( x_prime, ... )
            return( ( fx_prime - fx )/stepSizeVec[i] )
        }
    }
    else if ( method == "backward" ) {
        approx.gradf.index <- function(i, .x, func, fx, stepSizeVec, ...) {
            x_prime <- .x
            x_prime[i] <- .x[i] - stepSizeVec[i]
            stepSizeVec[i] <- .x[i] - x_prime[i]
            fx_prime <- func( x_prime, ... )
            return( ( fx - fx_prime )/stepSizeVec[i] )
        }
    }
    else if ( method == "central" ) {
        approx.gradf.index <- function(i, .x, func, fx, stepSizeVec, ...) {
            x_prime_f <- .x
            x_prime_b <- .x
            x_prime_f[i] <- .x[i] + 0.5 * stepSizeVec[i]    # forward
            x_prime_b[i] <- .x[i] - 0.5 * stepSizeVec[i]    # backward
            fx_prime_f <- func( x_prime_f, ... )
            fx_prime_b <- func( x_prime_b, ... )
            return( ( fx_prime_f - fx_prime_b )/stepSizeVec[i] )
        }    
    }
    else {
        stop( "method should be 'forward', 'backward', or 'central'")
    }
    
    # Do approximation.
    grad_fx <- sapply(indices, approx.gradf.index, .x=.x, func=func, fx=fx, stepSizeVec=stepSizeVec, ... )
    
    if ( returnList ) {
        # Return objective and gradient.
        return( list( "objective" = fx, "gradient" = grad_fx ) )
    } else {
        # Return gradient.
        return( grad_fx )
    }
}
