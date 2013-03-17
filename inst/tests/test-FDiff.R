library('testthat')
library('mpeccable')
context("FDiff")

source("helper-finite-diff.R")

test_that("Test FDiff constructor.", {
    # Check that combining two FDiffs with the same name, but
    # with different elements result in an error.
    x1. <- FDiff( c(1,2,3), 'x' )
    x2. <- FDiff( 4, 'x' )
    expect_that( x1. + x2., throws_error() )
} )

##
#
# Binary operators("FDiff","numeric")
#
##

test_that("Test Operator(FDiff, numeric) with wrong input.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4), and a vector v of different dimension. Should result in error.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- c(2.0,3.0)
    x.    <- FDiff( x=xval, name='x' )
    
    expect_that( x. + v, throws_error() )       # Operator('+')
    expect_that( x. - v, throws_error() )       # Operator('-')
    expect_that( x. * v, throws_error() )       # Operator('*')
    expect_that( x. / v, throws_error() )       # Operator('/')
    expect_that( x. ^ v, throws_error() )       # Operator('^')
} )

test_that("Test Operator(FDiff, numeric) where FDiff is scalar and numeric is a scalar.", {
    # Check Operator( x., v ) with one value at x.=1.0 and a scalar v.
    xval  <- 2.8
    v     <- 0.3
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- x. + v
    expect_that( fx.@F, equals( xval + v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x + v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- x. - v
    expect_that( fx.@F, equals( xval - v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x - v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- x. * v
    expect_that( fx.@F, equals( xval * v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x * v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- x. / v
    expect_that( fx.@F, equals( xval / v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x / v }, .x = xval ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( x. + 2 ) / v ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( xval + 2 ) / v ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Operator('^')
    fx.   <- x. ^ v
    expect_that( fx.@F, equals( xval ^ v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ) ) )
} )

test_that("Test Operator(FDiff, numeric) where FDiff is scalar and numeric is a vector.", {
    # Check Operator( x., v ) with one value at x.=1.0 and a vector v.
    xval  <- 2.8
    v     <- c(-0.1,0.3,2)
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- x. + v
    expect_that( fx.@F, equals( xval + v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x + v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- x. - v
    expect_that( fx.@F, equals( xval - v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x - v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- x. * v
    expect_that( fx.@F, equals( xval * v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x * v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- x. / v
    expect_that( fx.@F, equals( xval / v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x / v }, .x = xval ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( x. + 2 ) / v ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( xval + 2 ) / v ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('^')
    fx.   <- x. ^ v
    expect_that( fx.@F, equals( xval ^ v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ) ) )
} )

test_that("Test Operator(FDiff, numeric) where FDiff is a vector and numeric is a scalar.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4) and scalar v.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- 0.3
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- x. + v
    expect_that( fx.@F, equals( xval + v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x + v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- x. - v
    expect_that( fx.@F, equals( xval - v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x - v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- x. * v
    expect_that( fx.@F, equals( xval * v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x * v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- x. / v
    expect_that( fx.@F, equals( xval / v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x / v }, .x = xval ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( x. + 2 ) / v ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( xval + 2 ) / v ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Operator('^')
    fx.   <- x. ^ v
    expect_that( fx.@F, equals( xval ^ v ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
} )

test_that("Test Operator(FDiff, numeric) where FDiff is a vector and numeric is a vector.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4) and vector v.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- c(2.0,3.0,4.0,5.0)
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- x. + v
    expect_that( fx.@F, equals( xval + v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x + v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- x. - v
    expect_that( fx.@F, equals( xval - v ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { x - v }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- x. * v
    expect_that( fx.@F, equals( xval * v ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x * v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Operator('/')
    fx.   <- x. / v
    expect_that( fx.@F, equals( xval / v ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x / v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( x. + 2 ) / v ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( xval + 2 ) / v ) - 4 ) )
    expect_that( fx.@J, equals( 
        as( Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Operator('^')
    fx.   <- x. ^ v
    expect_that( fx.@F, equals( xval ^ v ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ), tolerance=1e-7 ) )
} )

##
#
# Binary operators("numeric","FDiff")
#
##

test_that("Test Operator(numeric, FDiff) with wrong input.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4), and a vector v of different dimension. Should result in error.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- c(2.0,3.0)
    x.    <- FDiff( x=xval, name='x' )
    
    expect_that( v + x., throws_error() )       # Operator('+')
    expect_that( v - x., throws_error() )       # Operator('-')
    expect_that( v * x., throws_error() )       # Operator('*')
    expect_that( v / x., throws_error() )       # Operator('/')
    expect_that( v ^ x., throws_error() )       # Operator('^')
} )

test_that("Test Operator(numeric, FDiff) where numeric is scalar and FDiff is a scalar.", {
    # Check Operator( v, x. ) with a scalar v and one value at x.=1.0.
    xval  <- 2.8
    v     <- 0.3
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- v + x.
    expect_that( fx.@F, equals( v + xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v + x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- v - x.
    expect_that( fx.@F, equals( v - xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v - x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- v * x.
    expect_that( fx.@F, equals( v * xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v * x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- v / x.
    expect_that( fx.@F, equals( v / xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( v + 2 ) / x. ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( v + 2 ) / xval ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Operator('^')
    expect_that( v ^ x., throws_error() )   # Not yet implemented.
    #fx.   <- v ^ x.
    #expect_that( fx.@F, equals( v ^ xval ) )
    #expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v ^ x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
} )

test_that("Test Operator(numeric, FDiff) where numeric is scalar and FDiff is a vector.", {
    # Check Operator( x., v ) with one value at x.=1.0 and a vector v.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- 0.3
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- v + x.
    expect_that( fx.@F, equals( v + xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v + x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- v - x.
    expect_that( fx.@F, equals( v - xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v - x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- v * x.
    expect_that( fx.@F, equals( v * xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v * x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- v / x.
    expect_that( fx.@F, equals( v / xval ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( v + 2 ) / x. ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( v + 2 ) / xval ) - 4 ) )
    expect_that( fx.@J, equals( 
        as( Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ), "dgCMatrix" ), tolerance=1e-7 ) )
    
    # Operator('^')
    expect_that( v ^ x., throws_error() )   # Not yet implemented.
    #fx.   <- v ^ x.
    #expect_that( fx.@F, equals( v ^ xval ) )
    #expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v ^ x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
} )

test_that("Test Operator(numeric, FDiff) where numeric is a vector and FDiff is a scalar.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4) and scalar v.
    xval  <- 2.8
    v     <- c(-0.1,0.3,2)
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- v + x.
    expect_that( fx.@F, equals( v + xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v + x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- v - x.
    expect_that( fx.@F, equals( v - xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v - x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- v * x.
    expect_that( fx.@F, equals( v * xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v * x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.   <- v / x.
    expect_that( fx.@F, equals( v / xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( v + 2 ) / x. ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( v + 2 ) / xval ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Operator('^')
    expect_that( v ^ x., throws_error() )   # Not yet implemented.
    #fx.   <- v ^ x.
    #expect_that( fx.@F, equals( v ^ xval ) )
    #expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v ^ x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
} )

test_that("Test Operator(numeric, FDiff) where numeric is a vector and FDiff is a vector.", {
    # Check Operator( x., v ) at vector of values x.=c(1,2,3,4) and vector v.
    xval  <- c(1.0,2.0,3.0,4.0)
    v     <- c(2.0,3.0,4.0,5.0)
    x.    <- FDiff( x=xval, name='x' )
    
    # Operator('+')
    fx.   <- v + x.
    expect_that( fx.@F, equals( v + xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v + x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.   <- v - x.
    expect_that( fx.@F, equals( v - xval ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v - x }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.   <- v * x.
    expect_that( fx.@F, equals( v * xval ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v * x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Operator('/')
    fx.   <- v / x.
    expect_that( fx.@F, equals( v / xval ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ), tolerance=1e-7 ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( v + 2 ) / x. ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( v + 2 ) / xval ) - 4 ) )
    expect_that( fx.@J, equals( 
        as( Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Operator('^')
    expect_that( v ^ x., throws_error() )   # Not yet implemented.
    #fx.   <- v ^ x.
    #expect_that( fx.@F, equals( v ^ xval ) )
    #expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v ^ x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
} )

test_that("Test Operator(FDiff, FDiff) where FDiff is a scalar and FDiff is a scalar.", {
    # Check Operator( x., y. ) where x. and y. are vectors.
    # finite.diff returns a vector, so we have to explicitly convert this
    # vector to a row matrix (when you do not specify the number of rows 
    # and columns, Matrix by default creates a column matrix).
    xval     <- 2.3
    yval     <- 0.7
    x.       <- FDiff( x=xval, name='x' )
    y.       <- FDiff( x=yval, name='y' )
    
    # Operator('+')
    fx.      <- x. + y.
    testfunc <- function(x) { return( x[1:length(xval)] + x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), nrow=1, ncol=2, sparse=TRUE ) ) )
    
    # Operator('-')
    fx.      <- x. - y.
    testfunc <- function(x) { return( x[1:length(xval)] - x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), nrow=1, ncol=2, sparse=TRUE ) ) )
    
    # Operator('*')
    fx.      <- x. * y.
    testfunc <- function(x) { return( x[1:length(xval)] * x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), nrow=1, ncol=2, sparse=TRUE ) ) )
    
    # Operator('/')
    fx.      <- x. / y.
    testfunc <- function(x) { return( x[1:length(xval)] / x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), nrow=1, ncol=2, sparse=TRUE ), tolerance=1e-7 ) )
    
    # Combined operators.
    fx.      <- 2.5 * ( ( x. + 2 ) / y. ) - 4
    testfunc <- function(x) { 
        xvec <- x[1:length(xval)]
        yvec <- x[(length(xval)+1):length(x)]
        return( 2.5 * ( ( xvec + 2 ) / yvec ) - 4 )
    }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), nrow=1, ncol=2, sparse=TRUE ), tolerance=1e-7 ) )
} )

test_that("Test Operator(FDiff, FDiff) where FDiff is a scalar and FDiff is a vector.", {
    # Check Operator( x., y. ) where x. and y. are vectors.
    xval     <- 2.3
    yval     <- c(2.0,3.0,4.0,5.0)
    x.       <- FDiff( x=xval, name='x' )
    y.       <- FDiff( x=yval, name='y' )
    
    # Operator('+')
    fx.      <- x. + y.
    testfunc <- function(x) { return( x[1:length(xval)] + x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.      <- x. - y.
    testfunc <- function(x) { return( x[1:length(xval)] - x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.      <- x. * y.
    testfunc <- function(x) { return( x[1:length(xval)] * x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.      <- x. / y.
    testfunc <- function(x) { return( x[1:length(xval)] / x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.      <- 2.5 * ( ( x. + 2 ) / y. ) - 4
    testfunc <- function(x) { 
        xvec <- x[1:length(xval)]
        yvec <- x[(length(xval)+1):length(x)]
        return( 2.5 * ( ( xvec + 2 ) / yvec ) - 4 )
    }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
} )

test_that("Test Operator(FDiff, FDiff) where FDiff is a vector and FDiff is a scalar.", {
    # Check Operator( x., y. ) where x. and y. are vectors.
    xval     <- c(1.0,2.0,3.0,4.0)
    yval     <- 0.7
    x.       <- FDiff( x=xval, name='x' )
    y.       <- FDiff( x=yval, name='y' )
    
    # Operator('+')
    fx.      <- x. + y.
    testfunc <- function(x) { return( x[1:length(xval)] + x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.      <- x. - y.
    testfunc <- function(x) { return( x[1:length(xval)] - x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.      <- x. * y.
    testfunc <- function(x) { return( x[1:length(xval)] * x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.      <- x. / y.
    testfunc <- function(x) { return( x[1:length(xval)] / x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ), tolerance=1e-7 ) )
    
    # Combined operators.
    fx.      <- 2.5 * ( ( x. + 2 ) / y. ) - 4
    testfunc <- function(x) { 
        xvec <- x[1:length(xval)]
        yvec <- x[(length(xval)+1):length(x)]
        return( 2.5 * ( ( xvec + 2 ) / yvec ) - 4 )
    }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ), tolerance=1e-7 ) )
} )

test_that("Test Operator(FDiff, FDiff) where FDiff is a vector and FDiff is a vector.", {
    # Check Operator( x., y. ) where x. and y. are vectors.
    xval     <- c(1.0,2.0,3.0,4.0)
    yval     <- c(2.0,3.0,4.0,5.0)
    x.       <- FDiff( x=xval, name='x' )
    y.       <- FDiff( x=yval, name='y' )
    
    # Operator('+')
    fx.      <- x. + y.
    testfunc <- function(x) { return( x[1:length(xval)] + x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('-')
    fx.      <- x. - y.
    testfunc <- function(x) { return( x[1:length(xval)] - x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('*')
    fx.      <- x. * y.
    testfunc <- function(x) { return( x[1:length(xval)] * x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Operator('/')
    fx.      <- x. / y.
    testfunc <- function(x) { return( x[1:length(xval)] / x[(length(xval)+1):length(x)] ) }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.      <- 2.5 * ( ( x. + 2 ) / y. ) - 4
    testfunc <- function(x) { 
        xvec <- x[1:length(xval)]
        yvec <- x[(length(xval)+1):length(x)]
        return( 2.5 * ( ( xvec + 2 ) / yvec ) - 4 )
    }
    expect_that( fx.@F, equals( testfunc( c(xval, yval) ) ) )
    expect_that( fx.@J, equals( Matrix( finite.diff( func = testfunc, .x = c(xval,yval) ), sparse=TRUE ) ) )
} )

test_that("Test log(FDiff).", {
    # Check log( x. ) with one value at x.=1.0.
    xval  <- 1.0
    x.    <- FDiff( x=xval, name='x' )
    logx. <- log( x. )
    expect_that( logx.@F, equals( log( xval ) ) )
    expect_that( logx.@J, equals( Matrix( finite.diff( func = log, .x = xval ), sparse=TRUE ) ) )
    
    # Check log( x. ) with one value at x.=2.0.
    xval  <- 2.0
    x.    <- FDiff( x=xval, name='x' )
    logx. <- log( x. )
    expect_that( logx.@F, equals( log( xval ) ) )
    expect_that( logx.@J, equals( Matrix( finite.diff( func = log, .x = xval ), sparse=TRUE ) ) )
    
    # Check log( x. ) at vector of values x.=c(1,2,3,4).
    xval  <- c(1.0,2.0,3.0,4.0)
    x.    <- FDiff( x=xval, name='x' )
    logx. <- log( x. )
    expect_that( logx.@F, equals( log( xval ) ) )
    expect_that( logx.@J, equals( as( Matrix( finite.diff( func = log, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Check composite function.
    xval    <- c(1.0,2.0,3.0,4.0)
    x.      <- FDiff( x=xval, name='x' )
    funcx.  <- log( 3*( log(x.) + 2 ) )
    tmpfunc <- function(x) { log(3*(log(x) + 2)) }
    expect_that( funcx.@F, equals( tmpfunc( xval ) ) )
    expect_that( funcx.@J, equals( as( Matrix( finite.diff( func = tmpfunc, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )  
} )
