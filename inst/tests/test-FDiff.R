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
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ) ) )
    
    # Operator('^')
    fx.   <- x. ^ v
    expect_that( fx.@F, equals( xval ^ v ) )
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
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
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( x + 2 ) / v ) - 4 }, .x = xval ), sparse=TRUE ) ) )
    
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
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { x ^ v }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
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
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
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
        as( Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
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
    expect_that( fx.@J, equals( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ) ) )
    
    # Combined operators.
    fx.   <- 2.5 * ( ( v + 2 ) / x. ) - 4
    expect_that( fx.@F, equals( 2.5 * ( ( v + 2 ) / xval ) - 4 ) )
    expect_that( fx.@J, equals( 
        Matrix( finite.diff( func = function(x) { 2.5 * ( ( v + 2 ) / x ) - 4 }, .x = xval ), sparse=TRUE ) ) )
    
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
    expect_that( fx.@J, equals( as( Matrix( finite.diff( func = function(x) { v / x }, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
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

test_that("Test log(FDiff).", {
    # Check log( x. ) with one value at x.=1.0.
    xval  <- 1.0
    x.    <- FDiff( x=xval, name='x' )
    logx. <- log( x. )
    expect_that( logx.@F, equals( log( xval ) ) )
    expect_that( logx.@J, equals( as( Matrix( finite.diff( func = log, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
    # Check log( x. ) with one value at x.=2.0.
    xval  <- 2.0
    x.    <- FDiff( x=xval, name='x' )
    logx. <- log( x. )
    expect_that( logx.@F, equals( log( xval ) ) )
    expect_that( logx.@J, equals( as( Matrix( finite.diff( func = log, .x = xval ), sparse=TRUE ), "dgCMatrix" ) ) )
    
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
