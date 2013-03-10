library('testthat')
library('mpeccable')
context("FDiff")

source("helper-finite-diff.R")

test_that("test log(FDiff)", {
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
