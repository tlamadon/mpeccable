library('testthat')
library('mpeccable')
library('ipoptr')
context("FDiff")

test_that("Check that ipoptr.sparse is correct", {


    M = Matrix(  pmax(runif(30^2)-0.8,0) ,30,30,sparse=TRUE)
    M1 = make.sparse(M)
    M2 = ipoptr.sparse(M)

    expect_true(all(unlist(M1) == unlist(M2)))
} )

