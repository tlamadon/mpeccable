library('testthat')
library('mpeccable')
library('ipoptr')
context("FDiff")

test_that("Check that F_splineInt1D works", {


  Na = 10
  asupp = seq(0, 1, l = Na)
  zvals = 1:3

  cc = expand.grid(a = seq(0.001,1,l=2*Na), z = zvals)

  V = F_SplineInt1D(asupp, zvals)
  g. = param0(V, "g.", 1)
  a_ = FDiff(cc$a/2, "a_")

  V(a_,cc$z,g.)




    expect_true(all(unlist(M1) == unlist(M2)))
} )

