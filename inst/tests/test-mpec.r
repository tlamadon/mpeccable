library('testthat')
library('mpeccable')
library('ipopt')
context("FDiff")

test_that("check if the mpec object behave correctly", {

  mpec = mpec.new()
  g. = FDiff(1:10,'g.')
  mpec = mpec.addVar(mpec,g.,lb=0)

  expect_true( 'g.' %in% names(mpec$vars) ,label='error in first adding variable to mpec object')
  expect_true( mpec$vars[['g.']] == 10    ,label='error in length of first added variable')

  u. = FDiff(-1:-20,'u.')
  mpec = mpec.addVar(mpec,u.,ub=0)

  expect_true( 'g.' %in% names(mpec$vars))
  expect_true( mpec$vars[['g.']] == 10)
  expect_true( 'u.' %in% names(mpec$vars))
  expect_true( mpec$vars[['u.']] == 20)

  expect_true(length(mpec.getVarsAsVector(mpec)) == 30)

  xx   = runif(30)
  mpec = mpec.setVarsFromVector(mpec,xx)
  Iu   = mpec.getVarRange(mpec,'u.')
  u.   = mpec.getVar(mpec,'u.')
  expect_true( all(u.@F == xx[Iu]) ,label='error in using setVarsFromVector or in getVarRange')

} )

