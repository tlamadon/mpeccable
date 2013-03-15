# This is very simple example that fits 3 splines
# on a very small grid. It is supposed to test
# the piping

require(mpeccable)
require(ggplot2)
require(ipoptr)

# trying a simple linear regression
# ----------------------
Nx = 20

# define a simple collocation with more points
cc = data.frame(x=seq(0,1,l=5*Nx),c=1)

# get a function to fit
cc$y = with(cc, x^2) + rnorm(nrow(cc),sd=0.1)

# define 2 parameters
F     = F_Spline(seq(0,1,l=Nx))
g. = param0(F,'g.',0)

# storing the list of variables 
vars <- mpec.vars.desc(list(g.))

# initial value
x0 = mpec.vars.init(vars,0)

# creating the function that computes the list of constraints
cFunc <- function(params) {

  # get the parameters
  g. = params[['g.']]

  # we want to minimize the distance between the spline and the values
  R =  F(cc$x,g.) - cc$y

  # return the list of constraints, sepecifying each type
  # for now it's only MSE constraints
  return(list(   C.ABS= R))
}

objFunc <- function(params) {
  # get the parameters
  g. = params[['g.']]

  # we want to regularize the spline so we add the sum of derivatives squared
  OBJ = 0.0001*sum(F(cc$x,g.,deriv=2)^2)

  return(list(OBJ=OBJ))
}


opts <- list("print_level"=5,
             "tol"=1.0e-8,
             "max_iter"=100)

# call the optimizer
res = mpeccable.asolve(cFunc=cFunc, x0=x0, vars = vars , opts = opts)

# extract results
g.         = res$solution[['g.']]
cc$y_hat = F(cc$x,g.)@F

# compute using actual spline method
# fit <- smooth.spline(cc$x,cc$y)
# g.@F = fit$fit$coef
# cc$y_hat2 = F(cc$x,g.)@F

# inverting the matrix directly
M = F(cc$x,g.)@J
g.@F = as.numeric(solve(t(M)%*%M ,t(M) %*% cc$y ))
cc$y_hat3 = F(cc$x,g.)@F

ggplot(cc,aes(x=x,y=y)) + 
  geom_point() + geom_line(aes(y=y_hat,color='mpec')) + 
  #geom_point() + geom_line(aes(y=y_hat2,color='smoothing spline')) + 
  geom_point() + geom_line(aes(y=y_hat3,color='matrix inversion')) + 
  theme_bw()





