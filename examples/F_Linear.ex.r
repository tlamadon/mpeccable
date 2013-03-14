# This is very simple example that fits 3 splines
# on a very small grid. It is supposed to test
# the piping

require(mpeccable)
require(ggplot2)
require(ipoptr)

# trying a simple linear regression
# ----------------------
Nx = 100

# define a simple collocation with more points
cc = data.frame(x=seq(0,1,l=Nx),c=1)
cc$x2 = cc$x^2

# get a function to fit
cc$y = with(cc, x + x2 + 1 ) + rnorm(nrow(cc),sd=0.1)

# define 2 parameters
MX    = F_Linear(cc[,c('x','x2','c')])
beta. = param0(MX,'beta.',0)

# storing the list of variables 
vars <- mpec.vars.desc(list(beta.))

# initial value
x0 = mpec.vars.init(vars,0)

# creating the function that computes the list of constraints
cFunc <- function(params) {

  # get the parameters
  beta. = params[['beta.']]

  # we want to minimize the distance between the spline and the values
  R = cc$y - MX(beta.)

  # return the list of constraints, sepecifying each type
  # for now it's only MSE constraints
  return(list(C.MSE=R))
}

opts <- list("print_level"=0,
             "file_print_level"=12,
             "output_file"="banana.out",
             "tol"=1.0e-8)

# call the optimizer
res = mpeccable.csolve( cFunc, x0, vars , opts)

# extract results
beta.         = res$solution[['beta.']]
cc$y_hat = MX(beta.)@F

ggplot(cc,aes(x=x,y=y)) + 
  geom_point() + geom_line(aes(y=y_hat,linetype='mpec')) + 
  theme_bw()



