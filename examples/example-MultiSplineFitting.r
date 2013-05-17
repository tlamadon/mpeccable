require(mpeccable)
require(ggplot2)
require(ipoptr)


# define the support of the function
Nx = 10
xsupp  = seq(0,1,l=Nx)
ivals  = 1:5

# simulate a sample with 5 times more points
N = 5 * Nx
cc   = expand.grid(x=seq(0,1,l=N),i=ivals)
cc$y = with(cc, x^(1/i)) + rnorm(nrow(cc),sd=0.05)

# create the functional representation
V  = F_SplineInt1D(xsupp,ivals)
v. = param0(V,'v.',0)

# create a constrain optimization problem
mpec = mpec.new()
mpec = mpec.addVar(mpec,v.)

# initial value
x0 = mpec.getVarsAsVector(mpec)

# creating the function that computes the list of constraints
cFunc <- function(mpec) {

  # get the parameters
  v. = mpec.getVar(mpec,'v.')

  # we want to minimize the distance between the spline and the values
  R    = V(cc$x,cc$i,v.) - cc$y
  mpec = mpec.addAbsConstraint(mpec,R,'spline.error')
  R2   = V(cc$x,cc$i,v.,deriv=2)
  mpec = mpec.addInequalityConstraint(mpec,R2,'spline.shape',ub=0)

  # return the mpec object
  return(mpec)
}

# simple test
R = V(cc$x,cc$i,v.)

opts <- list("print_level"=5,
             "tol"=1.0e-8,
			 "derivative_test"="first-order",
             "max_iter"=100)

# call the optimizer
res = mpec.solve(cFunc=cFunc, x0=x0, mpec = mpec , opts = opts)

# extract results
v.         = res$solution[['v.']]
cc$y_hat   = V(cc$x,cc$i,v.)@F

ggplot(cc,aes(y=y,x=x,color=factor(i))) + geom_point() + geom_line(aes(y=y_hat))


