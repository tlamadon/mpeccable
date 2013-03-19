# This is very simple example that fits 3 splines
# on a very small grid. It is supposed to test
# the piping

require(mpeccable)
require(ggplot2)
require(ipoptr)

# trying a simple fitting of a spline
# then we will try with constraint
# define the collocation 
# ----------------------
Nx = 15
xsupp  = seq(0,1,l=Nx)
ivals  = 1

# define a simple collocation with more points
cc = expand.grid(a=seq(0,1,l=Nx),z=ivals)

# get a function to fit
cc$values = with(cc, a^z ) + rnorm(nrow(cc),sd=0.01)

# create a functional and gets the function representation
V  = F_SplineInt1D(xsupp,ivals)
g. = param0(V,'g.',0.1)

# storing the list of variables 
vars <- mpec.vars.desc(list(g.))

# initial value
x0 = g.@F

# creating the function that computes the list of constraints
cFunc <- function(params) {

  # get the parameters
  g. = params[['g.']]

  # we want to minimize the distance between the spline and the values
  R = cc$values - V(cc$a,cc$z,g.)

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
g.  = res$solution[['g.']]
cc$values_mpec = V(cc$a,cc$z,g.)@F


# fitting the spline in a different way
rbest = c()
for (i in 1:3) {
  fit <- smooth.spline(cc$a[cc$z==i],cc$values[cc$z==i])
  rbest = c(rbest,fit$fit$coef)
}
g.@F = rbest
cc$values_spline = V(cc$a,cc$z,g.)@F

ggplot(cc,aes(x=a,y=values,color=factor(z))) + 
  geom_point() + geom_line(aes(y=values_mpec,linetype='mpec')) +
  geom_line(aes(y=values_spline,linetype='spline')) +
  theme_bw()



