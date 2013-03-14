# The user writes the following:
# ------------------------------
require(mpeccable)
require(ggplot2)
require(ipoptr)

# trying a simple fitting of a spline
# then we will try with constraint
# define the collocation 
# ----------------------
Nx = 10
xsupp  = seq(0,1,l=Nx)
ivals  = 1:3

# define a simple collocation
cc = expand.grid(a=xsupp,z=ivals)

# get a function to fit
cc$values = with(cc, a^z )
ggplot(cc,aes(x=a,y=values,color=factor(z))) + geom_line()

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

# small test
ps  = mpec.vars.collate(x0,vars,coloring=TRUE)
res = cFunc(ps)

res = mpeccable.solve( cFunc, x0, vars )
