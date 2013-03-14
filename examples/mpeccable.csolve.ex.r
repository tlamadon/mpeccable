# The user writes the following:
# ------------------------------
require(mpeccable)
require(ggplot2)

# trying a simple fitting of a spline
# then we will try with constraint
# define the collocation 
# ----------------------
Nx = 50
xsupp  = seq(0,1,l=Nx)
ivals  = 1:3

# define a simple collocation
cc = expand.grid(a=xsupp,z=1:3)

# get a function to fit
cc$values = with(cc, a^z )
ggplot(cc,aes(x=a,y=values,color=factor(z))) + geom_line()

# create a functional and gets the function representation
V  = F_SplineInt1D(xsupp,ivals)
g. = param0(V,'g.')
g2. = param0(V,'g2.')

vars <- variables(g.)

cFunc <- function(theta,coloring=FALSE) {

  g. = extract(theta)
  # we want to minimize the distance between the spline and the values
  R = cc$values - V(cc$z,cc$z,g.)

}


res = mpeccable.solve( constfun, params, x0 ,options)






# Functions
V = SplineInt1D(xsupp,ivals)  # this can be stored on disk/memory

params .... # find a way

constfun <- function(x,coloring=FALSE){
# input arguments
  
  # get control variables and function representations
  rr = computeVarRanges(vars)
  a_ = getVar(x,'a',vars,coloring)
  g. = getVar(x,'a',vars,coloring)

  # write the constraints of the model
  R1 = V(cc$a,cc$z,g.)  - U(a_) - beta * EE %*% V(a_,cc$z,g.) 
  R2 = U1(a_) - beta * EE %*% V(a_,cc$z,g.,deriv=1)

  R = rBind2(R1,R2)

  # return constraints level and jacobian
  return(R)
}

res = mpeccable.solve( constfun, params, x0 ,options)