# The user writes the following:
# ------------------------------
require(mpeccable)

# define the collocation 
# ----------------------
Nx = 10
xsupp  = seq(0,1,l=Nx)
ivals  = 1:3

# define a simple collocation
cc = expand.grid(a=xsupp,z=1:3)

# create a functional and gets the function representation
V  = F_SplineInt1D(xsupp,ivals)
g. = param0(V,'g.')

# define a control variable
a_ = FDiff(cc$a,'a')

# step1, get the structure of the Jacobian
options(mpeccable.coloring=TRUE)
R = V(cc$a,cc$z,g.)
R = V(a_,  cc$z,g.)




#options(mpeccable.coloring=FALSE)
# create a low of motion matrix, where 1->2 , 2->3, 3->1
z1 = (cc$z %% 3)+ 1 

# create a simple Euler Equation
R2 = a_ + V(cc$a,cc$z,g.) + 0.9 * V(a_,z1,g.)
image(R2@J)


# collocation

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