# The user writes the following:
# ------------------------------


# predefine objects

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