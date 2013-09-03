require(mpeccable)
require(ggplot2)
require(ipoptr)

Na =15 
asupp = seq(1, 20, l = Na)
zvals = 1:3

rho  = 0.2
r    = 0.02
beta = 0.95 

cc = expand.grid(a = seq(1,19,l=2*Na), z = zvals)

V = F_SplineInt1D(asupp, zvals)
g. = param0(V, "g.", 1)
a_ = FDiff(cc$a/2, "a_")
(1 + r) * cc$a  + cc$z/3 - a_
z1 = (cc$z%%3) + 1

# storing the list of variables 
mpec = mpec.new()
mpec = mpec.addVar(mpec,g.)
mpec = mpec.addVar(mpec,a_,lb=0,ub=cc$a)

# initial value
x0 = mpec.getVarsAsVector(mpec)

# creating the function that computes the list of constraints
cFunc <- function(mpec) {

  # get the parameters
  g. = mpec.getVar(mpec,'g.')
  a_ = mpec.getVar(mpec,'a_')

  # compute the utility and its derivative
  U   =     log( (1 + r) * cc$a  + cc$z/3 - a_)
  Ua_ =    - 1/( (1 + r) * cc$a  + cc$z/3 - a_)
  Ua  =  (1+r)/( (1 + r) * cc$a  + cc$z/3 - a_)

  # we append the level equation
  R =  V(cc$a,cc$z,g.)  - U - beta*V(a_,z1,g.)
  mpec = mpec.addAbsConstraint(mpec,R,'BE.level')

  # and the first order condition
  R2 =  Ua_ + beta*V(a_,z1,g.,deriv=1)
  mpec = mpec.addAbsConstraint(mpec,R2,'BE.foc')

  # we constrain consumption to be positive
  R3 = (1 + r) * cc$a  + cc$z/3  - a_
  mpec = mpec.addInequalityConstraint(mpec,R3,'c.pos',lb=0)

  # we add the envelope condition
  R4 = V(cc$a,cc$z,g.,deriv=1) - Ua
  mpec = mpec.addAbsConstraint(mpec,R4,'BE.env')

  # and finally a concavity restriction
  R5 = V(cc$a,cc$z,g.,deriv=2)
  mpec = mpec.addInequalityConstraint(mpec,R5,'V.shape',ub=0)


  # return the mpec object
  return(mpec)
}

opts <- list("print_level"=5,
             "tol"=1.0e-8,
			 "derivative_test"="first-order",
             "max_iter"=40)

# call the optimizer
res = mpec.solve(cFunc=cFunc, x0=x0, mpec = mpec , opts = opts)

# recover parameters

g.         = res$solution[['g.']]
a_         = res$solution[['a_']]
cc$a_opt   = a_@F
cc$consumption = (1 + r) * cc$a  + cc$z/3  - cc$a_opt
cc$V       = V(cc$a,cc$z,g.)@F

cc.l = melt(cc,id.vars=c('a','z'))

ggplot(cc.l,aes(x=a,y=value,color=factor(z))) + 
  geom_line() + facet_wrap(~variable,scales="free_y")



