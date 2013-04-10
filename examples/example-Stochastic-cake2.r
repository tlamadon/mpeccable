

rm(list=ls(all=T))
# make 5 period model out of this.
# note that z has no influence on payoff
# is just a counter variable

require(mpeccable)
require(ggplot2)
require(ipoptr)

Na =3 
asupp = seq(0.1, 20, l = Na)
zvals = 1:2

r    = 0.02
beta = 0.95 

cc = expand.grid(a = asupp, z = zvals)
cc = data.table(expand.grid(a = asupp, z = zvals))
cc[,z1 := z + 1]

V = F_SplineTime1D(asupp, zvals, degree=3,nbasis.funs=4)
# V = F_SplineInt1D(asupp, zvals)
g. = param0(V, "g.", 1)
c_ = FDiff(cc[,a/2], "c_")	# choose consumption


maxtime <- cc[,max(z1)]	# last period

#need to be able to tell mpec that there is a so-called 'last period' where the functional form is known, and no approximation is done.

# storing the list of variables 
mpec = mpec.new()
mpec = mpec.addVar(mpec,g.)
mpec = mpec.addVar(mpec,c_,lb=0.05,ub=cc[,a])	
# s = (a - c) * (1+r)
# positive savings: s>0 implies
# c < a


# initial value
x0 = mpec.getVarsAsVector(mpec)

# creating the function that computes the list of constraints
cFunc <- function(mpec) {

  # get the parameters
  g. = mpec.getVar(mpec,'g.')
  c_ = mpec.getVar(mpec,'c_')

  # compute the utility and its derivative
  U    = log( c_ )
  Uc_  = 1 / U

  # compute implied savings
  #   s_ = (1+r) * (cc$a - c_)

  # we append the level equation
  R =  V(cc$a,cc$z,g.)  - U - beta*V((1+r) * (cc$a - c_),cc$z1,g.)
  mpec = mpec.addAbsConstraint(mpec,R,'BE.level')

  # and the first order condition
  # aka Euler Equation
  R2 =  Uc_ - beta*(1+r)*V((1+r) * (cc$a - c_),cc$z1,g.,deriv=1)
  mpec = mpec.addAbsConstraint(mpec,R2,'BE.foc')

  # we constrain consumption to be positive
  #   R3 = (1 + r) * cc$a  + cc$z/3  - a_
  #   mpec = mpec.addInequalityConstraint(mpec,R3,'c.pos',lb=0)

  # we add the envelope condition
  #   R4 = V(cc$a,cc$z,g.,deriv=1) - (1+r)*beta*V(s,z1,g.,deriv=1)
  R4 = V(cc$a,cc$z,g.,deriv=1) - Uc_
  mpec = mpec.addAbsConstraint(mpec,R4,'BE.env')

  # and finally a concavity restriction
  #   R5 = V(cc$a,cc$z,g.,deriv=2)
  #   mpec = mpec.addInequalityConstraint(mpec,R5,'V.shape',ub=0)


  # return the mpec object
  return(mpec)
}

opts <- list("print_level"=5,
             "tol"=1.0e-8,
			 #              "derivative_test"="first-order",
             "max_iter"=40)

# call the optimizer
res = mpec.solve(cFunc=cFunc, x0=x0, mpec = mpec , opts = opts)

# recover parameters

g.         = res$solution[['g.']]
c_         = res$solution[['c_']]
cc$c_opt   = c_@F
cc$saving  = (1 + r) * (cc$a - cc$c_opt)
cc$V       = V(cc$a,cc$z,g.)@F

cc.l = melt(cc,id.vars=c('a','z'))

ggplot(subset(cc.l,variable != "z1"),aes(x=a,y=value,color=factor(z))) + 
  geom_line() + facet_grid(variable~.,scales="free_y")



