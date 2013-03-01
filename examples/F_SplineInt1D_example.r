# define the support of the function
Nx = 10
xsupp  = seq(0,1,l=Nx)
ivals  = 1:3

# define a simple collocation
N=10
cc = expand.grid(a=xsupp,z=1:3)

# input arguments
g. = FDiff(rep(1,length(ivals) * Nx),'g')
a_ = FDiff(cc$a,'a')

V = F_SplineInt1D(xsupp,ivals)

options(mpeccable.coloring=TRUE)
R = V(cc$a,cc$z,g.)
R = V(a_,  cc$z,g.)

#options(mpeccable.coloring=FALSE)
# create a low of motion matrix, where 1->2 , 2->3, 3->1
z1 = (cc$z %% 3)+ 1 

# create a simple Euler Equation
R2 = a_ + V(cc$a,cc$z,g.) + 0.9 * V(a_,z1,g.)
image(R2@J)


