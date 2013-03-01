# example 1: cake eating in stationary environent
#================================================

# problem:
# V(a) = max u(a - a') + beta V(a')

# mpec:
# min_x sum_i^N (e_i)^2 

# s.t.
# 1) F(a_j,gamma) - u(a_j - a') - beta F(a',gamma) - e_j = 0
# 2) u'(a_k - a') - beta F_a(a',gamma) - e_k = 0

# F is the approximation to V, parameterized by M basis functions B and coefficients gamma
# V(x) ~ B(x) %*% gamma

# A is the length of the collocation of a
# index j=1,...,A indexes constraint number 1
# index k=A+1,...,2A indexes constraint number 2
# N=2A is the number of constraints

# x is the vector of choice variables: x = c(e,a',gamma)
# length of x is N = 2A + A + M
# length(e) = 2A
# length(a') = A
# length(gamma) = M

# the constraint function returns a vector 2A by 1
# jacobian function returns a matrix 2A by N


# make a collocation
nA <- 10
a <- 1:nA

# choose basis functions and coefficients
M <- 13	# will use 13 basis functions

N <- 2*nA + nA + M

# make a hash table to find errors and controls easier
library(data.table)
library(Matrix)
h <- data.table(coloc=rep(1:nA,3), var=rep(c("err.1","err.2","saving"),each=nA))
h <- rbind(h,data.table(coloc=1:M, var="gamma"))
#h[,idx := 1:nrow(h)]
stopifnot(nrow(h)==N)
setkey(h,var)

# make a choice vector
x <- runif(N)

# define FDiff


# return a value 'F' and a jacobian 'J'
# dim(F) = 2A,1
# dim(J) = 2A,(A+M)

# constaint function
con <- function(level,x){
	return( level - x[ h[c("err.1","err.2")][,idx] ] )		
}

# jacobian function
jac <- function(J,x){
	return( cbind(Matrix(diag(-1,2*nA),sparse=TRUE),J) )
}

# objective
objective <- function(x) {
	return( sum(  x[ h[c("err.1","err.2")][,idx] ]^2   ))
}

# gradient of objective
grad <- function(x) {
	return( 2 * x[ h[c("err.1","err.2")][,idx] ] )
}











 

