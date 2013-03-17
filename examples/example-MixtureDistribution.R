# This example estimates the parameters of a mixture of two
# normal distributions.

# Load libraries.
library('mpeccable')
library('nloptr')

# Define parameters.
mu1.true    <-  2.0
sigma1.true <-  2.0
mu2.true    <- -1.0
sigma2.true <-  1.0
pi1.true    <-  0.3
num.obs     <- 1000

# Set random seed.
set.seed( 3141 )

# Generate data.
y1      <- rnorm( num.obs, mean=mu1.true, sd=sigma1.true )
y2      <- rnorm( num.obs, mean=mu2.true, sd=sigma2.true )
select1 <- 1*( runif( num.obs ) < pi1.true )
y       <- ifelse( select1, y1, y2 )

applyfun <- function(x, f, df) {
    x@J = Diagonal( x = df( x@F ) ) %*% x@J
    x@F = f( x@F )
    return( x )
}

dnorm2 <- function(x) {
    stopifnot( is( x, "FDiff" ) )
    return( applyfun( x=x, f=dnorm, df=function(x) { -x*dnorm(x) } ) )
}

# Define optimization options.
opts <- list( "algorithm"         = "NLOPT_LD_LBFGS",
              "xtol_rel"          = 1.0e-4,
              "check_derivatives" = TRUE )

# Define initial values.
params0 <- c( mu1.true, sigma1.true, mu2.true, sigma2.true, pi1.true )

# Define objective function.
eval_f_list <- function( params, y ) {
    stopifnot( length(params) == 5 )
    num.obs <- length( y )
    
    mean1.  <- FDiff( x=params[1], name='mean1' )
    sigma1. <- FDiff( x=params[2], name='sigma1' )
    mean2.  <- FDiff( x=params[3], name='mean2' )
    sigma2. <- FDiff( x=params[4], name='sigma2' )
    pi1.    <- FDiff( x=params[5], name='pi1' )

    vars    <- mpeccable:::mergevars(mean1.@vars, sigma1.@vars)
    vars    <- mpeccable:::mergevars(vars, mean2.@vars)
    vars    <- mpeccable:::mergevars(vars, sigma2.@vars)
    vars    <- mpeccable:::mergevars(vars, pi1.@vars)
    
    mean1.  <- mpeccable:::expandJacDomain(mean1., vars)
    sigma1. <- mpeccable:::expandJacDomain(sigma1., vars)
    mean2.  <- mpeccable:::expandJacDomain(mean2., vars)
    sigma2. <- mpeccable:::expandJacDomain(sigma2., vars)
    pi1.    <- mpeccable:::expandJacDomain(pi1., vars)
    
    loglik. <- sum( log( pi1.  * (1/(sqrt(2*pi)*sigma1.)) * dnorm2( ( y - mean1. )/sigma1. ) + 
                      (1-pi1.) * (1/(sqrt(2*pi)*sigma2.)) * dnorm2( ( y - mean2. )/sigma2. ) ) )
    
    # Take minus log-likelihood.
    objective. <- -loglik. / num.obs
    
    return( list( "objective" = objective.@F,
                  "gradient"  = as.vector( objective.@J ) ) )
}

# Find parameters that maximize the likelihood.
res <- nloptr( x0     = params0, 
               lb     = c( -Inf, 0.0001, -Inf, 0.0001, 0.01 ),
               ub     = c(  Inf,    Inf,  Inf,    Inf, 0.99 ),
               eval_f = eval_f_list,
               opts   = opts,
               y      = y )
print( res )


ll <- function( params, y ) {
    stopifnot( length(params) == 5 )
    num.obs <- length( y )
    
    mean1.  <- params[1]
    sigma1. <- params[2]
    mean2.  <- params[3]
    sigma2. <- params[4]
    pi1.    <- params[5]

    #loglik. <- sum( log( pi1. * dnorm( ( y - mean1. )/sigma1. ) + (1-pi1.) * dnorm( ( y - mean2. )/sigma2. ) ) )
    loglik. <- sum( log( pi1. * dnorm( y, mean=mean1., sd=sigma1. ) + (1-pi1.) * dnorm( y, mean=mean2., sd=sigma2. ) ) )
    
    # Take minus log-likelihood.
    objective. <- -loglik. / num.obs
    
    return( objective. )
}
optim( 
    par=params0, fn=ll, 
    lower=c( -Inf, 0.0001, -Inf, 0.0001, 0.01 ), 
    upper=c(  Inf,    Inf,  Inf,    Inf, 0.99 ), 
    method="L-BFGS-B", 
    y=y )
