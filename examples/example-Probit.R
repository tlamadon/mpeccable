# This example estimates the parameters of a probit model.

# Load libraries.
library('mpeccable')
library('nloptr')

# Define parameters.
beta.true <- c( 1, 2, -0.5 )
num.obs   <- 1000

# Set random seed.
set.seed( 3141 )

# Generate data.
X        <- cbind( rep( 1, num.obs ), matrix( rnorm( num.obs * (length(beta.true)-1) ), nrow=num.obs, ncol=length(beta.true)-1 ) )
eps      <- rnorm( num.obs )
y_latent <- as.vector( X %*% beta.true + eps )
y        <- 1*( y_latent > 0 )

applyfun <- function(x, f, df) {
    #x@J = Matrix( diag( df( x@F ), nrow=length(x@F), ncol=length(x@F) ), sparse=TRUE ) %*% x@J
    x@J = Diagonal( x = df( x@F ) ) %*% x@J
    x@F = f( x@F )
    return( x )
}

#pnorm2 <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
pnorm2 <- function(q) {
    stopifnot( is( q, "FDiff" ) )
    return( applyfun( x=q, f=pnorm, df=dnorm ) )
}

# Define optimization options.
opts <- list( "algorithm"         = "NLOPT_LD_LBFGS",
              "xtol_rel"          = 1.0e-4,
              "check_derivatives" = TRUE )

# Define initial values.
params0 <- beta.true

# Define objective function.
eval_f_list <- function( params, y, X ) {
    stopifnot( length(params) == ncol(X) )
    stopifnot( length(y) == nrow(X) )
    num.obs <- length( y )
    
    # Define parameters.
    beta.  <- FDiff( x=params, name='beta' )

    # Calculate log-likelihood.
    Xbeta.  <- X %*% beta.
    loglik. <- sum( y * log( pnorm2( Xbeta. ) ) + (1-y) * log( pnorm2( -Xbeta. ) ) )
    
    # Take minus log-likelihood.
    objective. <- -loglik. / num.obs
    
    return( list( "objective" = objective.@F,
                  "gradient"  = as.vector( objective.@J ) ) )
}

# Find parameters that maximize the likelihood.
time.start <- Sys.time();
res <- nloptr( x0     = params0, 
               lb     = rep( -Inf, length(beta.true) ),
               ub     = rep(  Inf, length(beta.true) ),
               eval_f = eval_f_list,
               opts   = opts,
               y      = y,
               X      = X )
time.stop <- Sys.time()
print(time.stop - time.start)
print( res )
print(glm( y ~ -1 + X, family=binomial(link=probit) ))
