# This example estimates the parameters of a normal linear model
#     y_i = X_i'beta + epsilon_i
#     epsilon_i ~ i.i.d. N(0,sigma^2)
# by an iterative maximum likelihood procedure. The same estimate
# can of course be obtained in closed form, but this is to for
# illustration purposes only.

# Load libraries.
library('mpeccable')
library('nloptr')

# Define parameters.
beta.true  <- c( 1, 2, -0.5 )
sigma.true <- 2
num.obs    <- 1000

# Set random seed.
set.seed( 3141 )

# Generate data.
X   <- cbind( rep( 1, num.obs ), matrix( rnorm( num.obs * (length(beta.true)-1) ), nrow=num.obs, ncol=length(beta.true)-1 ) )
eps <- sigma.true * rnorm( num.obs )
y   <- as.vector( X %*% beta.true + eps )

# Define optimization options.
opts <- list( "algorithm"         = "NLOPT_LD_LBFGS",
              "xtol_rel"          = 1.0e-4,
              "check_derivatives" = TRUE )

# Define initial values.
params0 <- c( beta.true, sigma.true )

# Define objective function.
eval_f_list <- function( params, y, X ) {
    stopifnot( length(params) == ncol(X) + 1 )
    stopifnot( length(y) == nrow(X) )
    num.obs <- length( y )
    
    beta.  <- FDiff( x=params[1:ncol(X)], name='beta' )
    sigma. <- FDiff( x=params[length(params)], name='sigma' )

    # Note: the order of evaluating the log-likelihood matters,
    # because the order in which the parameters occur, changes
    # the order in which the occur in columns of the Jacobian.
    # This obviously not a pretty way to do this.
    # loglik. <-  -0.5*num.obs*log(2*pi) - num.obs*log(sigma.) - sum( (y - X %*% beta.)^2 )
    loglik. <- -sum( ( (y - X %*% beta.)/sigma. )^2 ) - 0.5*num.obs*log(2*pi) - num.obs*log(sigma.)
    
    # Take minus log-likelihood.
    objective. <- -loglik. / num.obs
    
    return( list( "objective" = objective.@F,
                  "gradient"  = as.vector( objective.@J ) ) )
}

# Find parameters that maximize the likelihood.
res <- nloptr( x0     = params0, 
               lb     = c( rep( -Inf, length(beta.true) ), 0.0001 ),
               ub     = c( rep(  Inf, length(beta.true) ), Inf ),
               eval_f = eval_f_list,
               opts   = opts,
               y      = y,
               X      = X )
print( res )

coefficients( lm( y ~ -1 + X ) )
