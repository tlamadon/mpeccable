
require(splines)
require(Matrix)

F_crra <- function(xsupp,ivals){

	n <- length(xsupp)

	ff <- function(ain,zin,gin){
		# ain: control variable (consumption)
		# zin: exogenous state (period, agg shock)
		# gin: parameter of the function (exponent of power utility). 1x1.

		if (is.fdiff(gin)){
			M <- matrix(0,nrow=length(zin), ncol=1) # that's gonna be th jacobian
			F <- array(0,length(zin))	# function value at each state

			for (i in ivals){
				I <- which(zin==i)
				J <- ((i-1)*n+1) : ((i-1)*n+n)
				F[I] <- F[I] + ain@F[I]^gin@F[J]
				M[I,J] <- diag( log(ain@F[I]) * exp( gin@F[J]*log(ain@F[I]) )
			}
			vars <- list(v1 = length(ivals) * n )
			names(vars) <- names(gin@vars[[1]])

			R <- new("FDiff",F=c(F),J=Matrix(M,sparse=T),vars=gin@vars)

			if (is.fdiff(ain)){
				M = matrix(0,nrow=n*length(ivals) , ncol = length(ain@F))
				
				for (i in ivals) { # quite inneficient!
				  I  = which(zin==i)
				  J  = 1:n # we need all functional parameters
				  if (getOption('mpeccable.coloring')) {
					M[I,I] = diag(length(I))  
				  } else {
					M[I,J] <- diag( gin@F[J] * ain@F[I] ^ (gin@F[J] - 1) )
				  }
				}
				R = appendJac(R,Matrix(M,sparse=T),ain@vars)
				return(R)
			}
		}
	}
	class(ff) = 'frep'
attr(ff,'ng') = length(ivals) * n
return(ff)
}

