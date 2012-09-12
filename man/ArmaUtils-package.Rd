\name{ArmaUtils-package}
\alias{ArmaUtils-package}
\alias{ArmaUtils}
\docType{package}
\title{Fast RcppArmadillo Helper Functions}

\description{Implements 2 functions \code{armamax} and \code{armakron}. The first one is a fast implementation of \code{list(values=apply(A,1,max),indice=apply(A,1,which.max)} each row. (If you want it for each column, just transpose the matrix). 
One can optionally specify an index vector \code{inds}, that limits the search range in each row j, i.e. instead of comparing all indices, 
only the ones after \code{inds(j)} are taken into account.\cr 
The second is \code{armakron}, which is an efficient algorithm for calculating the product \code{A * y}, where \code{A} is the (large) kronecker product of M matrices \eqn{mat_1,...,mat_M}, each of dimension \eqn{n_i x n_i,i=1,\dots,M}, and \code{y} a compatible vector of length \eqn{N = n_1*n_2*\dots*n_M}}

\details{
\tabular{ll}{
Package: \tab ArmaUtils\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-11\cr
License: \tab What license is it under?\cr
}
See \code{"\link[=ArmaUtils]{armamax}"} and \code{"\link[=ArmaUtils]{armakron}"} for help on both functions.

}
\author{florian.oswald@gmail.com}

\references{C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173-182, 1979}
\keyword{ package }

\examples{
## not run
# armakron()

x1 <- matrix(rnorm(25),nrow=5,ncol=5)
x2 <- matrix(rnorm(9),nrow=3,ncol=3)
x3 <- matrix(rnorm(4),nrow=2,ncol=2)
y  <- rnorm(5*3*2)

R.kron   <- kronecker(x3,kronecker(x2,x1))	# form R kronecker product. impossible beyond a certain size
R.result <- R.kron %*% y

Arma.result <- armakron(y=y,list(x1,x2,x3))	# armakron() computes result directly without forming kronecker product

all.equal(R.result,Arma.result)		# TRUE

# armamax()
result <- armamax(x1)
all.equal(apply(x1,1,max),as.numeric(result$values))	# TRUE
all.equal(apply(x1,1,which.max),as.numeric(result$indices))	# TRUE

## end not run
}
