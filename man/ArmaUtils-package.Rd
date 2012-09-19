\name{ArmaUtils-package}
\alias{ArmaUtils-package}
\alias{ArmaUtils}
\docType{package}
\title{Fast RcppArmadillo and RcppEigen Helper Functions for Dynamic Programming problems}

\description{Implements several functions useful for discrete and continuous state dynamic programming. A prototypical problem would be to approximate
V(S) = max_x u(S,x) + b * V(S'|S,x) where S stands for the state space, x for control variables, V() and u() are real valued functions. The package is useful if S is large. 
It has the following functionality: \cr
\tabular{ll}{
compute u(): \tab \code{\link{utilfun}} efficiently computes different forms of u(S,x) at each point of kronecker(S,x)\cr
max u(S,x) + b * V(S'|S,x): \tab \code{\link{armamax}} finds the pointwise maximum at each S, optionally observing a restriction on the search space \cr
approximate V(S'): \tab \code{\link{armakron}} computes approximating B-spline coefficients on high-dimensional kronecker products of univariate basis functions \cr
evaluate approximations of V(S'): \tab \code{\link{krons}} evaluate high-dimensional kronecker products to obtain an approximation to V(S')\cr
}}

\details{
\tabular{ll}{
Package: \tab ArmaUtils\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-11\cr
License: \tab What license is it under?\cr
}

}
\author{florian.oswald@gmail.com}

\references{C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173-182, 1979\cr
M. Miranda and P. Fackler. Applied computational economics and finance. MIT Press, 2002}
\keyword{ package }

\examples{
\dontrun{
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
}
}
