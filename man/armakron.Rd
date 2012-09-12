
\name{armakron}
\alias{armakron}
\docType{package}
\title{RcppArmadillo algorithm to compute large Kronecker Products}

\description{Implements the DeBoor (1979) algorithm to compute the product \code{A * y}, where \code{A} is the (large) kronecker product of M matrices \eqn{mat_1,...,mat_M}, each of dimension \eqn{n_i x n_i,i=1,\dots,M}, and \code{y} is a compatible vector of length \eqn{N = n_1*n_2*\dots*n_M}}

\usage{armamax(y,matrices)}

\arguments{\item{y}{a numeric vector}\item{matrices}{A list of square matrices}}

\details{This function is useful for approximation high dimensional functional spaces with basis functions. Suppose we want to approximate an M-dimensional object \eqn{f} that maps R^M into R,  on the tensor product \code{A} of univariate grids \eqn{x_i} of length \eqn{n_i,i=1,\dots,M} each. Denoting the size \eqn{n_i x n_i} basis matrix for dimension \eqn{i} by \eqn{mat_i}, this kronecker product is \code{A <- kronecker(mat_M,kronecker(mat_M-1,...,kronecker(mat_2,mat_1))...)}.\cr
	Approximating coefficients \code{coefs} may be obtained by solving the system \code{A * coefs = y}, where \eqn{y_{i,j,...,k} = f(x_{1,i},x_{2,j},\dots,x_{M,k})} are function values at the grid. Beyond a certain number of dimensions M, or number of data points \eqn{n_i} this is infeasible. Firstly because \code{A} does not fit in memory, and secondly because solving the system is very costly. \code{armakron} implements the efficient deBoor (1979) algorithm to compute \code{coefs}. The kronecker product is never formed, and the calculation involves a series of repeated matrix multiplications. To obtain coefficients, it is most efficient to supply \emph{inverse} basis matrices. This algorithm is orders of magnitude faster than for example setting up \code{A} as a product of sparse matrices (as is the case with spline basis), and to use \code{solve(A,y)} from a sparse library.}

\value{a N x 1 matrix of approximating spline coefficients}

\author{<florian.oswald@gmail.com>}

\keyword{iteration,multivariate}

\references{C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173-182, 1979\cr
M. Miranda and P. Fackler. Applied computational economics and finance. MIT Press, 2002}

\examples{
\dontrun{
x1 <- matrix(rnorm(25),nrow=5,ncol=5)
x2 <- matrix(rnorm(9),nrow=3,ncol=3)
x3 <- matrix(rnorm(4),nrow=2,ncol=2)
y  <- rnorm(5*3*2)

R.kron   <- kronecker(x3,kronecker(x2,x1))	# form R kronecker product. impossible beyond a certain size
R.result <- R.kron \%*\% y

Arma.result <- armakron(y=y,list(x1,x2,x3))	# armakron() computes result directly without forming kronecker product

all.equal(R.result,Arma.result)		# TRUE

# spline example: estimate spline coefficients on kronecker product of 3 univariate spline basis x,y and z.
library(splines)

# degrees
dx <- 3
dy <- 3
dz <- 3

# knots
kx <- c(rep(0,times=dx),seq(0,1,length=17),rep(1,times=dx))
ky <- c(rep(1,times=dy),seq(1,10,length=8),rep(10,times=dy))
kz <- c(rep(-3,times=dz),seq(-3,3,length=11),rep(3,times=dz))

# evaluation points: choose such that # of datapoints is equal # of coefficients
# the algorithm works ONLY with square systems, i.e. number of data points must be equal number of spline coefficients!
x <- seq(0,1,length=length(kx)-dx-1)
y <- seq(1,10,length=length(ky)-dy-1)
z <- seq(-3,3,length=length(kz)-dz-1)

# basis matrices
X <- splineDesign(x=x,knots=kx,ord=dx+1)
Y <- splineDesign(x=y,knots=ky,ord=dy+1)
Z <- splineDesign(x=z,knots=kz,ord=dz+1)

# data generating process
dgp <- function(x) x[1]^2 + x[2]^3 + x[1]*x[2]*x[3]

# create sample of data: note that ordering of data is crucial. grid order and order implied by kronecker product must cooincide.
vals <- apply(expand.grid(x,y,z),1,dgp)

# plot a slice
persp(array(vals,dim=c(length(x),length(y),length(z)))[ , ,5],theta=100)

# estimate spline coefficients to be able to evaluate dgp off the grid
# the problem is: 
# kronprod \%*\% coefs = vals, where kronprod = kronecker(Z,kronecker(Y,X))
#                coefs = inv(kronprod) \%*\% vals
#
# often the forming of kronprod is infeasible because of memory limitations. If not that, then computing the inverse
# is very costly. The algorithm in armakron never forms the kronprod. Instead we give it the inverse matrices one by one:

library(MASS)	# adds ginv()
coefs <- armakron(y=vals,matrices=list(ginv(X),ginv(Y),ginv(Z)))

# check against R
kron <- kronecker(Z,kronecker(Y,X))
r.coefs <- solve(kron,vals)	# takes couple of seconds already.

all.equal(as.numeric(coefs),r.coefs)	# TRUE

# predict the grid. 
pred <- armakron(y=coefs,matrices=list(X,Y,Z))
all.equal(as.numeric(pred),vals)	# TRUE

}}
	

