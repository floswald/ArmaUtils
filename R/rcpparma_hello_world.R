

### Roxygen Documentation: Package description
#' Fast RcppArmadillo and RcppEigen Helper Functions for Dynamic Programming problems
#'
#' \tabular{ll}{
#' Package: \tab ArmaUtils\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2012-09-11\cr
#' License: \tab GPL (>=2)\cr
#' }
#'
#' Implements several functions useful for discrete and continuous state dynamic programming. A prototypical problem would be to approximate\cr
#' V(S) = max_x u(S,x) + b * V(S'|S,x) \cr
#' where S stands for the state space, x for control variables, V() and u() are real valued functions. The package is useful if S is large. \cr
#' It has the following functionality: 
#'
#' \tabular{ll}{
#' compute u(): \tab \code{\link{utilfun}} efficiently computes different forms of u(S,x) at each point of kronecker(S,x)\cr
#' max u(S,x) + b * V(S'|S,x): \tab \code{\link{armamax}} finds the pointwise maximum at each S, optionally observing a restriction on the search space \cr
#' approximate V(S'): \tab \code{\link{armakron}} computes approximating B-spline coefficients on high-dimensional kronecker products of univariate basis functions \cr
#' evaluate approximations of V(S'): \tab \code{\link{krons}} evaluate high-dimensional kronecker products to obtain an approximation to V(S')\cr
#' }
#'
#' @name ArmaUtils-package
#' @aliases ArmaUtils
#' @docType package
#' @title Fast RcppArmadillo and RcppEigen Helper Functions for Dynamic Programming problems
#' @author Florian Oswald \email{florian.oswald@@gmail.com}
#' @keyword package
#' @seealso \code{\link{utilfun}}, \code{\link{armamax}}, \code{\link{armakron}}, \code{\link{krons}}
NULL



#' RcppArmadillo row-wise Matrix Optimizer 
#'
#' Finds the row-wise maximum of a matrix \code{A}, and returns a list with optimal values and their position in each row. 
#' The returned value is similar to \code{list(values=apply(A,1,max),indices=apply(A,1,which.max)}, but much faster than this.\cr 
#' A search range delimiter may be specified.
#'
#' Takes a matrix \code{A} and returns the row-wise maximum and it's position in the row as a list. Optionally one can restrict the region that \cr
#' is searched over for each row by specifying an integer vector \code{b}. \code{b} must be of length \code{nrow(A)}, and at each position hold a \cr
#' value in \eqn{0,\dots,ncol(A)-1}. This is useful if there is some model restriction that defines values in \code{A[i, ]} to the left of, and \cr 
#' including \code{b[i]} as inadmissible. If optimal values and their position are required, \code{armamax} does only one calculation for each \cr
#' row (it finds the maximum and remembers it's position at the same time), where \code{list(values=max(A),indices=which.max(A)} needs two. \cr
#' This can be seen as an \code{R} version of the MATLAB idiom \code{[index, value] = max(A,[],2)}. On top of that the \code{RcppArmadillo} \cr
#' implementation of the loop over rows of \code{A} is much faster than the one implied by \code{apply()}. See \code{inst/tests.r} for some benchmarks.}
#'
#' @param A numeric matrix of dimension (n,m) for which the maximum value and index should be found on each row
#' @param b integer vector of length(n) holding index j in row i. Search over A[i,(j+1:m)]
#' @return a list with entries \item{values}{a vector of length \code{nrow(A)} with the row-wise maxima of \code{A}}\item{indices}{a vector of length \code{nrow(A)} with the row-wise maxima of \code{A}}
#' @author Florian Oswald \email{florian.oswald@@gmail.com}
#' @title RcppArmadillo row-wise Matrix Optimizer
#' @examples
#' \dontrun{
#' ## find the unrestricted row maximum
#' x1 <- matrix(rnorm(25),nrow=5,ncol=5)
#' result <- armamax(x1)
#' all.equal(apply(x1,1,max),as.numeric(result$values))	# TRUE
#' all.equal(apply(x1,1,which.max),as.numeric(result$indices))	# TRUE
#'
#' ## find the row maximum restricted to lie to the right of b in each row
#' A <- matrix(replicate(n=4,sample(1:5,size=5,replace=F)),nrow=4,ncol=5,byrow=T)
#' b <- sample(0:(ncol(A)-1),size=nrow(A),replace=TRUE)
#'
#' rfun <- function(A,b){
#' 	R.result <- list()
#' 	R.result$values  <- c()
#'	R.result$indices <- c()
#'	for (i in 1:nrow(A)){
#'		R.result$values[i]  <- max(A[i,(b[i]+1):ncol(A)])
#'		R.result$indices[i] <- which.max(A[i,(b[i]+1):ncol(A)]) + b[i]
#'	}
#'	return(R.result)
#' } 
#'
#' all.equal(R.result$values,as.numeric(armamax(A,b)$values))
#' all.equal(R.result$indices,as.numeric(armamax(A,b)$indices))}
armamax <- function(A,b=NULL){
	if (is.null(b)) {
		.Call( "armamax", A, PACKAGE = "ArmaUtils" )
	} else {
		.Call( "armamax2", A, b, PACKAGE = "ArmaUtils" )
	}
}


#' RcppArmadillo algorithm to compute large Kronecker Products 
#'
#' Implements the DeBoor (1979) algorithm to compute the product \code{A * y}, where \code{A} is the (large) kronecker product of M matrices \eqn{mat_1,...,mat_M}, each of dimension \eqn{n_i x n_i,i=1,\dots,M}, and \code{y} is a compatible vector of length \eqn{N = n_1*n_2*\dots*n_M}}
#'
#' This function is useful for approximation high dimensional functional spaces with basis functions. Suppose we want to approximate an M-dimensional object \eqn{f} that maps R^M into R,  on the tensor product \code{A} of univariate grids \eqn{x_i} of length \eqn{n_i,i=1,\dots,M} each. Denoting the size \eqn{n_i x n_i} basis matrix for dimension \eqn{i} by \eqn{mat_i}, this kronecker product is \code{A <- kronecker(mat_M,kronecker(mat_M-1,...,kronecker(mat_2,mat_1))...)}.\cr
#' Approximating coefficients \code{coefs} may be obtained by solving the system \code{A * coefs = y}, where \eqn{y_{i,j,...,k} = f(x_{1,i},x_{2,j},\dots,x_{M,k})} are function values at the grid. Beyond a certain number of dimensions M, or number of data points \eqn{n_i} this is infeasible. Firstly because \code{A} does not fit in memory, and secondly because solving the system is very costly. \code{armakron} implements the efficient deBoor (1979) algorithm to compute \code{coefs}. The kronecker product is never formed, and the calculation involves a series of repeated matrix multiplications. To obtain coefficients, it is most efficient to supply \emph{inverse} basis matrices. This algorithm is orders of magnitude faster than for example setting up \code{A} as a product of sparse matrices (as is the case with spline basis), and to use \code{solve(A,y)} from a sparse library.}
#'
#' @param matrices list of M \emph{square} matrices 
#' @param y numeric vector with \code{length(y) = prod(unlist(lapply(matrices,nrow)))}
#' @return length(y) x 1 matrix of approximating spline coefficients
#' @author Florian Oswald \email{florian.oswald@@gmail.com}
#' @keyword iteration,multivariate
#' @references C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173-182, 1979\cr
#' M. Miranda and P. Fackler. Applied computational economics and finance. MIT Press, 2002
#' @examples
#' \dontrun{
#' x1 <- matrix(rnorm(25),nrow=5,ncol=5)
#' x2 <- matrix(rnorm(9),nrow=3,ncol=3)
#' x3 <- matrix(rnorm(4),nrow=2,ncol=2)
#' y  <- rnorm(5*3*2)
#' 
#' R.kron   <- kronecker(x3,kronecker(x2,x1))	# form R kronecker product. impossible beyond a certain size
#' R.result <- R.kron \%*\% y
#' 
#' Arma.result <- armakron(y=y,list(x1,x2,x3))	# armakron() computes result directly without forming kronecker product
#' 
#' all.equal(R.result,Arma.result)		# TRUE
#' 
#' # spline example: estimate spline coefficients on kronecker product of 3 univariate spline basis x,y and z.
#' library(splines)
#' 
#' # degrees
#' dx <- 3
#' dy <- 3
#' dz <- 3
#' 
#' # knots
#' kx <- c(rep(0,times=dx),seq(0,1,length=17),rep(1,times=dx))
#' ky <- c(rep(1,times=dy),seq(1,10,length=8),rep(10,times=dy))
#' kz <- c(rep(-3,times=dz),seq(-3,3,length=11),rep(3,times=dz))
#' 
#' # evaluation points: choose such that # of datapoints is equal # of coefficients
#' # the algorithm works ONLY with square systems, i.e. number of data points must be equal number of spline coefficients!
#' x <- seq(0,1,length=length(kx)-dx-1)
#' y <- seq(1,10,length=length(ky)-dy-1)
#' z <- seq(-3,3,length=length(kz)-dz-1)
#' 
#' # basis matrices
#' X <- splineDesign(x=x,knots=kx,ord=dx+1)
#' Y <- splineDesign(x=y,knots=ky,ord=dy+1)
#' Z <- splineDesign(x=z,knots=kz,ord=dz+1)
#'
#' # data generating process
#' dgp <- function(x) x[1]^2 + x[2]^3 + x[1]*x[2]*x[3]
#' 
#' # create sample of data: note that ordering of data is crucial. grid order and order implied by kronecker product must cooincide.
#' vals <- apply(expand.grid(x,y,z),1,dgp)
#' 
#' # plot a slice
#' persp(array(vals,dim=c(length(x),length(y),length(z)))[ , ,5],theta=100)
#' 
#' # estimate spline coefficients to be able to evaluate dgp off the grid
#' # the problem is: 
#' # kronprod \%*\% coefs = vals, where kronprod = kronecker(Z,kronecker(Y,X))
#' #                coefs = inv(kronprod) \%*\% vals
#' #
#' # often the forming of kronprod is infeasible because of memory limitations. If not that, then computing the inverse
#' # is very costly. The algorithm in armakron never forms the kronprod. Instead we give it the inverse matrices one by one:
#' 
#' library(MASS)	# adds ginv()
#'coefs <- armakron(y=vals,matrices=list(ginv(X),ginv(Y),ginv(Z)))
#'
#' # check against R
#' kron <- kronecker(Z,kronecker(Y,X))
#' r.coefs <- solve(kron,vals)	# takes couple of seconds already.
#'
#' all.equal(as.numeric(coefs),r.coefs)	# TRUE
#' 
#' # predict the grid. 
#' pred <- armakron(y=coefs,matrices=list(X,Y,Z))
#' all.equal(as.numeric(pred),vals)	# TRUE
#' }
armakron <- function(y,matrices){
	.Call( "armakron", y, matrices, PACKAGE = "ArmaUtils" )
}





utilfun <- function( Resources, hsize, age, params, xi, own ){
	if (is.null(dim(Resources))) Resources <- matrix(Resources,nrow=length(Resources),ncol=1)
	.Call( "utilfun", Resources, hsize, age, params, xi, own, PACKAGE = "ArmaUtils" )
}


ufun_Atta <- function( Resources, hsize, params ){
	if (is.null(dim(Resources))) Resources <- matrix(Resources,nrow=length(Resources),ncol=1)
	if (is.null(dim(hsize))) hsize <- matrix(hsize,nrow=nrow(Resources),ncol=1)
	.Call( "ufun_Attanasio", Resources, hsize, params, PACKAGE = "ArmaUtils" )
}

# takes y and a list of matrices to compute kronecker of matrics times y. matrices can be sparse or dense (not implemented)
krons <- function( y, matrices ){
	stopifnot(is.list(matrices))
	n <- length(matrices)
	alldgC <- lapply(1:n,function(i) class(matrices[[i]])[1]=="dgCMatrix")
	stopifnot(all(unlist(alldgC)))
	if (n==4) .Call( "kron4", matrices[[1]], matrices[[2]],matrices[[3]],matrices[[4]], y, PACKAGE = "ArmaUtils" )
}


