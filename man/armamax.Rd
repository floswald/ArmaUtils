\name{armamax}
\alias{armamax}
\title{RcppArmadillo row-wise Matrix Optimizer}

\description{Finds the row-wise maximum of a matrix \code{A}, and returns a list with optimal values and their position in each row. 
	The returned value is similar to \code{list(values=apply(A,1,max),indices=apply(A,1,which.max)}, but much faster than this. A search range delimiter may be specified.}

\usage{armamax(A,b)}

\arguments{\item{A}{a numeric matrix}\item{b}{an optional vector of integer indices limiting the search range}}



\details{Takes a matrix \code{A} and returns the row-wise maximum and it's position in the row as a list. Optionally one can restrict the region that is searched over for each row by specifying an integer vector \code{b}. \code{b} must be of length \code{nrow(A)}, and at each position hold a value in \eqn{0,\dots,ncol(A)-1}. This is useful if there is some model restriction that defines values in \code{A[i, ]} to the left of, and including \code{b[i]} as inadmissible. If optimal values and their position are required, \code{armamax} does only one calculation for each row (it finds the maximum and remembers it's position at the same time), where \code{list(values=max(A),indices=which.max(A)} needs two. This can be seen as an \code{R} version of the MATLAB idiom \code{[index, value] = max(A,[],2)}. On top of that the \code{RcppArmadillo} implementation of the loop over rows of \code{A} is much faster than the one implied by \code{apply()}. See \code{inst/tests.r} for some benchmarks.}

\value{a list with entries \item{values}{a vector of length \code{nrow(A)} with the row-wise maxima of \code{A}}\item{indices}{a vector of length \code{nrow(A)} with the row-wise maxima of \code{A}}}

\author{<florian.oswald@gmail.com>}

\keyword{array,iteration}

\examples{
\dontrun{
# find the unrestricted row maximum
x1 <- matrix(rnorm(25),nrow=5,ncol=5)
result <- armamax(x1)
all.equal(apply(x1,1,max),as.numeric(result$values))	# TRUE
all.equal(apply(x1,1,which.max),as.numeric(result$indices))	# TRUE


# find the row maximum restricted to lie to the right of b in each row
A <- matrix(replicate(n=4,sample(1:5,size=5,replace=F)),nrow=4,ncol=5,byrow=T)
b <- sample(0:(ncol(A)-1),size=nrow(A),replace=TRUE)

rfun <- function(A,b){
	R.result <- list()
	R.result$values  <- c()
	R.result$indices <- c()
	for (i in 1:nrow(A)){
		R.result$values[i]  <- max(A[i,(b[i]+1):ncol(A)])
		R.result$indices[i] <- which.max(A[i,(b[i]+1):ncol(A)]) + b[i]
	}
	return(R.result)
}

all.equal(R.result$values,as.numeric(armamax(A,b)$values))
all.equal(R.result$indices,as.numeric(armamax(A,b)$indices))}}
	

