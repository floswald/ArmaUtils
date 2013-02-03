

\name{utilfun}
\alias{utilfun}
\title{Calculates Utility functions}

\description{Calculates CRRA utility of u(Resources, hsize, age, params, xi, own ) at each entry of Resources. Not used in Li, Meghir, Oswald.}

\usage{utilfun(Resources, hsize, age, params, xi, own )}

\arguments{\item{Resources}{a numeric matrix}
	       \item{hsize}{either a scalar or a vector of length \code{nrow(Resources)} of current house size}
		   \item{age}{integer age of individual}
		   \item{params}{list of scalar parameters:
		   \enumerate{\item theta: ownership premium
			          \item alpha: consumption share
				      \item cutoff: minimal consumption level c*. quadratic approximation for c<c*}}
		   \item{xi}{list of vector parameters}
		   \item{own}{boolean value for ownership}}

\details{computes CRRA utility of u(c,h,H,theta,t). see \link{doc/ufun.pdf}. Computes a quadratic approximation to u() for values of consumption such that c < c*, where c* > 0.}

\value{Numeric matrix of CRRA utility values}
\author{<florian.oswald@gmail.com>}

\examples{
\dontrun{
pars <- list(alpha=0.6,sigma=1.6,cutoff=0.05,theta=1.1)
maxage <- 10
fsize <- seq(1,3,le=maxage)
xi <- list()
xi[[1]]    <- fsize^(pars$sigma -1)
xi[[2]]    <- xi[[1]] / (1-pars$sigma)
xi[[3]]    <- pars$alpha*(1-pars$sigma)
xi[[4]]    <- (1-pars$alpha)*(1-pars$sigma)
xi[[5]]    <- pars$alpha*(1-pars$sigma) -1 

Res  <- Res <- outer(1:5,-4:8,FUN="+")

utilfun(Resources=Res,hsize=2,age=4,params=pars,xi=xi,own=TRUE)
utilfun(Resources=Res,hsize=2,age=4,params=pars,xi=xi,own=FALSE)
utilfun(Resources=Res,hsize=1:5,age=4,params=pars,xi=xi,own=FALSE)
utilfun(Resources=Res,hsize=1:6,age=4,params=pars,xi=xi,own=FALSE)	# error message: too many house sizes.
}}
	

