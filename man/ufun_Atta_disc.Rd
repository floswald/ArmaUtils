

\name{ufun_Atta_disc}
\alias{ufun_Atta_disc}
\title{CRRA Utility function with discrete labor supply and housing}

\description{Calculates CRRA utility from u(Resources, hsize, labour, params )}

\usage{ufun_Atta_disc(Resources, hsize, labour, params)}

\arguments{\item{Resources}{a numeric matrix}
	       \item{labour}{vector of length \code{nrow(Resources)} of discrete amounts of labour supplied}
	       \item{hsize}{either a scalar or a vector of length \code{nrow(Resources)} of current house size}
		   \item{params}{list of scalar parameters:
		   \enumerate{\item theta: elasticity of substitution c,h
			          \item alpha: utility weight of labour 
			          \item gamma: coefficient of relative risk aversion
			          \item phival: fraction of house-utility derived from flat
			          \item mu: ownership premium
				      \item cutoff: minimal consumption level c*. quadratic approximation for c<c*}}}
\details{computes augmented CRRA utility of u(c,l,h,params).

\deqn{\begin{array}{ll}
	u(c,l,h) &= \frac{\left(c \exp( \alpha l ) \right)^{1-\gamma}}{1-\gamma} \exp( \theta \phi(h) ) + \mu \phi(h) \\
	\phi(h) &= \left\{
		\begin{array}{ll}
		0                       &   h=0 \\
		\mbox{phival} \in (0,1) &   h=1 \\
		1                       &   h=2 
		\end{array}
		\right. \\
	\mbox{params} &= \left( \gamma, \alpha, \theta, \mbox{phival}, \mu \right)
	\end{array}}{ u(c,l,h) = (c * exp(alpha * l) )^1-gamma / (1-gamma) * exp( theta * phi(h) ) + mu * phi(h), phi(h) = 0 if h=0, phi(h) = phival if h=1, phi(h) = 1 if h=2 }
Computes a quadratic approximation to u() for values of consumption such that c < c*, where c* > 0. Resoures can be a scalar, a vector or a matrix. Utility is computed on all elements. Some care needs to be taken when supplying a matrix of consumption values: the corresponding vector of house sizes needs to be either a scalar (same house size at all cons values) or of length \code{nrow(Resources)}}

\value{Numeric Matrix of Utility values.}
\author{<florian.oswald@gmail.com>}
\references{Utility functions in 
	\tabular{l}{
	Low, Meghir, Pistaferri, "Wage Risk and Employment Risk over the Life Cycle", American Economic Review, 2010 and \cr
	Attanasio, Bottazzi, Low, Nesheim, Wakefield, "Modelling the Demand for Housing over the Lifecycle", Review of Economic Dynamics, 2012 
	}}

\examples{
pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(1:4,5:-1)	# consumption
lab <- c(0,1,0,1)	# discrete labor supply choices
s <- c(0,0,1,2)	# house sizes
ufun_Atta_disc(Res,s,lab,pars)
}


