

\name{ufun_Atta}
\alias{ufun_Atta}
\title{Calculates Utility functions from Wakefield et al REDyn 2012}

\description{Calculates CRRA utility from u(Resources, hsize, params )}

\usage{utilfun(Resources, hsize, params)}

\arguments{\item{Resources}{a numeric matrix}
	       \item{hsize}{either a scalar or a vector of length \code{nrow(Resources)} of current house size}
		   \item{params}{list of scalar parameters:
		   \enumerate{\item theta: elasticity of substitution c,h
			          \item gamma: coefficient of relative risk aversion
			          \item phival: fraction of house-utility derived from flat
			          \item mu: ownership premium
				      \item cutoff: minimal consumption level c*. quadratic approximation for c<c*}}}
\details{computes CRRA utility of u(c,h,params).

\deqn{\begin{array}{ll}
	u(c,h) &= \frac{c^{1-\gamma}}{1-\gamma} \exp( \theta \phi(h) ) + \mu \phi(h) \\
	\phi(h) &= \left\{
		\begin{array}{ll}
		0                       &   h=0 \\
		\mbox{phival} \in (0,1) &   h=1 \\
		1                       &   h=2 
		\end{array}
		\right. \\
	\mbox{params} &= \left( \gamma, \theta, \mbox{phival}, \mu \right)
	\end{array}}{ u(c,h) = c^1-gamma / (1-gamma) * exp( theta * phi(h) ) + mu * phi(h), phi(h) = 0 if h=0, phi(h) = phival if h=1, phi(h) = 1 if h=2 }
Computes a quadratic approximation to u() for values of consumption such that c < c*, where c* > 0. Resoures can be a scalar, a vector or a matrix. Utility is computed on all elements. Some care needs to be taken when supplying a matrix of consumption values: the corresponding vector of house sizes needs to be either a scalar (same house size at all cons values) or of length \code{nrow(Resources)}}

\value{Numeric matrix of CRRA utility values}
\author{<florian.oswald@gmail.com>}

\examples{
\dontrun{
	examples
}}


