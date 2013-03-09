

# R implementation of C functions

# ufun_labouR

ufun_labouR <- function(e,w,params){
	xi1       <- params$alpha * (1-params$gamma);
	xi2       <- (1-params$alpha)*(1-params$gamma);
	alphaxi1  <- params$alpha^xi1
	malphaxi2 <- (1-params$alpha)^xi2
	if (is.null(dim(e))) e <- matrix(e,length(e),1)
	w <- matrix(w,dim(e)[1],dim(e)[2])
	# return objects
	util <- matrix(0,dim(e)[1],dim(e)[2])
	cons <- matrix(0,dim(e)[1],dim(e)[2])
	labo <- cons
	# where work?
	work <- (1-params$alpha)*e < w
	# g(e,w,n)
	gw    <- with(params,( alphaxi1 * malphaxi2 * e^(1-gamma)) / w^xi2)
	gwx   <- with(params,( alphaxi1 * malphaxi2 * cutoff^(1-gamma)) / w^xi2)
	distw <- e-params$cutoff
	# g(e,w,-n)
	consmat <- e - w
	gnw     <- with(params,consmat^xi1)
	gnwx    <- with(params,cutoff^xi1)
	distnw  <- consmat-params$cutoff
	for (i in 1:nrow(e)){
		for (j in 1:ncol(e)){
			if (work[i,j]){
				# work. 
				if (e[i,j]>params$cutoff){
					# positive consumption?
					util[i,j] <- gw[i,j] / (1-params$gamma)
					cons[i,j] <- params$alpha * e[i,j]
					labo[i,j] <- 1 - (1-params$alpha) * e[i,j]/w[i,j]
				} else {
					grad      <- gwx[i,j] / params$cutoff
					hess      <- (-params$gamma/params$cutoff) * grad
					util[i,j] <- gwx[i,j] / (1-params$gamma) + grad*distw[i,j] + 0.5 * hess * distw[i,j]^2
					cons[i,j] <- e[i,j]
				}
			} else {
				# don't work. 
				if (consmat[i,j]>params$cutoff){
					# pos cons?
					util[i,j] <- gnw[i,j] / (1-params$gamma)
					cons[i,j] <- consmat[i,j]
				} else {
					grad      <- params$alpha * gnwx / params$cutoff
					hess      <- with(params,(alpha*(1-gamma)-1) * grad / cutoff)
					util[i,j] <- gnwx / (1-params$gamma) + grad*distnw[i,j] + 0.5*hess*distnw[i,j]^2
					cons[i,j] <- consmat[i,j]
				}
			}
		}
	}
	return(list(utility=util,consumption=cons,labour=labo))
}


# ufun_labouR_h
# same but with housing

ufun_labouR_h <- function(e,w,s,params){
	xi1       <- params$alpha * (1-params$gamma);
	xi2       <- (1-params$alpha)*(1-params$gamma);
	alphaxi1  <- params$alpha^xi1
	malphaxi2 <- (1-params$alpha)^xi2
	if (is.null(dim(e))) e <- matrix(e,length(e),1)
	w <- matrix(w,dim(e)[1],dim(e)[2])
	s <- matrix(s,dim(e)[1],1)
	# return objects
	util <- matrix(0,dim(e)[1],dim(e)[2])
	cons <- matrix(0,dim(e)[1],dim(e)[2])
	labo <- cons
	# where work?
	work <- (1-params$alpha)*e < w
	# g(e,w,n)
	gw    <- with(params,( alphaxi1 * malphaxi2 * e^(1-gamma)) / w^xi2)
	gwx   <- with(params,( alphaxi1 * malphaxi2 * cutoff^(1-gamma)) / w^xi2)
	distw <- e-params$cutoff
	# g(e,w,-n)
	consmat <- e - w
	gnw     <- with(params,consmat^xi1)
	gnwx    <- with(params,cutoff^xi1)
	distnw  <- consmat-params$cutoff
	for (i in 1:nrow(e)){
		for (j in 1:ncol(e)){
			if (work[i,j]){
				# work. 
				if (e[i,j]>params$cutoff){
					# positive consumption?
					util[i,j] <- gw[i,j] / (1-params$gamma)
					cons[i,j] <- params$alpha * e[i,j]
					labo[i,j] <- 1 - (1-params$alpha) * e[i,j]/w[i,j]
				} else {
					grad      <- gwx[i,j] / params$cutoff
					hess      <- (-params$gamma/params$cutoff) * grad
					util[i,j] <- gwx[i,j] / (1-params$gamma) + grad*distw[i,j] + 0.5 * hess * distw[i,j]^2
					cons[i,j] <- e[i,j]
					#                     labo[i,j] <- 1 - (1-params$alpha) * e[i,j]/w[i,j]
				}
			} else {
				# don't work. 
				if (consmat[i,j]>params$cutoff){
					# pos cons?
					util[i,j] <- gnw[i,j] / (1-params$gamma)
					cons[i,j] <- consmat[i,j]
				} else {
					grad      <- params$alpha * gnwx / params$cutoff
					hess      <- with(params,(alpha*(1-gamma)-1) * grad / cutoff)
					util[i,j] <- gnwx / (1-params$gamma) + grad*distnw[i,j] + 0.5*hess*distnw[i,j]^2
					cons[i,j] <- consmat[i,j]
				}
			}
		}
	}
	phivals <- c(0,params$phival,1)
	phivec <- rep(0,length(s))
	for (i in 1:length(s)){
		phivec[i] <- phivals[ s[i] + 1 ]
	}
	phimat <- matrix(phivec,length(s),ncol(e))
	hfac <- exp( params$theta * phimat )
	util <- util * hfac
	util <- util + params$mu * phimat
	return(list(utility=util,consumption=cons,labour=labo))
}

