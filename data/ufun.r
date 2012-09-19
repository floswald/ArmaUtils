


# data for ufun-example

pars <- list()
pars$alpha <- 0.6
pars$sigma <- 1.6
pars$cutoff <- 0.05
pars$theta  <- 1.1

fsize      <- seq(1,3,le=10)
xi <- list()
xi[[1]]    <- fsize^(pars$sigma -1)
xi[[2]]    <- xi[[1]] / (1-pars$sigma)
xi[[3]]    <- pars$alpha*(1-pars$sigma)
xi[[4]]    <- (1-pars$alpha)*(1-pars$sigma)
xi[[5]]    <- pars$alpha*(1-pars$sigma) -1 

Res <- outer(1:5,-4:8,FUN="+")
age <- 4

