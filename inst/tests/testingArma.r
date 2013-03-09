
rm(list=ls(all=T))
library(ArmaUtils)
library(rbenchmark)

# test file for ArmaUtils package
# ===============================

# Armamax
# -------

#### Useage 1: find the unrestricted maximal element for each row of a matrix

context("testing armamax")
A <- matrix(rnorm(200000),nrow=400,ncol=500)
result <- armamax(A)
test_that("A is a list",{
		  expect_that( is.list(result), is_true() )})
test_that("component 1 of A is a numeric",{
		  expect_that( is.numeric(result[[1]]), is_true() )})
test_that("component 2 of A is a integer",{
		  expect_that( is.integer(result[[2]]), is_true() )})
test_that("result$values is equal to apply(A,1,max)",{
		  expect_that( all.equal(apply(A,1,max),as.numeric(result$values)), is_true() )})
test_that("result$indices is equal to apply(A,1,which.max)",{
		  expect_that( all.equal(apply(A,1,which.max),as.numeric(result$indices)), is_true() )})


# benchmark against 2 calls to apply(), one for values, one for indices
bmk.cols <- c("test","replications","elapsed","relative","user.self","sys.self")
benchmark(R=list(apply(A,1,max),apply(A,1,which.max)),armamax=armamax(A),replications=50,columns=bmk.cols)

# benchmark against 1 call to apply() only
benchmark(R=apply(A,1,max),armamax=armamax(A),replications=50,columns=bmk.cols)



#### Useage 2: find the REstricted maximal element for each row of a matrix

# there is a vector b of length nrow(A), holding an integer in 0,...,ncol(A)-1 that indicates how many indices starting from the left are to be excluded from the maximization process. E.g. if b(i) = 3, this means that the routine does max(A[i,-c(1,2,3)])

# A <- matrix(replicate(n=4,sample(1:5,size=5,replace=F)),nrow=4,ncol=5,byrow=T)
# A
# 
# b <- sample(0:(ncol(A)-1),size=nrow(A),replace=TRUE)
# b
# 
# rfun <- function(A,b){
#     R.result <- list()
#     R.result$values  <- c()
#     R.result$indices <- c()
#     for (i in 1:nrow(A)){
#         R.result$values[i]  <- max(A[i,(b[i]+1):ncol(A)])
#         R.result$indices[i] <- which.max(A[i,(b[i]+1):ncol(A)]) + b[i]
#     }
#     return(R.result)
# }
# 
# all.equal(R.result$values,as.numeric(armamax(A,b)$values))
# all.equal(R.result$indices,as.numeric(armamax(A,b)$indices))
# 
# A <- matrix(replicate(n=4000,sample(1:500,size=500,replace=F)),nrow=4000,ncol=500,byrow=T)
# 
# b <- sample(0:(ncol(A)-1),size=nrow(A),replace=TRUE)
# 
# benchmark(R=rfun(A,b),arma=armamax(A,b),replications=5,columns=bmk.cols)
# 


# Armakron
# --------

x1 <- matrix(rnorm(25),nrow=5,ncol=5)
x2 <- matrix(rnorm(9),nrow=3,ncol=3)
x3 <- matrix(rnorm(4),nrow=2,ncol=2)
y  <- rnorm(5*3*2)

R.kron   <- kronecker(x3,kronecker(x2,x1))	# form R kronecker product. impossible beyond a certain size
R.result <- R.kron %*% y
Arma.result <- armakron(y=y,list(x1,x2,x3))	# armakron() computes result directly without forming kronecker product

context("testing Armakron")
test_that("armakron is equal to R kronecker times vector",{
		  expect_that( all.equal(R.result,Arma.result), is_true() )})



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

# create sample of data:
vals <- apply(expand.grid(x,y,z),1,dgp)

# plot a slice
# persp(array(vals,dim=c(length(x),length(y),length(z)))[ , ,5],theta=100)

# estimate spline coefficients to be able to evaluate dgp off the grid
# the problem is: 
# kronprod %*% coefs = vals, where kronprod = kronecker(Z,kronecker(Y,X))
#              coefs = inv(kronprod) %*% vals
#
# often the forming of kronprod is infeasible because of memory limitations. If not that, then computing the inverse
# is very costly. The algorithm in armakron never forms the kronprod. Instead we give it the inverse matrices one by one:

library(MASS)	# adds ginv()
coefs <- armakron(y=vals,matrices=list(ginv(X),ginv(Y),ginv(Z)))

# check against R
kron <- kronecker(Z,kronecker(Y,X))
r.coefs <- solve(kron,vals)	# takes couple of seconds already.

test_that("armakron is equal to R kronecker times vector, test 2",{
		  expect_that( all.equal(as.numeric(coefs),r.coefs), is_true() )})

# predict the grid. 
pred <- armakron(y=coefs,matrices=list(X,Y,Z))
test_that("predictions from armakron coefficients are the same as original values",{
		  expect_that( all.equal(as.numeric(pred),vals), is_true() )})








# test files
library(ArmaUtils)

# test utility function 1
# =======================

# pars <- list(alpha=0.6,sigma=1.6,cutoff=0.05,theta=1.1)
# maxage <- 10
# fsize <- seq(1,3,le=maxage)
# xi <- list()
# xi[[1]]    <- fsize^(pars$sigma -1)
# xi[[2]]    <- xi[[1]] / (1-pars$sigma)
# xi[[3]]    <- pars$alpha*(1-pars$sigma)
# xi[[4]]    <- (1-pars$alpha)*(1-pars$sigma)
# xi[[5]]    <- pars$alpha*(1-pars$sigma) -1 
# 
# Res  <- outer(1:4,5:-1)
# 
# utilfun(Resources=Res,hsize=2,age=4,params=pars,xi=xi,own=TRUE)
# utilfun(Resources=Res,hsize=2,age=4,params=pars,xi=xi,own=FALSE)
# utilfun(Resources=Res,hsize=1:4,age=4,params=pars,xi=xi,own=FALSE)
# utilfun(Resources=Res,hsize=1:5,age=4,params=pars,xi=xi,own=FALSE)	# error message: too many house sizes.



# test utility function 2
# =======================
pars <- list(theta=0.115,gamma=1.4,phival=0.9,mu=0.26,cutoff=0.1)
Res  <- outer(1:4,5:-1)

# on a datatable (i.e. a vector)
library(data.table)
dat <- data.table(expand.grid(res = 1:10, size=0:2))
dat[,util := ufun_Atta(Resources = res, hsize = size, params=pars)]

# on a matrix
h1 <- ufun_Atta(Resources=Res,hsize=2,params=pars)	
h2 <- ufun_Atta(Resources=Res,hsize=c(1,2),params=pars)
h3 <- ufun_Atta(Resources=Res,hsize=rep(2L,nrow(Res)),params=pars)

context("testing utility function 2")
test_that("ufun_Atta returns a numeric",{
		  expect_that( is.numeric(dat[,util]), is_true() )})
test_that("ufun_Atta returns a matrix for a scalar house size",{
		  expect_that( is.matrix(h1), is_true() )})
test_that("ufun_Atta returns a matrix for a vector house size",{
		  expect_that( is.matrix(h2), is_true() )})
test_that("ufun_Atta is different for different house sizes",{
		  expect_that( all.equal(h1,h2)==TRUE,is_false())})
test_that("ufun_Atta returns a numeric",{
		  expect_that(ufun_Atta(Resources=Res,hsize=rep(2L,nrow(Res)),params=c(1,2))
,throws_error())})



# test dependence on parameters
# mu
# muvals <- seq(from=0.1,to=2,length=30)
# uvals <- data.frame(mu = muvals)
# uvals$utility <- 0
# for (i in 1:30) uvals[i,]$utility <- ufun_Atta(Resources=1,hsize=1,params=list(theta=0.115,gamma=1.4,phival=0.9,mu=muvals[i],cutoff=0.1))
# plot(uvals$mu,uvals$utility,type="l")






# test kron
# =========

# 4 dimensions
library(ArmaUtils)
# make R data
a = matrix(c(1,2,3,4,5,0,0,6,0,0,0,7,8,0,9,0,0,0,10,0),nrow=5,ncol=4,byrow=T)
b <- matrix(c(1,2,0,0,0,3),nrow=2,ncol=3,byrow=T)
c <- matrix(c(0,1,0,3,4,0),nrow=3,ncol=2,byrow=T)
d <- matrix(c(0,1,2,0,2,0,0,1),nrow=2,ncol=4,byrow=T)
aa = as(a,"dgCMatrix")
bb = as(b,"dgCMatrix")
cc = as(c,"dgCMatrix")
dd = as(d,"dgCMatrix")
y <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)*ncol(dd)))

c.result <- krons(matrices=list(aa,bb,cc,dd),y=y)

# test with 4 dimensions. 
# ------------------------

# R kronecker product
r.result <- kronecker(a,kronecker(b,kronecker(c,d))) %*% y

# check output
context("testing krons")
test_that("4 dimensional kronecker product multiplying a vector is numeric",{
		  expect_that( is.numeric(c.result), is_true() )})
test_that("4dim krons is equal to R version",{
		  expect_that( all.equal(as.numeric(r.result),c.result), is_true() )})




# test with 3 dimensions
# ----------------------

y3 <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)))
c.result <- krons(matrices=list(aa,bb,cc),y=y3)
r.result <- kronecker(a,kronecker(b,c)) %*% y3
test_that("3dim krons is equal to R version",{
		  expect_that( all.equal(as.numeric(r.result),c.result), is_true() )})

# test with 2 dimensions
# ----------------------

y2 <- as.numeric(1:(ncol(aa)*ncol(bb)))
c.result <- krons(matrices=list(aa,bb),y=y2)
r.result <- kronecker(a,b) %*% y2
test_that("3dim krons is equal to R version",{
		  expect_that( all.equal(as.numeric(r.result),c.result), is_true() )})




# test ufun_labour_h
# ================

pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(1:4,5:-1)
w <- runif(4,min=1,max=4)
s <- c(0,0,1,2)
cc <- ufun_Atta_L(Res,w,s,pars)
source("../rtests.r")
rr <- ufun_labouR_h(Res,w,s,pars) 	# same function in R

context("testing ufun with labor supply, ufun_Atta_L")
test_that("ufun_Atta_L returns a list with three components",{
		  expect_that( is.list(cc) & length(cc)==3, is_true() )})
test_that("all components of ufun_Atta_L are numeric",{
		  expect_that( all(unlist(lapply(cc,is.numeric))), is_true() )})
test_that("c$utility is equal to R$utility",{
		  expect_that( all.equal(cc$utility,rr$utility), is_true() )})
test_that("c$consis equal to R$cons",{
		  expect_that( all.equal(cc$cons,rr$cons), is_true() )})
test_that("c$labo is equal to R$labo",{
		  expect_that( all.equal(cc$labo,rr$labo), is_true() )})
test_that("non-integer house size throws an error",{
		  expect_that( ufun_Atta_L(Res,w,0.2,pars), throws_error() )})
test_that("house size outside of c(0,1,2) throws an error",{
		  expect_that( ufun_Atta_L(Res,w,3,pars), throws_error() )})


# benchmark against ufun_Atta
# ===========================

# how much slower do we get if we also compute labour supply?

pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(seq(1,4,le=1000),seq(5,-1,le=100))
s    <- sample(0:2,size=nrow(Res),replace=TRUE)
w    <- runif(nrow(Res),min=1,max=4)

benchmark(atta = ufun_Atta(Res,s,pars), atta_L = ufun_Atta_L(Res,w,s,pars),replications=20)


