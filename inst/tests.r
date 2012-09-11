

# test file for ArmaUtils package

# armamax is a fast implementation of apply(X,1,max), that returns a list with 2 vectors 'values' and 'indices', containing maximal values and indices for each row.

A <- matrix(rnorm(20),nrow=4,ncol=5)


apply(A,1,max)




