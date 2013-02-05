#ifndef _RcppUtils_UFUNS_H
#define _RcppUtils_UFUNS_H

#include <RcppArmadillo.h>

// quadradtic approximation function for case of no work
double approx(double cons,double x,double xi1,double alpha,double imgamma){
	double g,grad,hess,diff;
	g = pow( x, xi1);
	grad = alpha*g / x;
	hess = (xi1 -1) * grad / x;
	diff = cons - x;
	return imgamma * g + (grad * diff) + 0.5 * hess * pow(diff,2);
}

// quadratic approx for work
double approxw(double cons, double x,double alphaxi1,double malphaxi2, double mgamma, double imgamma, double wage, double xi2, double gamma){
	double g,grad,hess,diff;
	g = alphaxi1 * malphaxi2 * pow( x, mgamma) / pow( wage, xi2) ;
	grad = g / x;
	hess = - gamma * grad / x;
	diff = cons - x;
	return imgamma * g + (grad * diff) + 0.5 * hess * pow(diff,2);
}




// Rcpp Exporting of functions
RcppExport SEXP utilfun(  SEXP Res_, SEXP s_, SEXP age_, SEXP par_, SEXP xi_, SEXP own_ ) ; 
RcppExport SEXP ufun_Attanasio( SEXP Res_, SEXP s_, SEXP par_) ; 
RcppExport SEXP ufun_Attanasio_L( SEXP Res_, SEXP w_, SEXP s_, SEXP par_) ; 

#endif
