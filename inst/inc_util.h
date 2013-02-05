
#ifndef _ArmaUtils_INC_UTIL_H
#define _ArmaUtils_INC_UTIL_H

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



#endif
