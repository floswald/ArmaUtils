#ifndef _ArmaUtils_UFUNS_H
#define _ArmaUtils_UFUNS_H

#include <RcppArmadillo.h>
#include "handler.h"

using namespace arma;
using namespace Rcpp;

// quadradtic approximation function for case of no work
double approx(double cons,double x,double xi1,double alpha,double imgamma){
    signal(SIGSEGV, handler);   // install our handler
	double g,grad,hess,diff;
	g = pow( x, xi1);
	grad = alpha*g / x;
	hess = (xi1 -1) * grad / x;
	diff = cons - x;
	return imgamma * g + (grad * diff) + 0.5 * hess * pow(diff,2);
}

// quadratic approx for work
double approxw(double cons, double x,double alphaxi1,double malphaxi2, double mgamma, double imgamma, double wage, double xi2, double gamma){
    signal(SIGSEGV, handler);   // install our handler
	double g,grad,hess,diff;
	g = alphaxi1 * malphaxi2 * pow( x, mgamma) / pow( wage, xi2) ;
	grad = g / x;
	hess = - gamma * grad / x;
	diff = cons - x;
	return imgamma * g + (grad * diff) + 0.5 * hess * pow(diff,2);
}



// no work, positive resources
vec u_no_pos(vec e, double xi, double imgamma){
    signal(SIGSEGV, handler);   // install our handler
	vec g = pow( e, xi);
	return imgamma * g;
}

// work, positive resources
vec u_work_pos(double x1, double x2, double xi2, vec e, double mgamma, vec w){
    signal(SIGSEGV, handler);   // install our handler
	vec g = x1 * x2 * pow( e, mgamma) / pow( w, xi2) ;
	return (1/mgamma) * g;
}

// quadratic approximation function for case of no work for vector
vec u_no_neg(vec cons,double x,double alpha,double gamma){
    signal(SIGSEGV, handler);   // install our handler
	vec diff;
	double g,grad,hess;
	double xi1 = alpha * (1-gamma);
    double imgamma = 1/(1-gamma);
	g = pow( x, xi1);
	grad = alpha*g / x;
	hess = (xi1 -1) * grad / x;
	diff = cons - x;
	return imgamma * g + (grad * diff) + 0.5 * hess * pow(diff,2);
}

// quadratic approximation function for case of no work for scalar
vec u_work_neg(vec cons,double x,double alphaxi1,double malphaxi2, double xi2, double gamma,vec w){
    signal(SIGSEGV, handler);   // install our handler
	vec diff,g,grad,hess;
	double mgamma = 1-gamma;
	double imgamma = 1/mgamma;
	g = alphaxi1 * malphaxi2 * pow( x, mgamma) / pow( w, xi2);
	grad = g / x;
	hess = (-gamma/x) * grad;
	diff = cons - x;
	return imgamma * g + (grad % diff) + 0.5 * hess % pow(diff,2);
}

// alternative implementation of the function in ufuns.cpp
Rcpp::List ufun_Atta_L2(mat Res, Rcpp::List par, uvec hsize, vec wage){
    signal(SIGSEGV, handler);   // install our handler

	int n = Res.n_rows;
	int m = Res.n_cols;
	int nh = hsize.n_elem;

	// input tests
	if (hsize.n_elem != Res.n_rows){
throw std::logic_error( "ufun_Attanasio: hsize and Res are not equal rows!" );
		return R_NilValue;
	}

	if ( wage.n_elem != n){
		throw std::logic_error("error. wage and Res are not conformable");
		return R_NilValue;
	}
	if ( hsize.n_elem != Res.n_rows){
		throw std::logic_error("error. hsize and Res are not conformable");
		return R_NilValue;
		}

	uvec neg = wage <= 0;
	if ( sum(neg) > 0 ){
		throw std::logic_error("error. wage must be positive");
		return R_NilValue;
	}

	// extract elements from list par and make some paramters
	double theta     = Rcpp::as<double>(par["theta"]);
	double phival    = Rcpp::as<double>(par["phival"]);
	double mu        = Rcpp::as<double>(par["mu"]);
	double gamma     = Rcpp::as<double>(par["gamma"]);
	double cutoff    = Rcpp::as<double>(par["cutoff"]);
	double alpha     = Rcpp::as<double>(par["alpha"]);
	double xi1       = alpha * (1-gamma);
	double xi2       = (1-alpha)*(1-gamma);
	double malpha    = 1-alpha;
	double mgamma    = 1-gamma;
	double imgamma   = 1/mgamma;
	double alphaxi1  = pow(alpha, xi1);
	double malphaxi2 = pow(malpha, xi2);

	// make vector that holds the right phi value
	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<nh; i++) {
		phivec(i) = phivals( hsize( i ) );
	}

	// repmat wage into wagemat
	mat wagemat = repmat(wage,1,Res.n_cols);
	mat consmat = Res - wagemat;

	// initiate return objects as zero
	mat util(Res), cons(Res), labour(Res);
	util.zeros();
	cons.zeros();
	labour.zeros();


	// split Res according to cases:
	// work
	//		- pos resouces
	//		- neg resources
	// no work
	//		- pos resouces
	//		- neg resources
	uvec iwork    = find(malpha*Res < wagemat);
	uvec inowork  = find(malpha*Res >= wagemat);
	vec workres   = Res.elem(iwork);
	vec noworkres = consmat.elem(inowork);
	vec workwage  = wagemat.elem(iwork);
	// collect results in thress vectors each
	vec uwork(workres), cwork(workres), lwork(workres);
	vec unowork(noworkres), cnowork(noworkres), lnowork(noworkres);

	// resources in each subcase
	// 1) work. pos and neg resources
	uvec iworkpos   = find( workres > cutoff );
	uvec iworkneg   = find( workres <= cutoff );
	vec workposres  = workres.elem( iworkpos );
	vec worknegres  = workres.elem( iworkneg );
	vec workposwage = workwage.elem( iworkpos );
	vec worknegwage = workwage.elem( iworkneg );
	// 2) NO work. pos and neg resources
	uvec inoworkpos  = find( noworkres > cutoff );
	uvec inoworkneg  = find( noworkres <= cutoff );
	vec noworkposres = noworkres.elem( inoworkpos );
	vec noworknegres = noworkres.elem( inoworkneg );

	// calculate utility in each case and
	// reassemble vectors from pos/neg
	// work, positive consupmtion
	uwork.elem(iworkpos) = u_work_pos(alphaxi1,malphaxi2,xi2,workposres,mgamma,workposwage);
	cwork.elem(iworkpos) = alpha * workposres;
	lwork.elem(iworkpos) = 1 - malpha * workposres / workposwage;
	if (!(iworkneg.is_empty())) {
	// work, negative consupmtion
		uwork.elem(iworkneg) = u_work_neg(worknegres,cutoff,alphaxi1,malphaxi2,xi2,gamma,worknegwage);
		cwork.elem(iworkneg) = worknegres;
		lwork.elem(iworkneg) = zeros<vec>(worknegres.n_elem);
	}

	// no work cases
	// positive
	unowork.elem(inoworkpos) = u_no_pos(noworkposres,xi1,imgamma);
	cnowork.elem(inoworkpos) = noworkposres;
	if (!(inoworkneg.is_empty())){
		// negative
		unowork.elem(inoworkneg) = u_no_neg(noworknegres,cutoff,alpha,gamma);
		cnowork.elem(inoworkneg) = noworknegres;
	}
	
	//reassemble from work/nowork
	util.elem(iwork)   = uwork;
	util.elem(inowork) = unowork;
	cons.elem(iwork)   = cwork;
	cons.elem(inowork) = cnowork;
	labour.elem(iwork) = lwork;

	// add additive premium for houseing
	mat phimat = repmat(phivec,1,m);
	mat hfac   = exp( theta * phimat);
	util       = util % hfac;
	util       = util + mu * phimat;
	
	Rcpp::List rlist = Rcpp::List::create( _["utility"] = util, _["consumption"] = cons , _["labour"] = labour);
	return rlist;
}


// utility with positive and negative consumption
vec u_pos(vec c,vec lab,vec h,double alph,double gamm){
    signal(SIGSEGV, handler);   // install our handler
	double mgamm = 1-gamm;
	double imgamm = 1/mgamm;
	double g,z,y;
	vec ret(c);
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab[i] ), mgamm);
		z      = y * h[i];
		g      = pow( c[i], mgamm);
		ret[i] = imgamm * g * z;
	}
	return ret;
}

// quadratic approximation function for case of no work for vector
vec u_neg(vec c,double cuto,vec lab,vec h,double alph,double gamm){
    signal(SIGSEGV, handler);   // install our handler
	vec diff;
	vec ret(c);
	double grad,hess,y,z;
	double mgamm = 1-gamm;
	double imgamm = 1/mgamm;
	diff = c - cuto;
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab[i] ), mgamm);
		z      = y * h[i];
		grad   = pow( cuto, -gamm) * z;
		hess   = -(gamm * grad) / cuto;
		ret[i] = imgamm * pow( cuto, 1-gamm) * z  + (grad * diff[i]) + 0.5 * hess * pow(diff[i],2);
	}
	return ret;
}


//' CRRA utility augmented with discrete labor supply and utility from housing
//'
//' same as other util functions but with discrete labor supply choice. no time separability, i.e. labor supply is not implied.
//' @param Resources numeric matrix of consumption levels
//' @param hsize numeric vector of house sizes
//' @param labor numeric vector of discrete labor supply choices
//' @param params list of parameters 
//' \itemize{
//' \item{theta}{elasticity of substitution between c and h}
//' \item{phival}{value of relative utility difference flat vs house}
//' \item{mu}{weight on additive utility premium}
//' \item{gamma}{coefficient of relative risk aversion}
//' \item{cutoff}{minimum level of consumption. below cutoff, u(c) is quadratically approximated.}
//' \item{alpha}{coefficient on labor}
//' }
//' @return numeric matrix of utility values 
// [[Rcpp::Export]]
SEXP ufun_discreteL(mat Res, Rcpp::List par, uvec hsize, vec labor){
    signal(SIGSEGV, handler);   // install our handler

	BEGIN_RCPP

	uword n = Res.n_rows;
	uword m = Res.n_cols;

    if (labor.size() != n){
		throw std::logic_error("ufun_discreteL:::error. labor vector not equal rows of cons matrix.");
		return R_NilValue;
	}

	if ( hsize.n_elem != Res.n_rows){
		throw std::logic_error("ufun_discreteL:::error. hsize and Res are not conformable");
		return R_NilValue;
	}

	// extract elements from list par and make some paramters
	double theta     = Rcpp::as<double>(par["theta"]);
	double phival    = Rcpp::as<double>(par["phival"]);
	double mu        = Rcpp::as<double>(par["mu"]);
	double gamma     = Rcpp::as<double>(par["gamma"]);
	double cutoff    = Rcpp::as<double>(par["cutoff"]);
	double alpha     = Rcpp::as<double>(par["alpha"]);

	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<hsize.n_elem; i++) {
		phivec(i) = phivals( hsize( i ) );
	}

	// initiate return objects as zero
	mat util(Res);
	util.zeros();

	// prepare additive housing premium
	mat phimat = repmat(phivec,1,m);
	mat hfac = exp( theta * phimat);

	mat labmat = repmat(labor,1,m);

	// split Res according to cases:
	//		- pos resouces
	//		- neg resources
	uvec ipos  = find( Res >= cutoff );
	uvec ineg  = find( Res < cutoff );
	vec posres = Res.elem( ipos );
	vec negres = Res.elem( ineg );
	vec posh   = hfac.elem( ipos );
	vec negh   = hfac.elem( ineg );
	vec posl   = labmat.elem( ipos );
	vec negl   = labmat.elem( ineg );

	// calculate utility in each case
	vec upos = u_pos(posres,posl,posh,alpha,gamma);
	vec uneg = u_neg(negres,cutoff,negl,negh,alpha,gamma);

	// reassemble vectors from pos/neg
	util.elem(ipos) = upos;

	if (!(ineg.is_empty())){
		util.elem(ineg) = uneg;
	}

	// add additive premium if houseing is of right size
	util = util + mu * phimat;

	return wrap(util);

	END_RCPP
}

// Rcpp Exporting of functions
RcppExport SEXP utilfun(  SEXP Res_, SEXP s_, SEXP age_, SEXP par_, SEXP xi_, SEXP own_ ) ; 
RcppExport SEXP ufun_Attanasio( SEXP Res_, SEXP s_, SEXP par_) ; 
RcppExport SEXP ufun_Attanasio_L( SEXP Res_, SEXP w_, SEXP s_, SEXP par_) ; 
RcppExport SEXP ufun_Attanasio_L2( SEXP Res_, SEXP w_, SEXP s_, SEXP par_) ; 
RcppExport SEXP ufun_discLabor( SEXP Res_, SEXP l_,SEXP s_, SEXP par_);

#endif
