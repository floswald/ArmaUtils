#include "ufuns.h"
#include <iostream>

using namespace Rcpp ;
using namespace arma ;

// u(c,h,l,H,n,theta,age,wage)
SEXP utilfun( SEXP Res_, SEXP s_, SEXP age_, SEXP par_, SEXP xi_, SEXP own_ )
{
	mat Res = Rcpp::as<arma::mat>(Res_);
	vec siz = Rcpp::as<arma::vec>(s_);
	int age = Rcpp::as<int>(age_);	// age
	bool own = Rcpp::as<bool>(own_);	// owner yes/no

	// now get parameters out of lists par xi
	Rcpp::List par( par_ ) ;
	Rcpp::List xi( xi_ ) ;	

	double theta = Rcpp::as<double>(par["theta"]);
	double cutoff = Rcpp::as<double>(par["cutoff"]);
	double alpha  = Rcpp::as<double>(par["alpha"]);

	// derived data
	// equation numbering in doc/ufun.pdf
	vec xi1 = xi[0]; // equation (1)
	vec xi2 = xi[1]; // equation (2)
	vec xivec3 = xi[2];	// equation (3)
	vec xivec4 = xi[3];	// equation (4)
	vec xivec5 = xi[4];	// equation (5)
	// use as scalars
	double xi3 = as_scalar(xivec3);
	double xi4 = as_scalar(xivec4);
	double xi5 = as_scalar(xivec5);

	//setup output matrix
	mat ret(Res);
	int n = ret.n_rows;
	int m = ret.n_cols;
	ret.zeros();
	double diff, g, u, grad, hess;

	// test if house size is a vector of length Res.rows().
	if (siz.n_elem == 1){
		siz = repmat(siz,n,1);
	} else if (siz.n_elem > 1){
		if (siz.n_elem != n){
			Rcpp::Rcout << "ERROR: length of house sizes is not equal length of cash" << std::endl;
			return wrap(0);
		}
	}

	// begin calculations
	// if owner, adjust house size by premium
	if (own) {
		siz *= theta;
	}

	vec sizefac = pow(siz,xi4);
	rowvec tmpvec(m);

	for (int i=0; i<n; i++){
		if( Res(i,m-1) < cutoff ){	// check if the last entry is below cutoff. if first is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<m; j++){
				if( Res(i,j) < cutoff ){ // compute approximated utility at all entries of row i below cutoff
					g        = pow(cutoff,xi3) * sizefac(i);	  // equation (7)
					grad     = xi1(age-1) * alpha * g / cutoff;   // equation (8)
					hess     = xi5 * grad / cutoff;				  // equation (9)
					diff     = Res(i,j) - cutoff;
					ret(i,j) = xi2(age-1) * (g - 1) + (grad * diff) + 0.5 * hess * pow(diff,2);	// equation (10)
				} else {
					g = pow(Res(i,j),xi3) * sizefac(i);	    // equation (6)
					ret(i,j) = xi2(age-1) * (g - 1);	    // equation (11)
				}
			}
		} else {
			tmpvec = pow(Res.row(i),xi3) * sizefac(i);   // equation (6)
			ret.row(i) = xi2(age-1) * (tmpvec - 1.0);    // equation (11)
		}
	}
	return wrap(ret);
}


// u(c,h) as in Attansio et al. Review of economic dynamics 2012
// in: 
// resources matrix
// hsize vector
// pars list
//
// out:
// matrix
SEXP ufun_Attanasio( SEXP Res_, SEXP s_, SEXP par_)
{

	mat Res = Rcpp::as<arma::mat>(Res_);
	uvec hsize = Rcpp::as<arma::uvec>(s_);
	if (hsize.n_elem != Res.n_rows){
		Rcpp::Rcout << "ERROR: length of house sizes is not equal length of cash" << std::endl;
		return wrap(0);
	}

	//par_ = list(gamma,theta,phival,mu,cutoff)

	// now get parameters out of lists par 
	Rcpp::List par( par_ ) ;

	double gamma = Rcpp::as<double>(par["gamma"]);
	double theta = Rcpp::as<double>(par["theta"]);
	double cutoff = Rcpp::as<double>(par["cutoff"]);
	double phival = Rcpp::as<double>(par["phival"]);
	double mu = Rcpp::as<double>(par["mu"]);
	double mgamma = 1-gamma;
	double imgamma = 1/mgamma;

	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<hsize.size(); i++) {
		phivec(i) = phivals( hsize( i ) );
	}

	//setup output matrix
	mat ret(Res);
	int n = ret.n_rows;
	int m = ret.n_cols;
	ret.zeros();
	double diff, g, u, grad, hess;
	rowvec tmpvec(m);

	// compute utility from non-durable consumption
	for (int i=0; i<n; i++){
		if( Res(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<m; j++){
				if( Res(i,j) < cutoff ){ // compute approximated utility at all entries of row i below cutoff
					g        = pow(cutoff,mgamma);	 
					grad     = g / cutoff;   
					hess     = -gamma * grad / cutoff;	
					diff     = Res(i,j) - cutoff;
					ret(i,j) = imgamma*g + (grad * diff) + 0.5 * hess * pow(diff,2);
				} else {
					g = pow(Res(i,j),mgamma);	  
					ret(i,j) = imgamma * g;	    
				}
			}
		} else {
			tmpvec = pow(Res.row(i),mgamma);   // equation (6)
			ret.row(i) = imgamma * tmpvec;    // equation (11)
		}
	}
	// add additive premium if houseing is of right size
	mat phimat = repmat(phivec,1,m);
	mat hfac = exp( theta * phimat);
	ret = ret % hfac;
	ret = ret + mu * phimat;
	return wrap(ret);
}


