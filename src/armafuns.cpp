#include "armafuns.h"

using namespace Rcpp ;
using namespace arma ;

SEXP armamax( SEXP X_){
	// map R objects
	arma::mat A  = Rcpp::as<arma::mat>(X_);	// map R matrix

	// allocate memory
	arma::uword j;
	arma::uvec iy(A.n_rows);
	arma::vec y(A.n_rows);
	arma::rowvec tmp(A.n_cols);

	// loop over rows of A and find maximal value and it's index for each row.
	for (int i=0; i<A.n_rows; i++) {
		tmp = A.row(i);
		y(i) = tmp.max( j );
		iy(i) = j + 1; // go to 1 based indices
	}
	Rcpp::List list = Rcpp::List::create( _["values"] = y, _["indices"] = iy );
	return list;
}

SEXP armamax2( SEXP X_, SEXP B_ ){
	// map R objects
	arma::mat A  = Rcpp::as<arma::mat>(X_);	// map R matrix
	arma::uvec b = Rcpp::as<arma::uvec>(B_);	// map R index vector

	// allocate memory
	arma::uword j;
	arma::uvec iy(A.n_rows);
	arma::vec y(A.n_rows);
	arma::rowvec tmp(A.n_cols);
	arma::rowvec tmp2;
	int m = A.n_cols - 1;

	// loop over rows of A and find maximal value and it's index for each row.
	for (int i=0; i<A.n_rows; i++) {
		if ( b( i ) == 0 ) {
			// no borrowing limit. max over entire row.
			tmp = A.row(i);
			y(i) = tmp.max( j );
			iy(i) = j + 1; // go to 1 based indices
		} else {
			tmp2 = A( i , arma::span( b(i), m ) );
			y(i) = tmp2.max( j );
			iy(i) = j + b(i) + 1; // go to 1 based indices
		}
	}
	Rcpp::List list = Rcpp::List::create( _["values"] = y, _["indices"] = iy );
	return list;
}

SEXP armakron( SEXP y_, SEXP R_ ){
	// the following statement enables to pass error messages from cpp to R without a crash.
	
	BEGIN_RCPP

	// we get a vector of values and a list of matrices from R
	vec y = Rcpp::as<vec>(y_)	;
	Rcpp::List mats(R_);
	int nmat =  mats.length();	// how many matrices
	int nall = y.size();
	vec cols(nmat);

	// need to map to arma::mat type
	std::vector<mat> list;
	for (size_t i=0; i<nmat; i++){
		list.push_back( Rcpp::as<mat>( mats[ i ] ) );
		// check matrices are square
		if ( !list.at(i).is_square() ){
			throw std::logic_error( "armakron: not all matrices are square!" );
			return R_NilValue;
		}
		// fill in dimensions
		cols(i) = list.at(i).n_cols;
	}

	if ( prod(cols) != nall ){
		throw std::logic_error( "armakron: prod(cols(matrices)) not equal length(y)" );
		return R_NilValue;
	}

	// TODO do a sanity check on dimensions of objects!!

	// product for first matrix
	vec y0 = y;
	vec y1(nall);
	y1.zeros();
	int n = list.at(0).n_rows;
	int m = nall/n;
	uvec lhs_out = linspace<uvec>(0,n-1,n);

	uvec lhs_in = lhs_out * m;
	uvec rhs(n);
	uvec lhs(n);

	for (int i=0; i<m; i++){
		lhs = lhs_in + i;
		rhs    = lhs_out + i*n;
		y1.elem( lhs ) = list.at( 0 ) * y0.elem( rhs ) ;
	}

	// process all other matrices
	if (nmat > 1){
		for(int imat=1; imat < nmat; imat++){
			y0 = y1;
			n  = list.at(imat).n_rows;
			m  = nall/n;
			lhs_out.resize(n);
			lhs_out = linspace<uvec>(0,n-1,n);
			uvec lhs_in = lhs_out * m;
			rhs.resize(n);
			lhs.resize(n);
			for (int i=0; i<m; i++){
				lhs = lhs_in + i;
				rhs    = lhs_out + i*n;
				y1.elem( lhs ) = list.at( imat ) * y0.elem( rhs ) ;
			}
		}
	}
	return Rcpp::wrap(y1);

	END_RCPP
		// turn off RCPP error handling
}

