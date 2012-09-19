rm(list=ls(all=T))



# utilfun inline development and tests

library(inline)
library(RcppArmadillo)
library(RcppEigen)

cppcode <- '
#include <iostream>
using namespace std;
using namespace arma;

mat Res = Rcpp::as<arma::mat>(Res_);
vec siz = Rcpp::as<arma::vec>(siz_);
int age = Rcpp::as<int>(age_);	// age
int own = Rcpp::as<int>(own_);	// owner yes/no

// now get parameters out of list par and convert to Eigen types
Rcpp::List par( par_ ) ;
Rcpp::List xi( xi_ ) ;	

double theta = Rcpp::as<double>(par["theta"]);
double cutoff = Rcpp::as<double>(par["cutoff"]);
double alpha  = Rcpp::as<double>(par["alpha"]);
Rcpp::Rcout << "theta = " << std::endl;
Rcpp::Rcout << theta << std::endl << std::endl;
Rcpp::Rcout << "alpha = " << std::endl;
Rcpp::Rcout << alpha << std::endl << std::endl;
Rcpp::Rcout << "cutoff = " << std::endl;
Rcpp::Rcout << cutoff << std::endl << std::endl;

// derived data
// equation numbering in doc/ufun.pdf

Rcpp::Rcout << "Res = " << std::endl;
Rcpp::Rcout << Res << std::endl << std::endl;
Rcpp::Rcout << "siz = " << std::endl;
Rcpp::Rcout << siz << std::endl << std::endl;
Rcpp::Rcout << "age = " << std::endl;
Rcpp::Rcout << age << std::endl << std::endl;
Rcpp::Rcout << "own = " << std::endl;
Rcpp::Rcout << own << std::endl << std::endl;

vec xi1 = xi[0]; // equation (1)
vec xi2 = xi[1]; // equation (2)
vec xivec3 = xi[2];	// equation (3)
vec xivec4 = xi[3];	// equation (4)
vec xivec5 = xi[4];	// equation (5)

// use as scalars
double xi3 = as_scalar(xivec3);
double xi4 = as_scalar(xivec4);
double xi5 = as_scalar(xivec5);

Rcpp::Rcout << "xi1 = " << std::endl;
Rcpp::Rcout << xi1 << std::endl << std::endl;
Rcpp::Rcout << "xi2 = " << std::endl;
Rcpp::Rcout << xi2 << std::endl << std::endl;
Rcpp::Rcout << "xi3 = " << std::endl;
Rcpp::Rcout << xi3 << std::endl << std::endl;
Rcpp::Rcout << "xi4 = " << std::endl;
Rcpp::Rcout << xi4 << std::endl << std::endl;
Rcpp::Rcout << "xi5 = " << std::endl;
Rcpp::Rcout << xi5 << std::endl << std::endl;

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
		throw std::runtime_error( "length of house sizes is not equal length of cash" );
	}
}

// begin calculations
// if owner, adjust house size by premium
if (own==1) {
	siz *= theta;
}

vec sizefac = pow(siz,xi4);
rowvec tmpvec(m);

for (int i=0; i<n; i++){
	if( Res(i,m-1) < cutoff ){	// check if the last entry is below cutoff. if last is not, none is.
	Rcpp::Rcout << "Res(i,m) < cutoff" << std::endl;
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
'

ufun <- cxxfunction(signature(Res_="matrix",siz_="numeric",age_="integer",own_="integer",par_="list",xi_="list"),body=cppcode,plugin="RcppArmadillo")

# data
pars <- list(alpha=0.6,sigma=1.6,cutoff=0.05,theta=1.1)
maxage <- 10
fsize <- seq(1,3,le=maxage)
xi <- list()
xi[[1]]    <- fsize^(pars$sigma -1)
xi[[2]]    <- xi[[1]] / (1-pars$sigma)
xi[[3]]    <- pars$alpha*(1-pars$sigma)
xi[[4]]    <- (1-pars$alpha)*(1-pars$sigma)
xi[[5]]    <- pars$alpha*(1-pars$sigma) -1 

Res  <- outer(1:4,5:-1)
s <- array(c(1,2,3,4),dim=c(4,1))
age <- 3L
own <- 1L

# call
ufun(Res_=Res,siz_=s,age_=age,own_=own,par_=pars,xi_=xi)

Res  <- outer(1:4,-1:5)
ufun(Res_=Res,siz_=s,age_=age,own_=own,par_=par_,xi_=xi)
ufun(Res_=Res,siz_=2,age_=age,own_=own,par_=par_,xi_=xi)	# can handle house size vectors and scalars.


# KronProdSPMat4: c++ kronecker product for 4 dimensions. geared towards spline basis matrices.
# =============================================================================================

include <- '
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat4(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3,
		Eigen::VectorXd y) {

	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows() * a3.rows()  );
	if ( y.rows() != a0.cols() * a1.cols() * a2.cols()  * a3.cols()  ) {
		cout << "KronProdMat5 error: y and matrices not conformable" << endl;
	}

	//loop rows a0
	for (int row_idx0=0; row_idx0<a0.outerSize(); ++row_idx0) {
		int row_offset1 = row_idx0;
		row_offset1    *= a1.rows();

		// loop rows a1
		for (int row_idx1=0; row_idx1<a1.outerSize(); ++row_idx1) {
			int row_offset2 = row_offset1 + row_idx1;
			row_offset2    *= a2.rows();

			// loop rows a2
			for (int row_idx2=0; row_idx2<a2.outerSize(); ++row_idx2) {
				int row_offset3 = row_offset2 + row_idx2;
				row_offset3    *= a3.rows();

				// loop rows a3
				for (int row_idx3=0; row_idx3<a3.outerSize(); ++row_idx3) {

					// loop cols a0 (non-zero elements only)
					for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it0(a0,row_idx0); it0; ++it0) {
						int col_offset1 = it0.index();
						col_offset1    *= a1.innerSize();
						double factor1 = it0.value();

						// loop cols a1
						for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it1(a1,row_idx1); it1; ++it1) {
							int col_offset2 = col_offset1 + it1.index();
							col_offset2    *= a2.innerSize();
							double factor2  = factor1 * it1.value();

							//loop cols a2
							for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it2(a2,row_idx2); it2; ++it2) {
								int col_offset3 = col_offset2 + it2.index();
								col_offset3    *= a3.innerSize();
								double factor3  = factor2 * it2.value();

								for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it3(a3,row_idx3); it3; ++it3){
										retvec( row_offset3 + row_idx3 ) += factor3 * it3.value() * y( col_offset3 + it3.index() );
								}
							}
						}
					}
				}
			}
		}
	}
	return retvec;
}
'

cpp.code <- '
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;
const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
const MappedSparseMatrix<double> a3(as<MappedSparseMatrix<double> >(aa3));
const MapVecd y(as<MapVecd>(yy));

VectorXd result( KronProdSPMat4( a0,a1,a2,a3,y ) );
return wrap( result );
'

# create c function

kroncpp4 <- cxxfunction( signature(aa0="dgCMatrix",aa1="dgCMatrix",aa2="dgCMatrix",aa3="dgCMatrix",yy="numeric"), body=cpp.code, includes = include, plugin="RcppEigen")

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


# call
kroncpp4(aa,bb,cc,dd,y)

