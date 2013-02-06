rm(list=ls(all=T))

setwd("~/git/ArmaUtils/inst")


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



# ufun_Atta

cpp.code <- '
#include <iostream>
using namespace std;
using namespace arma;

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

Rcpp::Rcout<< " Res " << Res << endl;
Rcpp::Rcout<< " hsize" << hsize << endl;
Rcpp::Rcout<< " phivec " << phivec << endl;


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
'

ufun_Atta <- cxxfunction(signature(Res_="matrix",s_="integer",par_="list"),body=cpp.code,plugin="RcppArmadillo")
pars <- list(theta=0.115,gamma=1.4,phival=0.9,mu=0.26,cutoff=0.1)
Res  <- outer(1:4,5:-1)
s <- array(c(1,2,0,0),dim=c(4,1))
ufun_Atta(Res,s,pars)

# ufun_labour 
# ===========

cpp.code <- '
#include <iostream>
using namespace std;
using namespace arma;

mat Res = Rcpp::as<arma::mat>(Res_);
int n = Res.n_rows;
int m = Res.n_cols;
vec wage = Rcpp::as<arma::vec>(w_);

if ( wage.n_elem != n){
	Rcpp::Rcout << "error. wage and Res are not conformable" << endl;
	return R_NilValue;
}

uvec neg = wage <= 0;
if ( sum(neg) > 0 ){
	Rcpp::Rcout << "error. wage must be positive" << endl;
	return R_NilValue;
}

//par_ = list(alpha,xi1,xi2,gamma,cutoff)

// now get parameters out of lists par 
Rcpp::List par( par_ ) ;

double gamma = Rcpp::as<double>(par["gamma"]);
double cutoff = Rcpp::as<double>(par["cutoff"]);
double alpha = Rcpp::as<double>(par["alpha"]);
	double xi1 = alpha * (1-gamma);
	double xi2 = (1-alpha)*(1-gamma);
	double malpha = 1-alpha;
	double mgamma = 1-gamma;
	double imgamma = 1/mgamma;
	double alphaxi1 = pow(alpha, xi1);
	double malphaxi2 = pow(malpha, xi2);

Rcpp::Rcout<< " Res " << Res << endl;
Rcpp::Rcout<< " wage " << wage << endl;
Rcpp::Rcout<< " alpha " << alpha << endl;
Rcpp::Rcout<< " malpha " << malpha << endl;
Rcpp::Rcout<< " alphaxi1 " << alphaxi1 << endl;
Rcpp::Rcout<< " malphaxi2 " << malphaxi2 << endl;
Rcpp::Rcout<< " xi1 " << xi1<< endl;
Rcpp::Rcout<< " xi2 " << xi2<< endl;
Rcpp::Rcout<< " gamma " << gamma << endl;

// repmat wage into wagemat
mat wagemat = repmat(wage,1,m);
mat idmat = malpha * (Res / wagemat);

Rcpp::Rcout<< " wagemat" << wagemat<< endl;
Rcpp::Rcout<< " idmat" << idmat<< endl;

// return matrices
mat retu = zeros<mat>(Res.n_rows,Res.n_cols);
mat retc = retu;
mat retn = retu;

//compute utilty. flow depends on
//a) whether mat work ==1 or ==0 and
//b) whether cons < 0

mat consmat = Res - wagemat; 	// resources in case of no work
double g;
rowvec tmpvec(m);

for( int i=0; i<n; i++){
	// in each row. check first work index
	uvec k = find(Res.row(i)*malpha < wagemat.row(i),1);
	cout << "k = " << k << endl;
	if( k.n_elem == 0 ){
		cout << "computing no work." << endl;
		//no work in this row
		if( consmat(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<m; j++){	// iterate over columns
				if( consmat(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approx(consmat(i,j),cutoff,xi1,alpha,imgamma);
					retc(i,j) = consmat(i,j);
					retn(i,j) = 0;
				} else {
					// entry j is not
					g = pow( consmat(i,j), xi1);
					retu(i,j) = imgamma * g;
					retc(i,j) = consmat(i,j);
				}
			}
		} else {
			tmpvec = pow(consmat.row(i),xi1);
			retu.row(i) = imgamma * tmpvec;
			retc.row(i) = consmat.row(i);
		}
	} else {
		// work starting at index k
		int ki = conv_to<int>::from(k);
		if( consmat(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<ki; j++){	// iterate over non work-columns
				// NO WORK
				if( consmat(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approx(consmat(i,j),cutoff,xi1,alpha,imgamma);
					retc(i,j) = consmat(i,j);
				} else {
					// entry j is not
					g = pow( consmat(i,j), xi1);
					retu(i,j) = imgamma * g;
					retc(i,j) = consmat(i,j);
				}
			}
			for (int j=ki; j<m; j++){  // iterate over work columns
				// WORK
				if( Res(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approxw(Res(i,j),cutoff,alphaxi1,malphaxi2, mgamma,imgamma, wagemat(i,j), xi2, gamma);
					retc(i,j) = Res(i,j) * alpha;
					retn(i,j) = 1- idmat(i,j);
				} else {
					// entry j is not
					g = alphaxi1 * malphaxi2 * pow( Res(i,j), mgamma) / pow( wagemat(i,j), xi2) ;
					retu(i,j) = imgamma * g;
					retc(i,j) = Res(i,j) * alpha;
					retn(i,j) = 1- idmat(i,j);
				}
			}
		} else {
			// no negative consumption in either work/non-work world in this row
			rowvec noworku = pow(consmat.row(i).cols(0,ki-1),xi1);
			rowvec worku = alphaxi1 * malphaxi2 * pow(Res.row(i).cols(ki,m),mgamma) / pow( wagemat.row(i).cols(ki,m), xi2);
			rowvec noworkc = consmat.row(i).cols(0,ki-1);
			rowvec workc = Res.row(i).cols(ki,m) * alpha;
			rowvec workn = 1 - idmat.row(i).cols(ki,m);
			retu.row(i) = join_cols(noworku,worku);
			retc.row(i) = join_cols(noworkc,workc);
			retn.row(i).cols(ki,m) = workn;
		}
	}
}
Rcpp::List rlist = Rcpp::List::create( _["utility"] = retu, _["consumption"] = retc , _["labour"] = retn);
return rlist;
'

ufun_labour <- cxxfunction(signature(Res_="matrix",w_="numeric",par_="list"),body=cpp.code,plugin="RcppArmadillo",includes="#include </Users/florianoswald/git/ArmaUtils/inst/inc_util.h>",verbose=TRUE)
pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1)
Res  <- outer(1:4,5:-1)
w <- runif(4,min=1,max=4)
cc <- ufun_labour(Res,w,pars)
source("rtests.r")
rr <- ufun_labouR(Res,w,pars) 	# same function in R
all.equal(cc$utility,rr$utility)	# TRUE
all.equal(cc$cons,rr$cons)	# TRUE
all.equal(cc$labo,rr$labo)	# TRUE


# ufun_labour_h
# ===========

cpp.code <- '
#include <iostream>
using namespace std;
using namespace arma;

// BEGIN_RCPP inline adds that

mat Res = Rcpp::as<arma::mat>(Res_);
int n = Res.n_rows;
int m = Res.n_cols;
vec wage = Rcpp::as<arma::vec>(w_);

uvec hsize = Rcpp::as<arma::uvec>(s_);
int nh = hsize.n_elem;
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

//par_ = list(alpha,xi1,xi2,gamma,cutoff)

// now get parameters out of lists par 
Rcpp::List par( par_ ) ;

double theta = Rcpp::as<double>(par["theta"]);
double phival = Rcpp::as<double>(par["phival"]);
double mu = Rcpp::as<double>(par["mu"]);
double gamma = Rcpp::as<double>(par["gamma"]);
double cutoff = Rcpp::as<double>(par["cutoff"]);
double alpha = Rcpp::as<double>(par["alpha"]);
double xi1 = alpha * (1-gamma);
double xi2 = (1-alpha)*(1-gamma);
double malpha = 1-alpha;
double mgamma = 1-gamma;
double imgamma = 1/mgamma;
double alphaxi1 = pow(alpha, xi1);
double malphaxi2 = pow(malpha, xi2);

Rcpp::Rcout<< " Res " << Res << endl;
Rcpp::Rcout<< " wage " << wage << endl;
Rcpp::Rcout<< " alpha " << alpha << endl;
Rcpp::Rcout<< " malpha " << malpha << endl;
Rcpp::Rcout<< " alphaxi1 " << alphaxi1 << endl;
Rcpp::Rcout<< " malphaxi2 " << malphaxi2 << endl;
Rcpp::Rcout<< " xi1 " << xi1<< endl;
Rcpp::Rcout<< " xi2 " << xi2<< endl;
Rcpp::Rcout<< " gamma " << gamma << endl;

vec phivals;
phivals << 0 << phival << 1 << endr;
vec phivec(hsize.size());
for (int i=0; i<nh; i++) {
	phivec(i) = phivals( hsize( i ) );
}


// repmat wage into wagemat
mat wagemat = repmat(wage,1,Res.n_cols);
mat idmat = malpha * Res / wagemat;

cout<< " wagemat" << wagemat<< endl;
cout<< " idmat" << idmat<< endl;

//labour indicator
//	umat work = find(idmat < 1);
//	cout<< " work" << work << endl;

//return matrices
mat retu = zeros<mat>(Res.n_rows,Res.n_cols);
mat retc = retu;
mat retn = retu;

//compute utilty. flow depends on
//a) whether mat work ==1 or ==0 and
//b) whether cons < 0

mat consmat = Res - wagemat; 	// resources in case of no work
cout << "consmat" << consmat << endl;
double g;
rowvec tmpvec(m);

for( int i=0; i<n; i++){
	// in each row. check first work index
	uvec k = find(Res.row(i)*malpha < wagemat.row(i),1);
	cout << "k = " << k << endl;
	if( k.n_elem == 0 ){
		cout << "computing no work." << endl;
		//no work in this row
		if( consmat(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<m; j++){	// iterate over columns
				if( consmat(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approx(consmat(i,j),cutoff,xi1,alpha,imgamma);
					retc(i,j) = consmat(i,j);
					retn(i,j) = 0;
				} else {
					// entry j is not
					g = pow( consmat(i,j), xi1);
					retu(i,j) = imgamma * g;
					retc(i,j) = consmat(i,j);
				}
			}
		} else {
			tmpvec = pow(consmat.row(i),xi1);
			retu.row(i) = imgamma * tmpvec;
			retc.row(i) = consmat.row(i);
		}
	} else {
		// work starting at index k
		int ki = conv_to<int>::from(k);
		if( consmat(i,m-1) < cutoff ){	// check if the last entry in row i is below cutoff. if last is not, none is. Res = c - savings, savings increase to the right.
			for (int j=0; j<ki; j++){	// iterate over non work-columns
				// NO WORK
				if( consmat(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approx(consmat(i,j),cutoff,xi1,alpha,imgamma);
					retc(i,j) = consmat(i,j);
				} else {
					// entry j is not
					g = pow( consmat(i,j), xi1);
					retu(i,j) = imgamma * g;
					retc(i,j) = consmat(i,j);
				}
			}
			for (int j=ki; j<m; j++){  // iterate over work columns
				// WORK
				if( Res(i,j) < cutoff ){
					// entry j is below cutoff
					retu(i,j) = approxw(Res(i,j),cutoff,alphaxi1,malphaxi2, mgamma,imgamma, wagemat(i,j), xi2, gamma);
					retc(i,j) = Res(i,j) * alpha;   // that is not the correct approx.
					retn(i,j) = 1- idmat(i,j);      // neither is this. but those neg cons values are not useful anyway.
				} else {
					// entry j is not
					g = alphaxi1 * malphaxi2 * pow( Res(i,j), mgamma) / pow( wagemat(i,j), xi2) ;
					retu(i,j) = imgamma * g;
					retc(i,j) = Res(i,j) * alpha;
					retn(i,j) = 1- idmat(i,j);
				}
			}
		} else {
			// no negative consumption in either work/non-work world in this row
			rowvec noworku = pow(consmat.row(i).cols(0,ki-1),xi1);
			rowvec worku = alphaxi1 * malphaxi2 * pow(Res.row(i).cols(ki,m),mgamma) / pow( wagemat.row(i).cols(ki,m), xi2);
			rowvec noworkc = consmat.row(i).cols(0,ki-1);
			rowvec workc = Res.row(i).cols(ki,m) * alpha;
			rowvec workn = 1 - idmat.row(i).cols(ki,m);
			retu.row(i) = join_cols(noworku,worku);
			retc.row(i) = join_cols(noworkc,workc);
			retn.row(i).cols(ki,m) = workn;
		}
	}
}
// add additive premium if houseing is of right size
mat phimat = repmat(phivec,1,m);
mat hfac = exp( theta * phimat);
retu = retu % hfac;
retu = retu + mu * phimat;

Rcpp::Rcout << "utility " << retu << endl;
Rcpp::Rcout << "consumption " << retc << endl;
Rcpp::Rcout << "labour " << retn << endl;
Rcpp::List rlist = Rcpp::List::create( _["utility"] = retu, _["consumption"] = retc , _["labour"] = retn);
return rlist;
'

ufun_labour_h <- cxxfunction(signature(Res_="matrix",w_="numeric",s_="numeric",par_="list"),body=cpp.code,plugin="RcppArmadillo",includes="#include </Users/florianoswald/git/ArmaUtils/inst/inc_util.h>")
pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(1:4,5:-1)
w <- runif(4,min=1,max=4)
s <- c(0,0,1,2)
cc <- ufun_labour_h(Res,w,s,pars)
rr <- ufun_labouR_h(Res,w,s,pars) 	# same function in R
all.equal(cc$utility,rr$utility)	# TRUE
all.equal(cc$cons,rr$cons)	# TRUE
all.equal(cc$labo,rr$labo)	# TRUE


# ufun_labour2
# ============


cpp.code <- '
#include <iostream>
using namespace std;
using namespace arma;

// BEGIN_RCPP inline adds that

mat Res = Rcpp::as<arma::mat>(Res_);
int n = Res.n_rows;
int m = Res.n_cols;
vec wage = Rcpp::as<arma::vec>(w_);

uvec hsize = Rcpp::as<arma::uvec>(s_);
int nh = hsize.n_elem;
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

//par_ = list(alpha,xi1,xi2,gamma,cutoff)

// now get parameters out of lists par 
Rcpp::List par( par_ ) ;

double theta = Rcpp::as<double>(par["theta"]);
double phival = Rcpp::as<double>(par["phival"]);
double mu = Rcpp::as<double>(par["mu"]);
double gamma = Rcpp::as<double>(par["gamma"]);
double cutoff = Rcpp::as<double>(par["cutoff"]);
double alpha = Rcpp::as<double>(par["alpha"]);
double xi1 = alpha * (1-gamma);
double xi2 = (1-alpha)*(1-gamma);
double malpha = 1-alpha;
double mgamma = 1-gamma;
double imgamma = 1/mgamma;
double alphaxi1 = pow(alpha, xi1);
double malphaxi2 = pow(malpha, xi2);

Rcpp::Rcout<< " Res " << Res << endl;
Rcpp::Rcout<< " wage " << wage << endl;
Rcpp::Rcout<< " alpha " << alpha << endl;
Rcpp::Rcout<< " malpha " << malpha << endl;
Rcpp::Rcout<< " alphaxi1 " << alphaxi1 << endl;
Rcpp::Rcout<< " malphaxi2 " << malphaxi2 << endl;
Rcpp::Rcout<< " xi1 " << xi1<< endl;
Rcpp::Rcout<< " xi2 " << xi2<< endl;
Rcpp::Rcout<< " gamma " << gamma << endl;

vec phivals;
phivals << 0 << phival << 1 << endr;
vec phivec(hsize.size());
for (int i=0; i<nh; i++) {
	phivec(i) = phivals( hsize( i ) );
}


	// repmat wage into wagemat
	mat wagemat = repmat(wage,1,Res.n_cols);
	mat consmat = Res - wagemat;
	mat idmat = malpha * Res / wagemat;

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
	uvec iwork = find(malpha*Res < wagemat);
	uvec inowork = find(malpha*Res >= wagemat);
	vec workres = Res.elem(iwork);
	vec noworkres = consmat.elem(inowork);
	vec workwage = wagemat.elem(iwork);
	// collect results in thress vectors each
	vec uwork(workres), cwork(workres), lwork(workres);
	vec unowork(noworkres), cnowork(noworkres), lnowork(noworkres);

	// resources in each subcase
	// 1) work. pos and neg resources
	uvec iworkpos = find( workres >= cutoff );
	uvec iworkneg = find( workres < cutoff );
	vec workposres = workres.elem( iworkpos );
	vec worknegres = workres.elem( iworkneg );
	vec workposwage = workwage.elem( iworkpos );
	vec worknegwage = workwage.elem( iworkneg );
	// 2) NO work. pos and neg resources
	uvec inoworkpos = find( noworkres >= cutoff );
	uvec inoworkneg = find( noworkres < cutoff );
	vec noworkposres = noworkres.elem( inoworkpos );
	vec noworknegres = noworkres.elem( inoworkneg );
	vec uworkneg  = zeros<vec>(iworkneg.n_elem);
	vec cworkneg  = zeros<vec>(iworkneg.n_elem);
	vec lworkneg  = zeros<vec>(iworkneg.n_elem);
	vec unoworkneg= zeros<vec>(inoworkneg.n_elem);
	vec cnoworkneg= zeros<vec>(inoworkneg.n_elem);

//	cout << "noworknegres " << noworknegres << endl;
//	cout << "iworkneg " << iworkneg << endl;
//	cout << "iworkpos " << iworkpos << endl;
//	cout << "inoworkneg " << inoworkneg << endl;
//	cout << "inoworkpos " << inoworkpos << endl;


	// calculate utility in each case
	vec uworkpos = u_work_pos(alphaxi1,malphaxi2,xi2,workposres,mgamma,workposwage);
	vec cworkpos = alpha * workposres;
	vec lworkpos = 1 - malpha * workposres / workposwage;
	if (!(iworkneg.is_empty())) {
		uworkneg = u_work_neg(worknegres,cutoff,xi1,alpha,imgamma);
		vec cworkneg = worknegres;
		vec lworkneg = zeros<vec>(worknegres.n_elem);
	}

	// no work cases
	vec unoworkpos = u_no_pos(noworkposres,xi1,mgamma);
	vec cnoworkpos = alpha * noworkposres;
	if (!(inoworkneg.is_empty())){
		unoworkneg = u_no_neg(noworknegres,cutoff,xi1,alpha,imgamma);
		cnoworkneg = noworknegres;
	}

	// reassemble vectors from pos/neg
	uwork.elem(iworkpos) = uworkpos;
	cwork.elem(iworkpos) = cworkpos;
	lwork.elem(iworkpos) = lworkpos;

	if (!(iworkneg.is_empty())){
	//	cout << "hello2!!!!!!" << endl;
	//	cout << "uworkneg" << uworkneg << endl;
	//	cout << "uwork" << uwork << endl;
	//	cout << "iworkneg" << iworkneg << endl;

		uwork.elem(iworkneg) = uworkneg;
	//	cout << "uwork" << uwork << endl;
	//	cout << "cwork" << cwork << endl;
	//	cout << "cworkneg" << worknegres << endl;


		cwork.elem(iworkneg) = worknegres;
		lwork.elem(iworkneg) = lworkneg;
	}
//	cout << unoworkpos << endl;
	unowork.elem(inoworkpos) = unoworkpos;
	cnowork.elem(inoworkpos) = cnoworkpos;
	if (!inoworkneg.is_empty()){
		unowork.elem(inoworkneg) = unoworkneg;
		cnowork.elem(inoworkneg) = noworknegres;
	}

//	cout << "you are here." << endl;

	//reassemble from work/nowork
	util.elem(iwork) = uwork;
	util.elem(inowork) = unowork;
	cons.elem(iwork) = cwork;
	cons.elem(inowork) = cnowork;
	labour.elem(iwork) = lwork;

	// add additive premium if houseing is of right size
	mat phimat = repmat(phivec,1,m);
	mat hfac = exp( theta * phimat);
	util = util % hfac;
	util = util + mu * phimat;
Rcpp::List rlist = Rcpp::List::create( _["utility"] = util, _["consumption"] = cons , _["labour"] = labour);
return rlist;
'

ufun_labour2 <- cxxfunction(signature(Res_="matrix",w_="numeric",s_="numeric",par_="list"),body=cpp.code,plugin="RcppArmadillo",includes="#include </Users/florianoswald/Dropbox/code-backup/eclipse-cpp/arma-tests/arma-tests/src/arma-head.h>")
pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(1:4,5:-1)
w <- runif(4,min=1,max=4)
s <- c(0,0,1,2)
cc <- ufun_labour2(Res,w,s,pars)
rr <- ufun_labouR_h(Res,w,s,pars) 	# same function in R


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




# testing with a header file

cpp.code <- '
using namespace std;
using namespace arma;

// map R objects

mat Res = Rcpp::as<arma::mat>(Res_);
vec wage = Rcpp::as<arma::vec>(w_);
uvec hsize = Rcpp::as<arma::uvec>(s_);
Rcpp::List par( par_ ) ;

// call C++ function
Rcpp::List result = ufun_labour2(Res,par,hsize,wage);

// return
return result;
'
ufun_labour2 <- cxxfunction(signature(Res_="matrix",w_="numeric",s_="numeric",par_="list"),body=cpp.code,plugin="RcppArmadillo",includes="#include </Users/florianoswald/git/ArmaUtils/inst/inc_util.h>")
pars <- list(alpha=0.3,gamma=1.4,cutoff=0.1,mu=0.9,theta=0.3,phival=0.8)
Res  <- outer(1:4,5:-1)
w <- runif(4,min=1,max=4)
s <- c(0,0,1,2)
cc <- ufun_labour2(Res,w,s,pars)
rr <- ufun_labouR_h(Res,w,s,pars) 	# same function in R




