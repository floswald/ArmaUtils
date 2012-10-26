#ifndef _RcppUtils_UFUNS_H
#define _RcppUtils_UFUNS_H

#include <RcppArmadillo.h>

RcppExport SEXP utilfun(  SEXP Res_, SEXP s_, SEXP age_, SEXP par_, SEXP xi_, SEXP own_ ) ; 
RcppExport SEXP ufun_Attanasio( SEXP Res_, SEXP s_, SEXP par_) ; 

#endif
