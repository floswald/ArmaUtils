#ifndef _ArmaUtils_ARMAFUNS_H
#define _ArmaUtils_ARMAFUNS_H

#include <RcppArmadillo.h>
#include "handler.h"


RcppExport SEXP armamax( SEXP X_) ;
RcppExport SEXP armamax2( SEXP X_, SEXP B_ ) ;
RcppExport SEXP armakron( SEXP y_, SEXP R_ ) ;
RcppExport SEXP boom( SEXP A_ ) ;

#endif
