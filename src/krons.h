
#ifndef _ArmaUtils_KRONS_H
#define _ArmaUtils_KRONS_H

#include <RcppEigen.h>
#include "handler.h"

RcppExport SEXP kron4( SEXP aa0, SEXP aa1, SEXP aa2, SEXP aa3, SEXP yy ) ; 
RcppExport SEXP kron3( SEXP aa0, SEXP aa1, SEXP aa2, SEXP yy ) ; 
RcppExport SEXP kron2( SEXP aa0, SEXP aa1, SEXP yy ) ; 

#endif
