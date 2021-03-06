


#include "krons.h"
#include <iostream>

using namespace Rcpp ;
using namespace Eigen ;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;

Eigen::VectorXd KronProdSPMat4(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3,
		Eigen::VectorXd y) {

    signal(SIGSEGV, handler);   // install our handler
	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows() * a3.rows()  );

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

Eigen::VectorXd KronProdSPMat3(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::VectorXd y) {

    signal(SIGSEGV, handler);   // install our handler
	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows()  );

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
						for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it2(a2,row_idx2); it2; ++it2){
								retvec( row_offset2 + row_idx2 ) += factor2 * it2.value() * y( col_offset2 + it2.index() );
						}
					}
				}
			}
		}
	}
	return retvec;
}


Eigen::VectorXd KronProdSPMat2(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::VectorXd y) {

    signal(SIGSEGV, handler);   // install our handler
	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows()  );

	//loop rows a0
	for (int row_idx0=0; row_idx0<a0.outerSize(); ++row_idx0) {
		int row_offset1 = row_idx0;
		row_offset1    *= a1.rows();

		// loop rows a1
		for (int row_idx1=0; row_idx1<a1.outerSize(); ++row_idx1) {

			// loop cols a0 (non-zero elements only)
			for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it0(a0,row_idx0); it0; ++it0) {
				int col_offset1 = it0.index();
				col_offset1    *= a1.innerSize();
				double factor1 = it0.value();

				for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it1(a1,row_idx1); it1; ++it1){
						retvec( row_offset1 + row_idx1 ) += factor1 * it1.value() * y( col_offset1 + it1.index() );
				}
			}
		}
	}
	return retvec;
}

SEXP kron4( SEXP aa0, SEXP aa1, SEXP aa2, SEXP aa3, SEXP yy )
{
    signal(SIGSEGV, handler);   // install our handler

	BEGIN_RCPP

	const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
	const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
	const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
	const MappedSparseMatrix<double> a3(as<MappedSparseMatrix<double> >(aa3));
	const MapVecd y(as<MapVecd>(yy));

	if ( y.size() != a0.cols() * a1.cols() * a2.cols()  * a3.cols()  ) {
		throw Rcpp::exception( "ArmaUtils::kronProdSPMat4 error: length(y) != a0.cols()*a1.cols()*a2.cols()*a3.cols()" );
		return R_NilValue;
	}
	VectorXd result( KronProdSPMat4( a0,a1,a2,a3,y ) );
	return wrap( result );

	END_RCPP

}

SEXP kron3( SEXP aa0, SEXP aa1, SEXP aa2, SEXP yy )
{

    signal(SIGSEGV, handler);   // install our handler
	BEGIN_RCPP

	const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
	const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
	const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
	const MapVecd y(as<MapVecd>(yy));

	if ( y.size() != a0.cols() * a1.cols() * a2.cols()) {
		throw Rcpp::exception( "ArmaUtils::kronProdSPMat3 error: length(y) != a0.cols()*a1.cols()*a2.cols()" );
		return R_NilValue;
	}

	VectorXd result( KronProdSPMat3( a0,a1,a2,y ) );
	return wrap( result );

	END_RCPP
}


SEXP kron2( SEXP aa0, SEXP aa1, SEXP yy )
{
    signal(SIGSEGV, handler);   // install our handler
 
	BEGIN_RCPP

	const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
	const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
	const MapVecd y(as<MapVecd>(yy));

	if ( y.size() != a0.cols() * a1.cols() ) {
		throw Rcpp::exception( "ArmaUtils::kronProdSPMat2 error: length(y) != a0.cols()*a1.cols()" );
		return R_NilValue;
	}
	VectorXd result( KronProdSPMat2( a0,a1,y ) );
	return wrap( result );

	END_RCPP
}

