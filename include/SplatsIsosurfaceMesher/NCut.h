/* 
	Normalized Cuts algorithm
	--> Slight modification of the implementation by Evan Herbst:
		web: http://homes.cs.washington.edu/~eherbst/useful-code/ncuts.php
*/

#ifndef NCUT_H
#define NCUT_H

/* INCLUDES */
// Std
#include <vector>
#include <algorithm>
// Eigen
#include <Eigen/Core>
#include <Eigen/Sparse>
// ARPACK++
#include <arssym.h>
#include <arrseig.h>



template <typename SparseMatrixT>
class arpackWeightMatrixMultiplierFunctor
{
	typedef Eigen::Matrix< double, 
						   Eigen::Dynamic, 
						   1 >						EVector ;
	typedef Eigen::Matrix< double,
						   ::Eigen::Dynamic,
						   ::Eigen::Dynamic>		EMatrix ;
	typedef Eigen::SparseMatrix< double >			SparseMatrix ;

	public:

		arpackWeightMatrixMultiplierFunctor(const SparseMatrixT& w) : m(w)
		{}
		virtual ~arpackWeightMatrixMultiplierFunctor() {}

		/*
			* b <- A * x
			*/
		void operator () (double* const x, double* const b)
		{
			//TODO make use of sparse
			EVector xv(m.cols());
			for(unsigned int i = 0; i < m.cols(); i++) xv[i] = x[i];
			EVector bv = m * xv;
			for(unsigned int i = 0; i < m.rows(); i++) b[i] = bv[i];
		}

	private:

		SparseMatrixT m;
};



template <class FT>
double sign(FT x)
{
	return (x < 0) ? -1 : (x > 0) ? 1 : 0;
}



template < class FT >
void ncut(	const Eigen::SparseMatrix<FT>& W,
			const unsigned int nbEigenValues,
			Eigen::Matrix< FT,
						   ::Eigen::Dynamic,
						   ::Eigen::Dynamic>& eigenVectors,
			Eigen::Matrix< FT, 
						   Eigen::Dynamic, 
						   1 >& eigenValues )
{
	// Defines
	typedef Eigen::Matrix< FT, 
						   Eigen::Dynamic, 
						   1 >						EVector ;
	typedef Eigen::Matrix< FT,
						   ::Eigen::Dynamic,
						   ::Eigen::Dynamic>		EMatrix ;
	typedef Eigen::SparseMatrix< FT >				SparseMatrix ;


	int numData = W.rows() ;
	SparseMatrix M( numData, numData ) ;
	M.reserve( numData ) ;
	for( unsigned int j = 0; j < W.rows(); j++ ) {
		FT val = 0 ;
		for( typename SparseMatrix::InnerIterator it( W, j ); it; ++it ) {
			val += fabs( it.value() ) ;
		}
		// std::cout << val << std::endl ;
		// M.insert( j, j ) = val ;
		// WARNING: We perform the inv( sqrt ) here!
		M.insert( j, j ) = 1 / CGAL::sqrt( val ) ;
	}
		
	SparseMatrix C = M*W*M ;

	// std::cout << C << std::endl ;

	EMatrix usEigenVectors( numData, nbEigenValues ) ; //one per column
	EVector usEigenValues( nbEigenValues ) ;

	/* ARPACK Solve */
	{
		const unsigned int maxIters = 50 ;
		const char* const whichEigvals = "LA"; //LA = largest algebraic
	//	const unsigned int numArnoldiVectors = MIN(2 * numEigenvals + 1, numData - 1); //value suggested by arpack
		const unsigned int numArnoldiVectors = std::min<int>(std::max<int>(35u, 2 * nbEigenValues), numData ); //saupp.h recommends to increase ncv if you continually get the "no shifts could be applied" error msg (but that does mean increasing space & time)
		const double relativeTol = 1e-4 ;

		// arpackWeightMatrixMultiplierFunctor< decltype(C) > wmmf( C ) ;
		arpackWeightMatrixMultiplierFunctor< SparseMatrix > wmmf( C ) ;
		//ARSymStdEig< FT, decltype( wmmf ) > solver( C.rows(), 
		ARSymStdEig< FT, arpackWeightMatrixMultiplierFunctor< SparseMatrix > > solver( C.rows(), 
													nbEigenValues, 
													&wmmf, 
													&arpackWeightMatrixMultiplierFunctor<SparseMatrix>::operator (), 
													const_cast<char*>(whichEigvals), 
													numArnoldiVectors, 
													relativeTol, 
													maxIters ) ;
		solver.FindEigenvectors() ;

		assert(solver.EigenvaluesFound());
		assert(solver.ConvergedEigenvalues() >= nbEigenValues);
		// cout << "eigenvalues (want " << nbEigenValues << "): <";
		for(unsigned int i = 0; i < nbEigenValues; i++)
		{
			// cout << " " << solver.Eigenvalue(i);
			usEigenValues[i] = solver.Eigenvalue(i);
		}
		// cout << " >" << endl;


		assert(solver.EigenvectorsFound());

		// OLD CODE:
//		// cout << " Eigenvectors found" << endl;
//		const std::vector<double>* concatenatedEigvecs = solver.StlEigenvectors(); //all concatenated into one vector
//		for(unsigned int i = 0; i < nbEigenValues; i++)
//		{
//			std::copy(concatenatedEigvecs->begin() + i * numData, concatenatedEigvecs->begin() + (i + 1) * numData, usEigenVectors.col(i).data());
//			// cout << "eigenvec " << usEigenVectors.col(i) << endl;
//		}

        // NEW CODE (following: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=230)
        for(unsigned int i = 0; i < nbEigenValues; i++) {
            double* evp = solver.RawEigenvector(i);
            for (unsigned int j = 0; j < numData; j++) {
                usEigenVectors(j, i) = evp[j];
            }
        }
	}

	/*
	* sort eigvals and rearrange eigvecs
	*/
	std::vector< std::pair< FT, unsigned int > > sds( usEigenValues.rows() ) ;
	for(unsigned int i = 0; i < usEigenValues.rows(); i++)
	{
		sds[i].first = usEigenValues[i];
		sds[i].second = i;
	}
	std::sort( sds.begin(), sds.end() ) ;
		
	eigenValues = EVector( usEigenValues.rows() ) ;
	for(unsigned int i = 0; i < usEigenValues.rows(); i++) 
		eigenValues[i] = sds[i].first;
	EMatrix vbar(usEigenVectors.rows(), nbEigenValues);
	for(unsigned int i = 0; i < usEigenValues.rows(); i++) 
		vbar.col(i) = usEigenVectors.col(sds[i].second);

	eigenVectors = M * vbar;

	for(unsigned int i = 0; i < eigenVectors.cols(); i++)
	{
		eigenVectors.col(i) *= sqrt( static_cast<double>( numData ) ) / eigenVectors.col(i).norm() ;
		if( eigenVectors(0, i) != 0 )
			eigenVectors.col(i) = eigenVectors.col(i) * -sign(eigenVectors(0, i)) ;
	}
		
}





/* Legacy... */
//std::tuple< EMatrix, EVector > 
//ncutFull( const EMatrix& W, 
//		  const unsigned int nbEigenValues )
//{
//	int numData = W.rows() ;
//	EMatrix M( numData, numData ) ;
//	for ( int i = 0; i < numData; i++ ) 
//	{
//		M( i, i ) = W.row(i).sum() ;
//	}

//	EMatrix B = M.array().sqrt() ; 
//	B = B.inverse() ;

//	EMatrix C = B*W*B ;

//	Eigen::EigenSolver< EMatrix > es( C ) ;
//	std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl ;
//	std::cout << "The second eigenvector is:" << std::endl << es.eigenvectors().col(1) << std::endl << endl ;

//	std::make_tuple( es.eigenvectors(), es.eigenvalues() ) ;
//}


#endif // NCUT_H