#ifndef SCALEESTIMATION_H
#define SCALEESTIMATION_H

#include <vector>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <iostream>
#include <CGAL/number_utils.h>

namespace ScaleEstimation {
	
template < class ScalarType >
ScalarType opSqrt ( ScalarType s ) { return sqrt(s) ; }



// Computes the median of a vector
// WARNING: When used for scale estimation, the residuals must NOT be squared!
template < class ScalarType >
double median( std::vector< ScalarType > data ) {

	if ( data.size() == 0 ) {
		return 0 ;
	}

	// std::cout << "[ScaleEstimation] Median computation started" << std::endl ;
	double med ;
	size_t size = data.size() ;
		
	std::sort( data.begin(), data.end() ) ;
	// std::cout << "[ScaleEstimation] Input vector sorted " << std::endl ;

	if ( size % 2 == 0 ) {		
		med = ( data[size / 2 - 1] + data[size / 2]) / 2.0 ;
	}
	else {		
		med = data[size / 2] ;
	}
	
	return med ;
}
             


template < class ScalarType >
double medianSE( const std::vector< ScalarType > &sqResiduals, int modelDimension ) {

	size_t size = sqResiduals.size() ;
	// Get median
	double med = median< ScalarType >( sqResiduals ) ;

	// Return Scale Estimation
	return 1.4826 * ( 1.0 + ( 5.0 / ( size - modelDimension - 1 ) ) ) * sqrt( med ) ;
	
}



template < class ScalarType >
double madSE( const std::vector< ScalarType > &sqResiduals, int modelDimension ) {

	size_t size = sqResiduals.size() ;
	std::vector< ScalarType > residuals ; // Note that madSE uses residuals, not squared...
	residuals.resize( size ) ;
	std::transform( sqResiduals.begin(), sqResiduals.end(), residuals.begin(), opSqrt<ScalarType> ) ;

	// Get median
	double med = median< ScalarType >( residuals ) ;

	// Compute differences between median and residuals
	typename std::vector< ScalarType >::const_iterator it ;
	std::vector< ScalarType > absDifs ;
	for( it = residuals.begin(); it != residuals.end(); ++it ) {
		absDifs.push_back( abs( *it - med ) ) ;
	}

	double medDif = median< ScalarType >( absDifs ) ;

	// Return Scale Estimation
	return 1.4826*medDif ;	

}



template < class ScalarType >
double msSE( std::vector< ScalarType > sqResiduals, int modelDimension, double quantile = 0.2, double T = 2.5 ) {

	size_t size = sqResiduals.size() ;
	//std::cout << "residuals size = " << size << std::endl ;
		
	// Sort residuals
	std::sort( sqResiduals.begin(), sqResiduals.end() ) ;

	// compute the last inlier according to quantile
	int q = (int)floor( size * quantile ) - 1 ; // 0-based indexing!
	// std::cout << "q = " << q << std::endl ;
	if ( q < 0 ) {
		return sqResiduals[ sqResiduals.size()-1 ] ;
	}

	// Starting from the last inlier, increment the index iteratively till a big jump is found
	bool bigJump = false ;	
	double sumSqRes = std::accumulate( sqResiduals.begin(), sqResiduals.begin() + q, 0.0 ) ;
	double den = ( q - modelDimension + 1 ) ;
	if ( den <= 0 ) { // Workaround for when q < modelDimension, is there a better way?
		//std::cout << "q < modelDimension!" << std::endl ;
		return sqResiduals[q] ; // Worst case (?) estimate
	}
	double sqSigmaPrev = sumSqRes / den ; // First sigma computation (refer to the original paper for the definition of sigma)
	//std::cout << "sumSqRes = " << sumSqRes << std::endl ;
	//std::cout << "sqSigmaPrev = " << sqSigmaPrev << std::endl ;
	//std::cout << "den = " << den << std::endl ;
	//std::cout << "q = " << den << std::endl ;
	double sqSigmaCur = 0.0 ;
	q++ ;
	while ( !bigJump && q < size ) {
		
		sumSqRes += sqResiduals[q] ;
		sqSigmaCur = sumSqRes / ( q - modelDimension + 1 ) ;

		//std::cout << "sqResiduals[" << q << "] = " << sqResiduals[q] << std::endl ;
		//std::cout << "sumSqRes = " << sumSqRes << std::endl ;
		//std::cout << "sqSigmaCur = " << sqSigmaCur << std::endl ;

		// Is there a big jump?
		double den2 = ( q - 1 ) - modelDimension + 1.0 ;
		if ( den2 <= 0 ) // Workaround for when q < modelDimension, is there a better way?
			return sqResiduals[q] ; // Worst case (?) estimate
		bigJump = ( sqSigmaCur / sqSigmaPrev ) > ( 1.0 + ( ( T*T - 1.0 ) / ( den2 ) ) ) ;
		
		if ( !bigJump ) {
			// Update for next iteration
			sqSigmaPrev = sqSigmaCur ;
			q++ ;
		}
		// else
		//	std::cout << "bigJump!" << std::endl ;
	}
	
	double scale = sqrt( sqSigmaPrev ) ;
	
	// std::cout << "scale = " << sqrt( sqSigmaPrev ) << std::endl ;

	// Return Scale Estimation
	return scale ;

}

} // End namespace ScaleEstimation

#endif // #ifndef SCALEESTIMATION_H