#include "LocalSurfaceFit/WLBQModelEstimator.h"
#include <CGAL/linear_least_squares_fitting_3.h>



WLBQModelEstimator::WLBQModelEstimator( const Point_3 &refPoint, 
										const FT &distThres, 
										const FT &gaussianH ) : ModelEstimator< Point_3, FT >( 6 ), 
																m_sqDistanceThreshold( distThres*distThres ), 
																m_refPoint( refPoint ),
																m_gaussianH( gaussianH )
{}



void WLBQModelEstimator::fit(	const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
																std::vector< FT > &parameters ) const {
	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	// Change Point_3* data vector into a Point_3 vector
	std::vector< Point_3 > dataRef ;
	std::vector< typename std::vector< Point_3 >::const_pointer >::const_iterator itDataPtr ;
	for ( itDataPtr = data.begin(); itDataPtr != data.end(); ++itDataPtr ) {
		dataRef.push_back( **itDataPtr ) ;
	}

	// Compute PCA basis from the points
	Plane_3 plane ;
	Point_3 centroid ;
	double res = linear_least_squares_fitting_3( dataRef.begin(), dataRef.end(), plane, centroid, CGAL::Dimension_tag<0>() ) ;

	Vector_3 planeNormal = plane.orthogonal_vector() ;

	// Compute weights according to reference point
	// TODO: Weights should also modify the PCA?
	std::vector< FT > weights ;
	for( std::vector< Point_3 >::iterator itP = dataRef.begin(); itP != dataRef.end(); ++itP ) {
		FT w = LocalBivariateQuadric::weightGaussian( CGAL::sqrt( CGAL::squared_distance( *itP, m_refPoint ) ),
													  m_gaussianH ) ;
		weights.push_back( w ) ;
	}

	// Compute the LBQ on the computed plane
	LocalBivariateQuadric lbq( centroid, planeNormal ) ;
	lbq.weightedFit( dataRef, weights ) ;

	parameters = LBQ2parameters( lbq ) ;

}



void WLBQModelEstimator::leastSquaresFit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
								  std::vector< FT > &parameters ) const {

	// Fitting is already a least squares process...
	fit( data, parameters ) ;	
	
}



bool WLBQModelEstimator::compatible(	const std::vector< double > &parameters, 
									const Point_3 &data ) const {

	// Rebuild the LBQ from the parameters
	LocalBivariateQuadric lbq = parameters2LBQ( parameters ) ;

	// Projection!!
	FT sqDistance = CGAL::squared_distance( data, lbq.projection( data ) ) ;

	return sqDistance < m_sqDistanceThreshold ;

}



std::vector< double > WLBQModelEstimator::sqResiduals( const std::vector< FT > &parameters, 
													  const std::vector< Point_3 > &data ) const {
	
	// Rebuild the LBQ from the parameters
	LocalBivariateQuadric lbq = parameters2LBQ( parameters ) ;

	std::vector< double > sqResiduals ;	
	std::vector< Point_3 >::const_iterator it ;
	for( it = data.begin(); it != data.end(); ++it ) {
		// Orthogonal distance! There is no closed form to get the euclidean distance from a point to a quadric...
		double sqDistance = CGAL::squared_distance( *it, lbq.projection( *it ) ) ;
		
		// --- Debug (Start) ---
		// std::cout << "SqDistance = " << sqDistance << std::endl ;
		// --- Debug  (End)  ---

		sqResiduals.push_back( abs( sqDistance ) ) ;
	}

	return sqResiduals ;

}



bool WLBQModelEstimator::degenerate( const std::vector< std::vector< Point_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}

	return false ;

}



// Rebuild the LBQ from parameters vector
WLBQModelEstimator::LocalBivariateQuadric	WLBQModelEstimator::parameters2LBQ( const std::vector< FT > &parameters ) {
		
	Point_3 origin	( parameters[ 0 ], parameters[ 1 ], parameters[ 2 ] ) ;
	Vector_3 normal	( parameters[ 3 ], parameters[ 4 ], parameters[ 5 ] ) ;
	Vector_3 xAxis	( parameters[ 6 ], parameters[ 7 ], parameters[ 8 ] ) ;
	Vector_3 yAxis	( parameters[ 9 ], parameters[ 10 ], parameters[ 11 ] ) ;
	std::vector<FT> coeffs ;
	coeffs.resize(6) ;
	coeffs[0] = parameters[12] ;
	coeffs[1] = parameters[13] ;
	coeffs[2] = parameters[14] ;
	coeffs[3] = parameters[15] ;
	coeffs[4] = parameters[16] ;
	coeffs[5] = parameters[17] ;
	LocalBivariateQuadric lbq( origin, normal, xAxis, yAxis, coeffs ) ;

	return lbq ;

}



// Build the parameters vector from the LBQ
std::vector< FT > WLBQModelEstimator::LBQ2parameters( const LocalBivariateQuadric& lbq ) {

	std::vector< FT > parameters ;

	// Origin
	parameters.push_back( lbq.origin().x() ) ;
	parameters.push_back( lbq.origin().y() ) ;
	parameters.push_back( lbq.origin().z() ) ;
	// Normal
	parameters.push_back( lbq.axisZ().x() ) ;
	parameters.push_back( lbq.axisZ().y() ) ;
	parameters.push_back( lbq.axisZ().z() ) ;
	// X Axis
	parameters.push_back( lbq.axisX().x() ) ;
	parameters.push_back( lbq.axisX().y() ) ;
	parameters.push_back( lbq.axisX().z() ) ;
	// Y Axis
	parameters.push_back( lbq.axisY().x() ) ;
	parameters.push_back( lbq.axisY().y() ) ;
	parameters.push_back( lbq.axisY().z() ) ;
	// Coefficients
	std::vector< FT > coeffs = lbq.coefficients() ;
	parameters.push_back( coeffs[0] ) ;
	parameters.push_back( coeffs[1] ) ;
	parameters.push_back( coeffs[2] ) ;
	parameters.push_back( coeffs[3] ) ;
	parameters.push_back( coeffs[4] ) ;
	parameters.push_back( coeffs[5] ) ;

	return parameters ;

}
