#include "RobustStatistics/Models/FitLocalSmoothSurface/LocalSmoothSurfaceModelEstimator.h"

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dimension.h>

#include <iostream>

LocalSmoothSurfaceModelEstimator::LocalSmoothSurfaceModelEstimator( unsigned int minSamples, unsigned int degree, Point_3 refPoint, FT distThres ) : ModelEstimator< Point_3, FT >( minSamples ), 
																																   m_sqDistanceThreshold( distThres*distThres ) {

	m_degree = degree ;
	if ( m_degree <= 0 )
		m_degree = 1 ;
	if ( m_degree > 4 )
		m_degree = 4 ;

	m_refPoint = refPoint ;

}  



void LocalSmoothSurfaceModelEstimator::fit(	const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
											std::vector< FT > &parameters ) const {

	// There is no closed form solution for a minimum number of points, using the least-squares approximation	
	// leastSquaresFit( data, parameters ) ;
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

	// Compute the monge form of the local surface
	Monge_via_jet_fitting monge_fit ;
	Monge_form monge_form = monge_fit( dataRef.begin(), dataRef.end(), m_degree, m_degree ) ;
	
	/* Fill parameters */
	// Origin
	Point_3 origin = monge_form.origin() ;
	parameters.push_back( origin.x() ) ;
	parameters.push_back( origin.y() ) ;
	parameters.push_back( origin.z() ) ;
	// Maximal direction
	Vector_3 d1 = monge_form.maximal_principal_direction() ;
	parameters.push_back( d1.x() ) ;
	parameters.push_back( d1.y() ) ;
	parameters.push_back( d1.z() ) ;
	// Minimal direction
	Vector_3 d2 = monge_form.minimal_principal_direction() ;
	parameters.push_back( d2.x() ) ;
	parameters.push_back( d2.y() ) ;
	parameters.push_back( d2.z() ) ;
	// Normal direction
	Vector_3 n = monge_form.normal_direction() ;
	parameters.push_back( n.x() ) ;
	parameters.push_back( n.y() ) ;
	parameters.push_back( n.z() ) ;
	// Coefficients
	if ( m_degree > 1 ) {
		parameters.push_back( monge_form.coefficients()[0] ) ;
		parameters.push_back( monge_form.coefficients()[1] ) ;
	}
	if ( m_degree > 2 ) {
		parameters.push_back( monge_form.coefficients()[2] ) ;
		parameters.push_back( monge_form.coefficients()[3] ) ;
		parameters.push_back( monge_form.coefficients()[4] ) ;
		parameters.push_back( monge_form.coefficients()[5] ) ;
	}
	if ( m_degree > 3 ) {
		parameters.push_back( monge_form.coefficients()[6] ) ;
		parameters.push_back( monge_form.coefficients()[7] ) ;
		parameters.push_back( monge_form.coefficients()[8] ) ;
		parameters.push_back( monge_form.coefficients()[9] ) ;
		parameters.push_back( monge_form.coefficients()[10] ) ;
	}

}



void LocalSmoothSurfaceModelEstimator::leastSquaresFit(	const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
														std::vector< FT > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	// Change Point_3* data vector into a Point_3 vector
	std::vector< Point_3 > dataRef ;
	std::vector< std::vector< Point_3 >::const_pointer >::const_iterator itDataPtr ;
	for ( itDataPtr = data.begin(); itDataPtr != data.end(); ++itDataPtr ) {
		dataRef.push_back( **itDataPtr ) ;
	}

	// Put the reference point the first (if exists)
	std::vector< Point_3 >::iterator refPIt = std::find( dataRef.begin(), dataRef.end(), m_refPoint ) ;
	if ( refPIt != dataRef.begin() && refPIt != dataRef.end() ) {		
		std::cout << "Reference point swapt!!!!" << std::endl ;
		std::swap( *refPIt, *( dataRef.begin() ) ) ;
	}

	// Compute the monge form of the local surface
	Monge_via_jet_fitting monge_fit ;
	Monge_form monge_form = monge_fit( dataRef.begin(), dataRef.end(), m_degree, m_degree ) ;
	
	/* Fill parameters */
	// Origin
	Point_3 origin = monge_form.origin() ;
	parameters.push_back( origin.x() ) ;
	parameters.push_back( origin.y() ) ;
	parameters.push_back( origin.z() ) ;
	// Maximal direction
	Vector_3 d1 = monge_form.maximal_principal_direction() ;
	parameters.push_back( d1.x() ) ;
	parameters.push_back( d1.y() ) ;
	parameters.push_back( d1.z() ) ;
	// Minimal direction
	Vector_3 d2 = monge_form.minimal_principal_direction() ;
	parameters.push_back( d2.x() ) ;
	parameters.push_back( d2.y() ) ;
	parameters.push_back( d2.z() ) ;
	// Normal direction
	Vector_3 n = monge_form.normal_direction() ;
	parameters.push_back( n.x() ) ;
	parameters.push_back( n.y() ) ;
	parameters.push_back( n.z() ) ;
	// Coefficients
	if ( m_degree > 1 ) {
		parameters.push_back( monge_form.coefficients()[0] ) ;
		parameters.push_back( monge_form.coefficients()[1] ) ;
	}
	if ( m_degree > 2 ) {
		parameters.push_back( monge_form.coefficients()[2] ) ;
		parameters.push_back( monge_form.coefficients()[3] ) ;
		parameters.push_back( monge_form.coefficients()[4] ) ;
		parameters.push_back( monge_form.coefficients()[5] ) ;
	}
	if ( m_degree > 3 ) {
		parameters.push_back( monge_form.coefficients()[6] ) ;
		parameters.push_back( monge_form.coefficients()[7] ) ;
		parameters.push_back( monge_form.coefficients()[8] ) ;
		parameters.push_back( monge_form.coefficients()[9] ) ;
		parameters.push_back( monge_form.coefficients()[10] ) ;
	}
	
}


	
bool LocalSmoothSurfaceModelEstimator::compatible(	const std::vector< FT > &parameters, 
													const Point_3 &data ) const {

	// Rebuild the monge
	Point_3 origin( parameters[0], parameters[1], parameters[2] ) ;
	Vector_3 d1( parameters[3], parameters[4], parameters[5] ) ;
	Vector_3 d2( parameters[6], parameters[7], parameters[8] ) ;
	Vector_3 n( parameters[9], parameters[10], parameters[11] ) ;
	std::vector< FT > coefs ;
	if ( m_degree > 1 ) {
		coefs.push_back( parameters[12] ) ;
		coefs.push_back( parameters[13] ) ;			
	}
	if ( m_degree > 2 ) {
		coefs.push_back( parameters[14] ) ;
		coefs.push_back( parameters[15] ) ;	
		coefs.push_back( parameters[16] ) ;
		coefs.push_back( parameters[17] ) ;		
	}
	if ( m_degree > 3 ) {
		coefs.push_back( parameters[18] ) ;		
		coefs.push_back( parameters[19] ) ;		
		coefs.push_back( parameters[20] ) ;		
		coefs.push_back( parameters[21] ) ;		
		coefs.push_back( parameters[22] ) ;	
	}	

	Monge_form monge = Monge_form( origin, d1, d2, n, coefs ) ;											
											
	// Orthogonal distance! There is no closed form to get the euclidean distance from a point to a quadric...
	FT sqDistance = squared_distance( data, monge.project( data ) ) ;

	// --- Debug (Start) ---
	// std::cout << "Distance  = " << squared_distance( data, monge.project( data ) ) << ",  m_sqDistanceThreshold = " << m_sqDistanceThreshold << std::endl ;
	// --- Debug  (End)  ---

	
	return sqDistance < m_sqDistanceThreshold ;

}

std::vector< double > LocalSmoothSurfaceModelEstimator::sqResiduals(  const std::vector< FT > &parameters, 
																	  const std::vector< Point_3 > &data ) const {

	// Rebuild the monge
	Point_3 origin( parameters[0], parameters[1], parameters[2] ) ;
	Vector_3 d1( parameters[3], parameters[4], parameters[5] ) ;
	Vector_3 d2( parameters[6], parameters[7], parameters[8] ) ;
	Vector_3 n( parameters[9], parameters[10], parameters[11] ) ;
	std::vector< FT > coefs ;
	if ( m_degree > 1 ) {
		coefs.push_back( parameters[12] ) ;
		coefs.push_back( parameters[13] ) ;			
	}
	if ( m_degree > 2 ) {
		coefs.push_back( parameters[14] ) ;
		coefs.push_back( parameters[15] ) ;	
		coefs.push_back( parameters[16] ) ;
		coefs.push_back( parameters[17] ) ;		
	}
	if ( m_degree > 3 ) {
		coefs.push_back( parameters[18] ) ;		
		coefs.push_back( parameters[19] ) ;		
		coefs.push_back( parameters[20] ) ;		
		coefs.push_back( parameters[21] ) ;		
		coefs.push_back( parameters[22] ) ;	
	}	
	
	Monge_form monge = Monge_form( origin, d1, d2, n, coefs ) ;											
											
	std::vector< Point_3 >::const_iterator it ;
	std::vector< double > sqResiduals ;
	
	for( it = data.begin(); it != data.end(); ++it ) {
		// Orthogonal distance! There is no closed form to get the euclidean distance from a point to a quadric...
		double sqDistance = CGAL::squared_distance( *it, monge.project( *it ) ) ;
		
		// --- Debug (Start) ---
		// std::cout << "SqDistance = " << sqDistance << std::endl ;
		// --- Debug  (End)  ---

		sqResiduals.push_back( abs( sqDistance ) ) ;
	}

	return sqResiduals ;
	
}

bool LocalSmoothSurfaceModelEstimator::degenerate( const std::vector< typename std::vector< Point_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}
	else {
		// TODO: How to tell if a monge is degenerate?
		return false ;
	}

}