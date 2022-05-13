#include "QuadricModelEstimator.h"

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dimension.h>

QuadricModelEstimator::QuadricModelEstimator( unsigned int minSamples, FT distThres ) : ModelEstimator< Point_3, FT >( minSamples ), m_sqDistanceThreshold( distThres*distThres ) {}  

void QuadricModelEstimator::fit(	const std::vector< std::vector< Point_3 >::const_pointer > &data,
								std::vector< FT > &parameters ) const {

	// There is no closed form solution for a minimum number of points, using the least-squares approximation
	leastSquaresFit( data, parameters ) ;

}



void QuadricModelEstimator::leastSquaresFit(	const std::vector< std::vector< Point_3 >::const_pointer > &data,
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

	// Compute the monge form of the local surface
	Monge_via_jet_fitting monge_fit ;
	Monge_form monge_form = monge_fit( dataRef.begin(), dataRef.end(), 2, 2 ) ;
	
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
	parameters.push_back( monge_form.coefficients()[0] ) ;
	parameters.push_back( monge_form.coefficients()[1] ) ;
	
}


	
bool QuadricModelEstimator::compatible(	const std::vector< FT > &parameters, 
										const Point_3 &data ) const {

	// Rebuild the monge
	Point_3 origin( parameters[0], parameters[1], parameters[2] ) ;
	Vector_3 d1( parameters[3], parameters[4], parameters[5] ) ;
	Vector_3 d2( parameters[6], parameters[7], parameters[8] ) ;
	Vector_3 n( parameters[9], parameters[10], parameters[11] ) ;
	std::vector< FT > coefs ;
	coefs.push_back( parameters[12] ) ;
	coefs.push_back( parameters[13] ) ;	

	Monge_form monge = Monge_form( origin, d1, d2, n, coefs ) ;											
											
	// Orthogonal distance! There is no closed form to get the euclidean distance from a point to a quadric...
	double sqDistance = squared_distance( data, monge.project( data ) ) ;
	
	return sqDistance < ( m_sqDistanceThreshold ) ;

}



bool QuadricModelEstimator::degenerate( const std::vector< std::vector< Point_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}
	else {
		// TODO: How to tell if a monge is degenerate?
		return false ;
	}

}