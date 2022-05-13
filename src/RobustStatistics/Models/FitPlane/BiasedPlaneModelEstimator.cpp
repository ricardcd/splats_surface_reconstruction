#include "RobustStatistics/Models/FitPlane/BiasedPlaneModelEstimator.h"

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dimension.h>

BiasedPlaneModelEstimator::BiasedPlaneModelEstimator( double distThres, Point_3 forcedInlier ) : ModelEstimator< Point_3, double >( 3 ), m_sqDistanceThreshold( distThres*distThres ), m_forcedInlier( forcedInlier ) {}  

void BiasedPlaneModelEstimator::fit(	const std::vector< std::vector<Point_3>::const_pointer > &data,
								std::vector< double > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	Plane_3 plane( *data[0], *data[1], *data[2] ) ;

	parameters.push_back( plane.a() ) ;
	parameters.push_back( plane.b() ) ;
	parameters.push_back( plane.c() ) ;
	parameters.push_back( plane.d() ) ;

}



void BiasedPlaneModelEstimator::leastSquaresFit(	const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
											std::vector< double > &parameters ) const {

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

	// Do the least squares fitting (CGAL)
	Plane_3 plane ;
	Point_3 centroid ;
	double res = linear_least_squares_fitting_3( dataRef.begin(), dataRef.end(), plane, centroid, CGAL::Dimension_tag<0>() ) ;

	// Fill parameters
	parameters.push_back( plane.a() ) ;
	parameters.push_back( plane.b() ) ;
	parameters.push_back( plane.c() ) ;
	parameters.push_back( plane.d() ) ;

}


	
bool BiasedPlaneModelEstimator::compatible(	const std::vector< double > &parameters, 
										const Point_3 &data ) const {

	Plane_3 plane( parameters[0], parameters[1], parameters[2], parameters[3] ) ;
	double sqDistance = squared_distance( data, plane ) ;
	
	return sqDistance < ( m_sqDistanceThreshold ) ;

}



bool BiasedPlaneModelEstimator::degenerate( const std::vector< typename std::vector< Point_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}

	Plane_3 plane( *data[0], *data[1], *data[2] ) ;
	
	//std::cout << "[RANSAC] First Point = " << data[0]->x() << data[0]->y() << data[0]->z() << std::endl ;
	//std::cout << "[RANSAC] Plane = " << plane.a() << " " << plane.b() << " " << plane.c() << " " << plane.d() << " " << std::endl ;
	/*if ( plane.is_degenerate() ) 
		std::cout << "[RANSAC] Degenerate." << std::endl ;
	else
		std::cout << "[RANSAC] Good sample set." << std::endl ;*/

	// This point must be an inlier of the model for this plane to be taken into account
	bool bias = squared_distance( m_forcedInlier, plane ) < m_sqDistanceThreshold ;

	return plane.is_degenerate() && bias ;

}