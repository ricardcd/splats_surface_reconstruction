#include "RobustStatistics/Models/FitCentroid/CentroidModelEstimator.h"
#include <math.h>
#include <CGAL/centroid.h>

CentroidModelEstimator::CentroidModelEstimator( double distThres ) : ModelEstimator< Point_3, double >( 1 ), m_distanceThreshold( distThres ) {} 



void CentroidModelEstimator::fit( const std::vector< std::vector<Point_3>::const_pointer > &data,
							  std::vector< double > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	parameters.push_back( data[0]->x() ) ;
	parameters.push_back( data[0]->y() ) ;
	parameters.push_back( data[0]->z() ) ;
	
}



void CentroidModelEstimator::leastSquaresFit( const std::vector< std::vector<Point_3>::const_pointer > &data,
										  std::vector< double > &parameters ) const {

	
	parameters.clear();
	if(data.size()<this->minSamples())
		return;

	// Change Point_3* data vector into a Point_3 vector
	std::vector< Point_3 > dataRef ;
	std::vector< std::vector< Point_3 >::const_pointer >::const_iterator itDataPtr ;
	for ( itDataPtr = data.begin(); itDataPtr != data.end(); ++itDataPtr ) {
		dataRef.push_back( **itDataPtr ) ;
	}

	Point_3 c = CGAL::centroid( dataRef.begin(), dataRef.end() ) ;

	parameters.push_back( c.x() ) ;
	parameters.push_back( c.y() ) ;
	parameters.push_back( c.z() ) ;

}



bool CentroidModelEstimator::compatible( const std::vector< double > &parameters, 
									 const Point_3 &data ) const {

	double sqDist = CGAL::squared_distance( Point_3( parameters[0], parameters[1], parameters[2] ), data ) ; 
	return sqDist < m_distanceThreshold*m_distanceThreshold ;

}



std::vector< double > CentroidModelEstimator::sqResiduals(  const std::vector< double > &parameters, 
													  const std::vector< Point_3 > &data ) const {

	// NOT USED, NOT IMPLEMENTED

	return std::vector< double >() ;
}



bool CentroidModelEstimator::degenerate( const std::vector< std::vector<Point_3>::const_pointer > &data ) const {
	
	return false ;

}