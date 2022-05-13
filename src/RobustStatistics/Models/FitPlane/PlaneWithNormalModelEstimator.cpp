#include "RobustStatistics/Models/FitPlane/PlaneWithNormalModelEstimator.h"

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Dimension.h>

PlaneWithNormalModelEstimator::PlaneWithNormalModelEstimator( double distThres, double angleThres ) : ModelEstimator< PointWithNormal_3, double >( 3 ), m_sqDistanceThreshold( distThres*distThres ), m_angleThreshold( angleThres ) {}  

void PlaneWithNormalModelEstimator::fit(	const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data,
								std::vector< double > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}
	
	Plane_3 plane( (*data[0]).p, (*data[1]).p, (*data[2]).p ) ;

	parameters.push_back( plane.a() ) ;
	parameters.push_back( plane.b() ) ;
	parameters.push_back( plane.c() ) ;
	parameters.push_back( plane.d() ) ;
		
}



void PlaneWithNormalModelEstimator::leastSquaresFit(	const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data,
											std::vector< double > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	// Change Point_3* data vector into a Point_3 vector
	std::vector< Point_3 > dataRef ;
	std::vector< std::vector< PointWithNormal_3 >::const_pointer >::const_iterator itDataPtr ;
	for ( itDataPtr = data.begin(); itDataPtr != data.end(); ++itDataPtr ) {
		dataRef.push_back( (**itDataPtr).p ) ;
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


	
bool PlaneWithNormalModelEstimator::compatible(	const std::vector< double > &parameters, 
										const PointWithNormal_3 &data ) const {

	Plane_3 plane( parameters[0], parameters[1], parameters[2], parameters[3] ) ;
	double sqDistance = squared_distance( data.p, plane ) ;
	
	Vector_3 normal = normalizeVector( plane.orthogonal_vector() ) ;

	return sqDistance < ( m_sqDistanceThreshold ) && 
		compatibleVectors( data.n, normal ) &&
		dihedralAngleNormals( data.n, normal ) < m_angleThreshold ;

}



bool PlaneWithNormalModelEstimator::degenerate( const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}

	Plane_3 plane( (*data[0]).p, (*data[1]).p, (*data[2]).p ) ;
	
	//std::cout << "[RANSAC] First Point = " << data[0]->x() << data[0]->y() << data[0]->z() << std::endl ;
	//std::cout << "[RANSAC] Plane = " << plane.a() << " " << plane.b() << " " << plane.c() << " " << plane.d() << " " << std::endl ;
	/*if ( plane.is_degenerate() ) 
		std::cout << "[RANSAC] Degenerate." << std::endl ;
	else
		std::cout << "[RANSAC] Good sample set." << std::endl ;*/

	// Also compute compatibility between normals

	//std::cout << "[RANSAC::degenerate] points: " << std::endl 
	//											 << "  " << (*data[0]).p << std::endl << "  " << (*data[1]).p << std::endl << "  " << (*data[2]).p << std::endl ;
	//std::cout << "[RANSAC::degenerate] normals: " << std::endl 
	//											  << "  " << (*data[0]).n << std::endl << "  " << (*data[1]).n << std::endl << "  " << (*data[2]).n << std::endl ;
	//std::cout << "[RANSAC] Diedral angles (unoriented): " << dihedralAngleUnorientedNormals( (*data[0]).n, (*data[1]).n ) << ", " <<
	//			 dihedralAngleUnorientedNormals( (*data[0]).n, (*data[2]).n ) << ", " <<
	//			 dihedralAngleUnorientedNormals( (*data[1]).n, (*data[2]).n ) << std::endl ;

	return plane.is_degenerate() ||
		   !compatibleVectors( (*data[0]).n, (*data[1]).n ) ||
		   !compatibleVectors( (*data[0]).n, (*data[2]).n ) ||
		   !compatibleVectors( (*data[1]).n, (*data[2]).n ) ||
		   dihedralAngleNormals( (*data[0]).n, (*data[1]).n ) > m_angleThreshold ||
		   dihedralAngleNormals( (*data[0]).n, (*data[2]).n ) > m_angleThreshold ||
		   dihedralAngleNormals( (*data[1]).n, (*data[2]).n ) > m_angleThreshold ;

}



/* Auxiliary functions */
double PlaneWithNormalModelEstimator::angleWrap( double a ) const {
	if ( a > CGAL_PI ) {
		return a - ( 2*CGAL_PI ) ;
	}
	else if ( a < -CGAL_PI ) {
		return a + ( 2*CGAL_PI ) ;
	}
	else {
		return a ;
	}
}



double PlaneWithNormalModelEstimator::dihedralAngleNormals( const Vector_3& v1, const Vector_3& v2 ) const {
	
	Vector_3 crossProd = CGAL::cross_product( v1, v2 ) ;
	double normCrossProd = CGAL::sqrt( crossProd.squared_length() ) ;
	double dotProd = v1*v2 ;
	//return angleWrap( atan2( normCrossProd, dotProd ) ) ;
	return atan2( normCrossProd, dotProd ) ;

}



Vector_3 PlaneWithNormalModelEstimator::normalizeVector( const Vector_3& v ) const {

	double length = CGAL::sqrt( v.squared_length() ) ;
	return v/length ;

}