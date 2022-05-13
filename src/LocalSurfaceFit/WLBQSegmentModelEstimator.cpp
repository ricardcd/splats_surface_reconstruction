#include "LocalSurfaceFit/WLBQSegmentModelEstimator.h"

WLBQSegmentModelEstimator::WLBQSegmentModelEstimator( const Segment_3& segment, 
													  const FT &distThres,
													  const FT &gaussianH ) : 
ModelEstimator< Point_3, FT >( 6 ), 
m_sqDistanceThreshold( distThres*distThres ), 
m_segment( segment ),
m_gaussianH( gaussianH )
{}

void WLBQSegmentModelEstimator::fit(	const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
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

	// Compute centroid, to fix the fitting basis along the segment
	Point_3 centroid = CGAL::centroid( dataRef.begin(), dataRef.end(), CGAL::Dimension_tag<0>() ) ;	

	// --- Debug (Start) ---
	//std::cout << "Segment = " << m_segment << std::endl ;
	//std::cout << "Centroid before projection" << centroid << std::endl ;
	//Point_3 cprev = centroid ;
	// --- Debug  (End)  ---

	// Project centroid to the fitting line
	centroid = m_segment.supporting_line().projection( centroid ) ; // WARNING: Should test if falls in the segment?? or it does not matter??
	
	// --- Debug (Start) ---
	//std::cout << "Centroid after projection" << centroid << std::endl ;
	//FT dot = Vector_3( centroid, cprev ) * m_segment.to_vector() ;
	//std::cout << "Dot product = " << dot << std::endl ;
	// --- Debug  (End)  ---

	// Compute weights according to the points' distance to the segment
	std::vector< FT > weights ;
	for( std::vector< Point_3 >::iterator itP = dataRef.begin(); itP != dataRef.end(); ++itP ) {
		// Distance from segment
		FT w = LocalBivariateQuadric::weightGaussian( CGAL::sqrt( CGAL::squared_distance( *itP, m_segment ) ),
													  m_gaussianH ) ;
		// Distance from projected centroid
		/*FT w = LocalBivariateQuadric::weightGaussian( CGAL::sqrt( CGAL::squared_distance( *itP, centroid ) ),
													  m_gaussianH ) ;*/
		//std::cout << "Dist = " << CGAL::sqrt( CGAL::squared_distance( *itP, m_segment ) ) << ", w = " << w << std::endl ;
		weights.push_back( w ) ;
	}
	// std::cout << std::endl ; 

	// Compute the LBQ on the computed plane
	LocalBivariateQuadric lbq( centroid, m_segment.to_vector() ) ;
	lbq.weightedFit( dataRef, weights ) ;

	parameters = LBQ2parameters( lbq ) ;

}



void WLBQSegmentModelEstimator::leastSquaresFit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
								  std::vector< FT > &parameters ) const {

	// Fitting is already a least squares process...
	fit( data, parameters ) ;	

}



bool WLBQSegmentModelEstimator::compatible(	const std::vector< double > &parameters, 
									const Point_3 &data ) const {

	// Rebuild the LBQ from the parameters
	LocalBivariateQuadric lbq = parameters2LBQ( parameters ) ;

	// Projection!!
	FT sqDistance = CGAL::squared_distance( data, lbq.projection( data ) ) ;

	return sqDistance < m_sqDistanceThreshold ;

}



std::vector< double > WLBQSegmentModelEstimator::sqResiduals( const std::vector< FT > &parameters, 
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



bool WLBQSegmentModelEstimator::degenerate( const std::vector< std::vector< Point_3 >::const_pointer > &data ) const {

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return true ;
	}

	return false ;

}



// Rebuild the LBQ from parameters vector
LocalBivariateQuadric	WLBQSegmentModelEstimator::parameters2LBQ( const std::vector< FT > &parameters ) {
		
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
std::vector< FT > WLBQSegmentModelEstimator::LBQ2parameters( const LocalBivariateQuadric& lbq ) {

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
