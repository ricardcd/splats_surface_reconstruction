/*
Local Bivariate Quadric estimator (using weighted Least squares), restricted on segment

*/

#ifndef WLBQSEGMENTMODELESTIMATOR_H
#define WLBQSEGMENTMODELESTIMATOR_H

#include "RobustStatistics/ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/centroid.h>
#include "LocalBivariateQuadric.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel		K ;
typedef K::Point_3												Point_3 ;
typedef K::Vector_3												Vector_3 ;
typedef K::Segment_3											Segment_3 ;
typedef K::Line_3												Line_3 ;
typedef K::FT													FT ;
typedef Fitting::LocalBivariateQuadric<K>						LocalBivariateQuadric ;



class WLBQSegmentModelEstimator : public ModelEstimator< Point_3, FT > {
	
public:
	WLBQSegmentModelEstimator( const Segment_3 &segment, const FT &distThres, const FT &gaussianH ) ;

	virtual void fit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
					  std::vector< FT > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
								  std::vector< FT > &parameters ) const ;
	
	virtual bool compatible( const std::vector< FT > &parameters, 
							 const Point_3 &data ) const ;

	// Used in LeastKthSquares method	
	virtual std::vector< double > sqResiduals(  const std::vector< FT > &parameters, 
											  const std::vector< Point_3 > &data ) const ;

	virtual bool degenerate( const std::vector< typename std::vector< Point_3 >::const_pointer > &data ) const ;

	void setSqDistanceThreshold( double sqDistanceThreshold ) { m_sqDistanceThreshold = sqDistanceThreshold ; }

	static LocalBivariateQuadric	parameters2LBQ( const std::vector< FT > &parameters ) ;
	static std::vector< FT >		LBQ2parameters( const LocalBivariateQuadric& lbq ) ;

private:
	FT m_sqDistanceThreshold ;
	Segment_3 m_segment ;
	FT m_gaussianH ;

} ;

#endif // WLBQSEGMENTMODELESTIMATOR_H