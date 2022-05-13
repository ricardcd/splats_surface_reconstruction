#ifndef _CENTROID_MODEL_ESTIMATOR_H_
#define _CENTROID_MODEL_ESTIMATOR_H_

#include "RobustStatistics/ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_3		Point_3 ;

class CentroidModelEstimator : public ModelEstimator< Point_3, double > {

public:
	CentroidModelEstimator( double distThres ) ;

	virtual void fit( const std::vector< std::vector<Point_3>::const_pointer > &data,
					  std::vector< double > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< std::vector<Point_3>::const_pointer > &data,
								  std::vector< double > &parameters ) const ;
	
	// Used in RANSAC method
	virtual bool compatible( const std::vector< double > &parameters, 
							 const Point_3 &data ) const ;

	// Used in LeastKthSquares method
	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
											  const std::vector< Point_3 > &data ) const ;

	virtual bool degenerate( const std::vector< std::vector<Point_3>::const_pointer > &data ) const ;

private:
	double m_distanceThreshold ;
} ;


#endif // _CENTROID_MODEL_ESTIMATOR_H_