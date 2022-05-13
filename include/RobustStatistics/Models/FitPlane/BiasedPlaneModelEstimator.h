#ifndef _BIASED_PLANE_MODEL_ESTIMATOR_H_
#define _BIASED_PLANE_MODEL_ESTIMATOR_H_

#include "ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_3		Point_3 ;
typedef K::Plane_3		Plane_3 ;


class BiasedPlaneModelEstimator : public ModelEstimator< Point_3, double > {

public:
	BiasedPlaneModelEstimator( double distThres, Point_3 forcedInlier ) ;

	virtual void fit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
					  std::vector< double > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< typename std::vector< Point_3 >::const_pointer > &data,
								  std::vector< double > &parameters ) const ;
	
	// Used in RANSAC method
	virtual bool compatible( const std::vector< double > &parameters, 
							 const Point_3 &data ) const ;

	// Used in LeastKthSquares method
	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
												const std::vector< Point_3 > &data ) const {
		// TODO: Implement!!
		return std::vector<double>() ;
	}

	virtual bool degenerate( const std::vector< typename std::vector< Point_3 >::const_pointer > &data ) const ;

private:
	double m_sqDistanceThreshold ;
	Point_3 m_forcedInlier ;

} ;


#endif // _BIASED_PLANE_MODEL_ESTIMATOR_H_