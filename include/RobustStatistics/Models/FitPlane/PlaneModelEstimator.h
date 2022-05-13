#ifndef _PLANE_MODEL_ESTIMATOR_H_
#define _PLANE_MODEL_ESTIMATOR_H_

#include "RobustStatistics/ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_3		Point_3 ;
typedef K::Plane_3		Plane_3 ;


class PlaneModelEstimator : public ModelEstimator< Point_3, double > {

public:
	PlaneModelEstimator( double distThres ) ;

	virtual void fit( const std::vector< std::vector< Point_3 >::const_pointer > &data,
					  std::vector< double > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< std::vector< Point_3 >::const_pointer > &data,
								  std::vector< double > &parameters ) const ;
	
	virtual bool compatible( const std::vector< double > &parameters, 
							 const Point_3 &data ) const ;

	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
												const std::vector< Point_3 > &data ) const ;

	virtual bool degenerate( const std::vector< std::vector< Point_3 >::const_pointer > &data ) const ;

private:
	double m_sqDistanceThreshold ;
} ;


#endif // _PLANE_MODEL_ESTIMATOR_H_