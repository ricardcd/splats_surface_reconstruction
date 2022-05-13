#ifndef _LOCAL_SMOOTH_SURFACE_MODEL_ESTIMATOR_H_
#define _LOCAL_SMOOTH_SURFACE_MODEL_ESTIMATOR_H_

#include "RobustStatistics/ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
typedef CGAL::Monge_via_jet_fitting<K, K>			Monge_via_jet_fitting ;
typedef Monge_via_jet_fitting::Monge_form			Monge_form ;

typedef K::Point_3		Point_3 ;
typedef K::Vector_3		Vector_3 ;
typedef K::FT			FT ;


class LocalSmoothSurfaceModelEstimator : public ModelEstimator< Point_3, FT > {

public:
	LocalSmoothSurfaceModelEstimator( unsigned int minSamples, unsigned int degree, Point_3 refPoint, FT distThres ) ;

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

private:
	FT m_sqDistanceThreshold ;
	unsigned int m_degree ;
	Point_3 m_refPoint ;

} ;


#endif // _LOCAL_SMOOTH_SURFACE_MODEL_ESTIMATOR_H_