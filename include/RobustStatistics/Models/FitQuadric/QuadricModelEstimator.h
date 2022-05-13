#ifndef _QUADRIC_MODEL_ESTIMATOR_H_
#define _QUADRIC_MODEL_ESTIMATOR_H_

#include "ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

#include <CGAL/Monge_via_jet_fitting.h> // Remember to include the path to the modified version first during compilation
typedef CGAL::Monge_via_jet_fitting<K, K>			Monge_via_jet_fitting ;
typedef Monge_via_jet_fitting::Monge_form			Monge_form ;

typedef K::Point_3		Point_3 ;
typedef K::Vector_3		Vector_3 ;
typedef K::FT			FT ;


class QuadricModelEstimator : public ModelEstimator< Point_3, FT > {

public:
	QuadricModelEstimator( unsigned int minInliers, FT distThres ) ;

	virtual void fit( const std::vector< std::vector< Point_3 >::const_pointer > &data,
					  std::vector< FT > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< std::vector< Point_3 >::const_pointer > &data,
								  std::vector< FT > &parameters ) const ;
	
	virtual bool compatible( const std::vector< FT > &parameters, 
							 const Point_3 &data ) const ;

	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
												const std::vector< Point_3 > &data ) const {
		// TODO: Implement!!
		return std::vector<double>() ;
	}

	virtual bool degenerate( const std::vector< std::vector< Point_3 >::const_pointer > &data ) const ;

private:
	double m_sqDistanceThreshold ;
} ;


#endif // _QUADRIC_MODEL_ESTIMATOR_H_