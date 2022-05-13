#ifndef _LINE_MODEL_ESTIMATOR_H_
#define _LINE_MODEL_ESTIMATOR_H_

#include "ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2		Point_2 ;
typedef K::Line_2		Line_2 ;


class LineModelEstimator : public ModelEstimator< Point_2, double > {

public:
	LineModelEstimator( double distThres ) ;

	virtual void fit( const std::vector< std::vector<Point_2>::const_pointer > &data,
					  std::vector< double > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< std::vector<Point_2>::const_pointer > &data,
								  std::vector< double > &parameters ) const ;
	
	// Used in RANSAC method
	virtual bool compatible( const std::vector< double > &parameters, 
							 const Point_2 &data ) const ;

	// Used in LeastKthSquares method
	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
											  const std::vector< Point_2 > &data ) const ;

	virtual bool degenerate( const std::vector< std::vector<Point_2>::const_pointer > &data ) const ;

private:
	double m_distanceThreshold ;
} ;


#endif // _LINE_MODEL_ESTIMATOR_H_