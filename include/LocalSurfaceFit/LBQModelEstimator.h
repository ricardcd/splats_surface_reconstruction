/*
Local Bivariate Quadric estimator
*/

#ifndef LBQMODELESTIMATOR_H
#define LBQMODELESTIMATOR_H

#include "RobustStatistics/ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

#include "LocalBivariateQuadric.h"

typedef K::Point_3		Point_3 ;
typedef K::Plane_3		Plane_3 ;
typedef K::Vector_3		Vector_3 ;
typedef K::FT			FT ;

class LBQModelEstimator : public ModelEstimator< Point_3, FT > {

	typedef Fitting::LocalBivariateQuadric<K>		LocalBivariateQuadric ;

public:
	LBQModelEstimator( const FT &distThres ) ;

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

} ;

#endif // LBQSEGMENTMODELESTIMATOR_H