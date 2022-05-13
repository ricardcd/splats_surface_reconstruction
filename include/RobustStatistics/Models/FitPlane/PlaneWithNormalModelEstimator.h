#ifndef _PLANE_WITH_NORMAL_MODEL_ESTIMATOR_H_
#define _PLANE_WITH_NORMAL_MODEL_ESTIMATOR_H_

#include "ModelEstimator.h"
#include <iostream> // Debug purposes

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_3		Point_3 ;
typedef K::Plane_3		Plane_3 ;
typedef K::Vector_3		Vector_3 ;

// Create struct defining the point and the normal
struct PointWithNormal_3 {
	Point_3 p ;
	Vector_3 n ;
} ;



class PlaneWithNormalModelEstimator : public ModelEstimator< PointWithNormal_3, double > {

public:
	PlaneWithNormalModelEstimator( double distThres, double angleThres ) ;

	virtual void fit( const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data,
					  std::vector< double > &parameters ) const ;

	virtual void leastSquaresFit( const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data,
								  std::vector< double > &parameters ) const ;
	
	virtual bool compatible( const std::vector< double > &parameters, 
							 const PointWithNormal_3 &data ) const ;

	virtual std::vector< double > sqResiduals(  const std::vector< double > &parameters, 
												const std::vector< PointWithNormal_3 > &data ) const {
		// TODO: Implement!!
		return std::vector<double>() ;
	}

	virtual bool degenerate( const std::vector< std::vector< PointWithNormal_3 >::const_pointer > &data ) const ;

	inline bool compatibleVectors( const Vector_3& v1, const Vector_3& v2 ) const { return ( v1*v2 ) >= 0 ; }

	double angleWrap( double a ) const ;

	double dihedralAngleNormals( const Vector_3& v1, const Vector_3& v2 ) const ;
			
	Vector_3 normalizeVector( const Vector_3& v ) const ;

private:
	double m_sqDistanceThreshold ;
	double m_angleThreshold ;

} ;


#endif // _PLANE_WITH_NORMAL_MODEL_ESTIMATOR_H_