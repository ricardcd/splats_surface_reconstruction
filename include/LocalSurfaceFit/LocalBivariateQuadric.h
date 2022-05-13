#ifndef LOCALBIVARIATEQUADRIC_H
#define LOCALBIVARIATEQUADRIC_H

// Eigen-related includes
#undef Success
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/SVD>

//CGAL
#include <CGAL/Object.h>

// Std
#include <iostream>
#include <fstream>



namespace Fitting {

template <class K > // K stands for Kernel...
class LocalBivariateQuadric {

	typedef typename K::FT							FT ;
	typedef typename K::Point_3						Point_3 ;
	typedef typename K::Vector_3					Vector_3 ;
	typedef typename K::Plane_3						Plane_3 ;
	typedef typename K::Segment_3					Segment_3 ;
	typedef typename K::Ray_3						Ray_3 ;
	typedef Eigen::Matrix< FT, 
						   Eigen::Dynamic, 
						   1 >						EigenVector ;
	typedef Eigen::Matrix< FT,
						   ::Eigen::Dynamic,
						   ::Eigen::Dynamic>		EigenMatrix ;
	typedef CGAL::Aff_transformation_3<K>			Aff_transformation ;
	typedef CGAL::Object							Object ;
	


	// Class Attributes
	Plane_3								m_fitPlane ;	
	Vector_3							m_xAxis, m_yAxis, m_zAxis ;
	Point_3								m_origin ;
	std::vector<FT>						m_coefficients ;
	Aff_transformation					m_world2fit, m_fit2world ;



public:

	// Default constructor
	LocalBivariateQuadric() {
		m_fitPlane = Plane_3() ;	
		m_xAxis = Vector_3() ;
		m_yAxis = Vector_3() ;
		m_zAxis = Vector_3() ;
		m_origin = Point_3() ;
		m_coefficients = std::vector<FT>() ;
		m_world2fit = Aff_transformation() ;
		m_fit2world = Aff_transformation() ;
	}

	// Constructor
	// Define the local fitting plane using the provided normal vector
	LocalBivariateQuadric( const Point_3& origin, const Vector_3& normal ) {
		m_zAxis = normal ;
		m_origin = origin ;
		m_fitPlane = Plane_3( origin, normal ) ;
		m_xAxis = m_fitPlane.base1() ;
		m_yAxis = m_fitPlane.base2() ;

		// Unitize axis?
		if ( m_zAxis.squared_length() != 1 ) {
			// Unitize
			m_zAxis = m_zAxis / CGAL::sqrt( m_zAxis.squared_length() ) ;
		}
		if ( m_xAxis.squared_length() != 1 ) {
			// Unitize
			m_xAxis = m_xAxis / CGAL::sqrt( m_xAxis.squared_length() ) ;
		}
		if ( m_yAxis.squared_length() != 1 ) {
			// Unitize
			m_yAxis = m_yAxis / CGAL::sqrt( m_yAxis.squared_length() ) ;
		}		

		// Compute transformations to change basis World <--> Fit
		Aff_transformation rotation( m_xAxis.x(), m_xAxis.y(), m_xAxis.z(), 0, 
									 m_yAxis.x(), m_yAxis.y(), m_yAxis.z(), 0, 
									 m_zAxis.x(), m_zAxis.y(), m_zAxis.z(), 0 ) ;
	
		Point_3 orig( 0.0, 0.0, 0.0 ) ;
		Vector_3 v( orig - m_origin ) ;
		Aff_transformation translation( CGAL::TRANSLATION, v ) ;
		
		m_world2fit = rotation * translation ;
		m_fit2world = m_world2fit.inverse() ;

		// Coefficients will be fit later
		m_coefficients.resize( 6 ) ;
		m_coefficients[0] = m_coefficients[1] = m_coefficients[2] = m_coefficients[3] = m_coefficients[4] = m_coefficients[5] = 0.0 ;
	}



	// Constructor with all the information defining a LBQ
	LocalBivariateQuadric( const Point_3& origin, const Vector_3& normal, const Vector_3& xAxis, const Vector_3& yAxis, const std::vector<FT>& coeffs ) {
		m_zAxis = normal ;
		m_origin = origin ;
		m_fitPlane = Plane_3( origin, normal ) ;
		m_xAxis = xAxis ;
		m_yAxis = yAxis ;

		// Unitize axis?
		if ( m_zAxis.squared_length() != 1 ) {
			// Unitize
			m_zAxis = m_zAxis / CGAL::sqrt( m_zAxis.squared_length() ) ;
		}
		if ( m_xAxis.squared_length() != 1 ) {
			// Unitize
			m_xAxis = m_xAxis / CGAL::sqrt( m_xAxis.squared_length() ) ;
		}
		if ( m_yAxis.squared_length() != 1 ) {
			// Unitize
			m_yAxis = m_yAxis / CGAL::sqrt( m_yAxis.squared_length() ) ;
		}	

		// WARNING: No orthogonality check is performed on the reference system axis
		// Compute transformations to change basis World <--> Fit
		Aff_transformation rotation( m_xAxis.x(), m_xAxis.y(), m_xAxis.z(), 0, 
									 m_yAxis.x(), m_yAxis.y(), m_yAxis.z(), 0, 
									 m_zAxis.x(), m_zAxis.y(), m_zAxis.z(), 0 ) ;
	
		Point_3 orig( 0.0, 0.0, 0.0 ) ;
		Vector_3 v( orig - m_origin ) ;
		Aff_transformation translation( CGAL::TRANSLATION, v ) ;
		
		m_world2fit = rotation * translation ;
		m_fit2world = m_world2fit.inverse() ;

		// Coefficients
		m_coefficients = coeffs ;

	}



	// Getters
	Point_3 origin() const { return m_origin ; }
	Vector_3 axisX() const { return m_xAxis ; }
	Vector_3 axisY() const { return m_yAxis ; }
	Vector_3 axisZ() const { return m_zAxis ; }		
	std::vector<FT> coefficients() const { return m_coefficients; }
	void print( std::ofstream &os ) const { 
		os << m_origin << " " << m_xAxis << " " << m_yAxis << " " << m_zAxis << " " << m_coefficients[ 0 ] << " " << m_coefficients[ 1 ] << " " << m_coefficients[ 2 ] << " " << m_coefficients[ 3 ] << " " << m_coefficients[ 4 ] << " " << m_coefficients[ 5 ] ;
	}
	void show() const { 
		std::cout << m_origin << " " << m_xAxis << " " << m_yAxis << " " << m_zAxis << " " << m_coefficients[ 0 ] << " " << m_coefficients[ 1 ] << " " << m_coefficients[ 2 ] << " " << m_coefficients[ 3 ] << " " << m_coefficients[ 4 ] << " " << m_coefficients[ 5 ] << std::endl ;
	}
	bool valid() const { return m_coefficients.size() > 0 ; }



	// Fitting a quadric to a set of points (in world coordinates!)
	FT fit( const std::vector< Point_3 >& points ) {

		int numPts = points.size() ;

		// Transform the points to the local fitting coordinate system
		std::vector< Point_3 > pointsFit ;
		typename std::vector< Point_3 >::const_iterator itP ;
		for ( itP = points.begin(); itP != points.end(); ++itP ) {
			pointsFit.push_back( m_world2fit( *itP ) ) ;
		}

		// Build system of equations Ax=B
		EigenVector B( numPts * 6 ) ;
		EigenMatrix A( numPts * 6, 6 ) ;
		int i = 0 ;
		for ( itP = pointsFit.begin(); itP != pointsFit.end(); ++itP, i++ ) {
			FT x = itP->x() ;
			FT y = itP->y() ;
			FT z = itP->z() ;

			FT x2 = x*x ;
			FT x3 = x2*x ;
			FT x4 = x3*x ;
			FT y2 = y*y ;
			FT y3 = y2*y ;
			FT y4 = y3*y ;

			A( 6*i, 0 )	  = x4 ;			A( 6*i, 1 )   = x2*y2 ;			A( 6*i, 2 )   = x3 ;		A( 6*i, 3 )   = x2*y ;		A( 6*i, 4 )   = x3*y ;			A( 6*i, 5 )   = x2 ;
			A( 6*i+1, 0 ) = x2*y2 ;			A( 6*i+1, 1 ) = y4 ;			A( 6*i+1, 2 ) = x*y2 ;		A( 6*i+1, 3 ) = y3 ;		A( 6*i+1, 4 ) = x*y3 ;			A( 6*i+1, 5 ) = y2 ;
			A( 6*i+2, 0 ) = x3 ;			A( 6*i+2, 1 ) = x*y2 ;			A( 6*i+2, 2 ) = x2 ;		A( 6*i+2, 3 ) = x*y ;		A( 6*i+2, 4 ) = x2*y ;			A( 6*i+2, 5 ) = x ;
			A( 6*i+3, 0 ) = x2*y ;			A( 6*i+3, 1 ) = y3 ;			A( 6*i+3, 2 ) = x*y ;		A( 6*i+3, 3 ) = y2 ;		A( 6*i+3, 4 ) = x*y2 ;			A( 6*i+3, 5 ) = y ;
			A( 6*i+4, 0 ) = x3*y ;			A( 6*i+4, 1 ) = x*y3 ;			A( 6*i+4, 2 ) = x2*y ;		A( 6*i+4, 3 ) = x*y2 ;		A( 6*i+4, 4 ) = x2*y2 ;			A( 6*i+4, 5 ) = x*y ;
			A( 6*i+5, 0 ) = x2 ;			A( 6*i+5, 1 ) = y2 ;			A( 6*i+5, 2 ) = x ;			A( 6*i+5, 3 ) = y ;			A( 6*i+5, 4 ) = x*y ;			A( 6*i+5, 5 ) = 1 ;

			B( 6*i ) = z*x2 ;				B( 6*i+1 ) = z*y2 ;				B( 6*i+2 ) = z*x ;			B( 6*i+3 ) = z*y ;			B( 6*i+4 ) = z*x*y ;			B( 6*i+5 ) = z ;
		}

		// Solve the system of equations
		Eigen::JacobiSVD< EigenMatrix > jacobiSvd( A,::Eigen::ComputeThinU | ::Eigen::ComputeThinV ) ;
		EigenVector coeffs( 6 ) ;
		coeffs = jacobiSvd.solve( EigenVector( B ) ) ;

		m_coefficients[0] = coeffs(0) ;
		m_coefficients[1] = coeffs(1) ;
		m_coefficients[2] = coeffs(2) ;
		m_coefficients[3] = coeffs(3) ;
		m_coefficients[4] = coeffs(4) ;
		m_coefficients[5] = coeffs(5) ;

		// Return condition number
		return jacobiSvd.singularValues().array().abs().maxCoeff() /
			   jacobiSvd.singularValues().array().abs().minCoeff() ;
				
	}



	// Fitting a quadric to a set of points
	bool weightedFit( const std::vector< Point_3 >& points, const std::vector< FT >& weights ) {

		// There must be a weight specified for each input point
		if ( points.size() != weights.size() ) {
			return false ;
		}

		int numPts = points.size() ;

		// Transform the points to the local fitting coordinate system
		std::vector< Point_3 > pointsFit ;
		typename std::vector< Point_3 >::const_iterator itP ;
		for ( itP = points.begin(); itP != points.end(); ++itP ) {
			pointsFit.push_back( m_world2fit( *itP ) ) ;
		}

		// Build system of equations Ax=B
		EigenVector B( numPts * 6 ) ;
		EigenMatrix A( numPts * 6, 6 ) ;
		int i = 0 ;
		for ( itP = pointsFit.begin(); itP != pointsFit.end(); ++itP, i++ ) {
			FT x = itP->x() ;
			FT y = itP->y() ;
			FT z = itP->z() ;

			FT x2 = x*x ;
			FT x3 = x2*x ;
			FT x4 = x3*x ;
			FT y2 = y*y ;
			FT y3 = y2*y ;
			FT y4 = y3*y ;

			// Apply weight of current observation
			FT w = weights[i] ;
			
			A( 6*i, 0 )	  = w*x4 ;			A( 6*i, 1 )   = w*x2*y2 ;		A( 6*i, 2 )   = w*x3 ;			A( 6*i, 3 )   = w*x2*y ;		A( 6*i, 4 )   = w*x3*y ;		A( 6*i, 5 )   = w*x2 ;
			A( 6*i+1, 0 ) = w*x2*y2 ;		A( 6*i+1, 1 ) = w*y4 ;			A( 6*i+1, 2 ) = w*x*y2 ;		A( 6*i+1, 3 ) = w*y3 ;			A( 6*i+1, 4 ) = w*x*y3 ;		A( 6*i+1, 5 ) = w*y2 ;
			A( 6*i+2, 0 ) = w*x3 ;			A( 6*i+2, 1 ) = w*x*y2 ;		A( 6*i+2, 2 ) = w*x2 ;			A( 6*i+2, 3 ) = w*x*y ;			A( 6*i+2, 4 ) = w*x2*y ;		A( 6*i+2, 5 ) = w*x ;
			A( 6*i+3, 0 ) = w*x2*y ;		A( 6*i+3, 1 ) = w*y3 ;			A( 6*i+3, 2 ) = w*x*y ;			A( 6*i+3, 3 ) = w*y2 ;			A( 6*i+3, 4 ) = w*x*y2 ;		A( 6*i+3, 5 ) = w*y ;
			A( 6*i+4, 0 ) = w*x3*y ;		A( 6*i+4, 1 ) = w*x*y3 ;		A( 6*i+4, 2 ) = w*x2*y ;		A( 6*i+4, 3 ) = w*x*y2 ;		A( 6*i+4, 4 ) = w*x2*y2 ;		A( 6*i+4, 5 ) = w*x*y ;
			A( 6*i+5, 0 ) = w*x2 ;			A( 6*i+5, 1 ) = w*y2 ;			A( 6*i+5, 2 ) = w*x ;			A( 6*i+5, 3 ) = w*y ;			A( 6*i+5, 4 ) = w*x*y ;			A( 6*i+5, 5 ) = w ;

			B( 6*i ) = w*z*x2 ;				B( 6*i+1 ) = w*z*y2 ;			B( 6*i+2 ) = w*z*x ;			B( 6*i+3 ) = w*z*y ;			B( 6*i+4 ) = w*z*x*y ;			B( 6*i+5 ) = w*z ;

		}

		// Solve the system of equations
		Eigen::JacobiSVD< EigenMatrix > jacobiSvd( A,::Eigen::ComputeThinU | ::Eigen::ComputeThinV ) ;
		EigenVector coeffs( 6 ) ;
		coeffs = jacobiSvd.solve( EigenVector( B ) ) ;

		m_coefficients[0] = coeffs(0) ;
		m_coefficients[1] = coeffs(1) ;
		m_coefficients[2] = coeffs(2) ;
		m_coefficients[3] = coeffs(3) ;
		m_coefficients[4] = coeffs(4) ;
		m_coefficients[5] = coeffs(5) ;

		// Return condition number
		return jacobiSvd.singularValues().array().abs().maxCoeff() /
			   jacobiSvd.singularValues().array().abs().minCoeff() ;


	}



	// Gaussian weight function
	static FT weightGaussian( const FT& dist, const FT& h = 1.0 ) {
		return exp( -( ( dist*dist ) / ( h*h ) ) ) ;
	}



	// Wendland weight function
	static FT weightWendland( const FT& dist, const FT& h = 1.0 ) {
		return ( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*(4*dist / h + 1 ) ;
	}



	Point_3 projection( const Point_3& pt ) const {

		// Change to fitting basis
		Point_3 ptFit = m_world2fit( pt ) ;
		
		// Evaluate the function
		FT x = ptFit.x() ;
		FT y = ptFit.y() ;
		FT z = x*x*m_coefficients[0] + y*y*m_coefficients[1] + x*m_coefficients[2] + y*m_coefficients[3] + x*y*m_coefficients[4] + m_coefficients[5] ;

		// Change to world basis
		return m_fit2world( Point_3( x, y, z ) ) ;

	}



	// Returns the intersection points with a ray
	bool intersection( Ray_3 ray, Object& int1, Object& int2 ) {
		
		// get normalized direction vector
		Vector_3 dir = ray.to_vector() ;
		if ( dir.squared_length() != 1.0 ) {		
			dir = dir / CGAL::sqrt( dir.squared_length() ) ;
		}

		/* Solve the closed form */
		
		// Change basis (move to fitting basis)
		Point_3 startRayFit = m_world2fit( ray.source() ) ;
		Vector_3 dirRayFit = m_world2fit( dir ) ;

		// Elements of the general equation
		FT x0 = startRayFit.x() ;
		FT y0 = startRayFit.y() ;
		FT z0 = startRayFit.z() ;
		FT xd = dirRayFit.x() ;
		FT yd = dirRayFit.y() ;
		FT zd = dirRayFit.z() ;
		std::vector< FT > coeffs = coefficients() ;
		FT a = coeffs[0] ;
		FT b = coeffs[1] ;
		FT c = coeffs[2] ;
		FT d = coeffs[3] ;
		FT e = coeffs[4] ;
		FT f = coeffs[5] ;

		// --- Debug (Start) ---
		//std::cout << "startRay(fit) = " << startRayFit << std::endl ;
		//std::cout << "dirRay(fit) = " << dirRayFit << std::endl ;
		//std::cout << "a = " << a << std::endl ;
		//std::cout << "b = " << b << std::endl ;
		//std::cout << "c = " << c << std::endl ;
		//std::cout << "d = " << d << std::endl ;
		//std::cout << "e = " << e << std::endl ;
		//std::cout << "f = " << f << std::endl ;
		// --- Debug  (End)  ---
						
		// A, B and C coefficients of general quadratic equation Ax + By + C = 0
		FT A = a*xd*xd + b*yd*yd + e*xd*yd ;
		FT B = 2*a*x0*xd + 2*b*y0*yd + c*xd + d*yd + e*x0*yd + e*xd*y0 - zd ;
		FT C = a*x0*x0 + b*y0*y0 + c*x0 + d*y0 + e*x0*y0 + f - z0 ;
		
		// --- Debug (Start) ---
		//std::cout << "A = " << A << std::endl ;
		//std::cout << "B = " << B << std::endl ;
		//std::cout << "C = " << C << std::endl ;
		// --- Debug  (End)  ---

		// Get solutions
		FT discriminant = (B*B) - 4*A*C ;
		
		if ( discriminant < 0 ) {
			// Result is an imaginary number, no intersection in the real plane

			int1 = Object() ;
			int2 = Object() ;
			return false ;
		}
		else {

			// Ray intersects the surface at exactly one point (tangent to the quadric) OR the ray intersects the sphere at 2 points 
			// (get the minimum, since it correspond to the closest intersection) 

			// Check if A == 0, then is a linear equation!
			FT t1, t2 ;
			if ( abs( A ) < 1e-20 ) {
				t1 = -C/B ;
				t2 = -999999.9 ; // dummy, just a single intersection
			}
			else {
				t1 = ( -B + sqrt( discriminant ) ) / ( 2*A ) ;
				t2 = ( -B - sqrt( discriminant ) ) / ( 2*A ) ;
			}

			// --- Debug (Start) ---
			/*std::cout << "discriminant = " << discriminant << std::endl ;
			std::cout << "t1 = " << t1 << std::endl ;
			std::cout << "t2 = " << t2 << std::endl ;*/
			// --- Debug  (End)  ---

			if ( t1 >= 0 ) {
				// t1 solution is inside the ray
				int1 = CGAL::make_object( m_fit2world( startRayFit + ( dirRayFit * t1 ) ) ) ;
				if ( t2 >= 0 ) {
					// t2 solution is inside the ray
					int2 = CGAL::make_object( m_fit2world( startRayFit + ( dirRayFit * t2 ) ) ) ;
					return true ;
				}
				else {
					// Second solution is outside the ray, do not return it
					int2 = Object() ;
					return true ;
				}
			}
			else {
				// t1 solution is NOT inside the ray, check t2
				int2 = Object() ; // Solution 2 will be left blank, if there is a solution will be placed in "int1"

				if ( t2 >= 0 ) {
					// t2 solution is inside the ray
					int1 = CGAL::make_object( m_fit2world( startRayFit + ( dirRayFit * t2 ) ) ) ;					
					return true ;
				}
				else {
					// Second solution is outside the ray, do not return it
					int1 = Object() ;					
					return false ;
				}
			}		
		}

	}
	


	// Returns the intersection points with a segment
	bool intersection( const Segment_3& seg, Object& int1, Object& int2 ) {

		// Compute the ray-lbq intersection
		Ray_3 ray( seg.start(), seg.to_vector() ) ;
		intersection( ray, int1, int2 ) ;

		// Check the validity of the results
		Point_3 dummy ;
		if ( CGAL::assign( dummy, int1 ) ) {
			// Check that is within the range
			if ( CGAL::squared_distance( dummy, seg.start() ) > seg.squared_length() ) {
				// "int1" is NOT inside the segment, check "int2"
				if ( CGAL::assign( dummy, int2 ) ) {
					if ( CGAL::squared_distance( dummy, seg.start() ) <= seg.squared_length() ) {
						// int2 is within the range, save it in place of "int1"
						int1 = int2 ;
						int2 = Object() ;
						return true ;
					}
					else {
						int1 = Object() ;
						int2 = Object() ;
						return false ;
					}
				}
				else {
					int1 = Object() ;
					int2 = Object() ;
					return false ;
				}
			}
			else {
				// "int1" is inside the segment, check "int2"
				if (  CGAL::assign( dummy, int2 ) && 
					( CGAL::squared_distance( dummy, seg.start() ) > seg.squared_length() ) ) {
					// Outside segment
					int2 = Object() ;					
				}
				return true ;
			}
		}
		else {
			// Else, no intersection found (single intersections are stored in "int1")
			return false ;
		}

	}



} ;

} // namespace Fitting

#endif // LOCALBIVARIATEQUADRIC_H