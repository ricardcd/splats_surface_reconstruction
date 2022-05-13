#ifndef CGAL_SPLAT_3_RAY_3_INTERSECTION_H
#define CGAL_SPLAT_3_RAY_3_INTERSECTION_H

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>

#include <iostream> // Debug

namespace CGAL {

namespace internal {

template <class K>
Object
intersection( const Splat_3<K> &s, 
			  const typename K::Ray_3 &r,
			  const K & k )
{
	typedef typename K::Point_3			Point_3 ;
	typedef typename K::Ray_3			Ray_3 ;
	
	// Compute the intersection between the ray and the Splat_3
	bool do_intersect = false ;
	Point_3 result ;
	if ( s.monge().getDegree() == 1 ) {
		Object int1 = s.monge().intersection_deg1( r ) ;
		if ( const Point_3 *ipoint = CGAL::object_cast< Point_3 >( &int1 ) ) {
			result = *ipoint ;
			do_intersect = true ;
		}
	}
	else {
		if ( s.monge().getDegree() == 2 ) {
			Object int1, int2 ;
			s.monge().intersection_deg2( r, int1, int2 ) ;
			if ( const Point_3 *ipoint1 = CGAL::object_cast< Point_3 >( &int1 ) ) {
				do_intersect = true ;
				if ( const Point_3 *ipoint2 = CGAL::object_cast< Point_3 >( &int2 ) ) {
					// Two intersections, select the one closest to the query point
					if ( CGAL::squared_distance( r.start(), *ipoint1 ) < CGAL::squared_distance( r.start(), *ipoint2 ) )
						result = *ipoint1 ;
					else
						result = *ipoint2 ;
				}
				else {
					result = *ipoint1 ;
				}
			}
			else {
				do_intersect = false ;
			}
		}
		else {
			double residuals[1] ;
			do_intersect = s.monge().intersection( r.start(), r.to_vector(), result, residuals ) ;
		}
	}

	if ( do_intersect ) {
		
		// Intersection is a point, check if its inside the circle
		if ( CGAL::squared_distance( result, s.center() ) < s.squared_radius() ) {
			return CGAL::make_object( result ) ;
		}
		else {
			//std::cout << "Intersection with the plane, but not inside circle..." << std::endl ;
			return Object() ;
		}
		
	}
	else {
		return Object() ;
	}
}

template <class K>
inline
Object
intersection( const typename K::Ray_3 &r, 
			  const Splat_3<K> &c,
			  const K & k )
{
	return internal::intersection( c, r ) ;
}

} // end namespace internal



template <class K>
inline
Object
intersection(const Splat_3<K> &c, const Ray_3<K> &r)
{
  // return typename K::Intersect_3()(c, r);
  return internal::intersection(c,r,K());
}

template <class K>
inline
Object
intersection(const Ray_3<K> &r, const Splat_3<K> &c)
{
  return internal::intersection(c,r,K());
}

} // end namespace CGAL

#endif // CGAL_SPLAT_3_RAY_3_INTERSECTION_H