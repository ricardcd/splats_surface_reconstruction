#ifndef CGAL_SPLAT_3_SEGMENT_3_INTERSECTION_H
#define CGAL_SPLAT_3_SEGMENT_3_INTERSECTION_H

#include <iostream>

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include "Splat_3/Splat_3.h"
#include <math.h>



namespace CGAL {

namespace internal {

template <class K>
Object
intersection( const Splat_3<K> &splat,
			  const typename K::Segment_3 &s,
			  const K& k )
{
	typedef typename K::Point_3			Point_3 ;
	typedef typename K::Segment_3		Segment_3 ;
	typedef typename K::FT				FT ;

	
	// Compute the intersection between the plane and the line
	// Point_3 result ;
	// double residuals[1] ;
	// bool do_intersect = splat.monge().intersection( s.start(), s.to_vector(), result, residuals, 0.0, CGAL::sqrt<FT>( s.squared_length() ) ) ;
	bool do_intersect = false ;
	Point_3 result ;
	if ( splat.monge().getDegree() == 2 ) {
		Object int1, int2 ;
		splat.monge().intersection_deg2( s, int1, int2 ) ;
		if ( const Point_3 *ipoint1 = CGAL::object_cast< Point_3 >( &int1 ) ) {
			do_intersect = true ;
			if ( const Point_3 *ipoint2 = CGAL::object_cast< Point_3 >( &int2 ) ) {
				// Two intersections, select the one closest to the query point
				if ( CGAL::squared_distance( s.start(), *ipoint1 ) < CGAL::squared_distance( s.start(), *ipoint2 ) )
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
		do_intersect = splat.monge().intersection( s.start(), s.to_vector(), result, residuals, 0.0, CGAL::sqrt<FT>( s.squared_length() ) ) ;
	}

	if ( do_intersect ) {
		
		// Intersection is a point, check if its inside the circle
		if ( CGAL::squared_distance( result, splat.center() ) < splat.squared_radius() ) {
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
intersection( const typename K::Segment_3 &s, 
			  const Splat_3<K> &c, 
			  const K& k )
{
	return internal::intersection( c, s, k ) ;
}

} // end namespace internal


template <class K>
inline
Object
intersection(const Splat_3<K> &c, const Segment_3<K> &s)
{
  // return typename K::Intersect_3()(c,s);
  return internal::intersection(c,s,K());
}

template <class K>
inline
Object
intersection(const Segment_3<K> &s, const Splat_3<K> &c)
{
  // return typename K::Intersect_3()(c,s);
  return internal::intersection(c,s,K());
}

} // end namespace CGAL

#endif // CGAL_SPLAT_3_SEGMENT_3_INTERSECTION_H