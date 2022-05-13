#ifndef CGAL_SPLAT_3_LINE_3_INTERSECTION_H
#define CGAL_SPLAT_3_LINE_3_INTERSECTION_H

#include <iostream>

#include <CGAL/kernel_basic.h>
#include <CGAL/intersections.h>
#include "Splat_3/Splat_3.h"

namespace CGAL {

namespace internal {

template <class K>
Object
intersection( const Splat_3<K> &s, 
			  const typename K::Line_3 &l,
			  const K & k )
{
	typedef typename K::Point_3			Point_3 ;
	typedef typename K::Line_3			Line_3 ;
		
	// Compute the intersection between the plane and the line
	Point_3 aPoint = l.point(0) ;
	// Move the start/end point far away (simulate infinity of the line by making a big )
	Point_3 startPoint = aPoint + ( -999999 * l.to_vector() ) ;
	// Point_3 endPoint = aPoint + ( 999999 * l.to_vector() ) ;
	Point_3 result ;

	double residuals[1] ;
	bool do_intersect = s.monge().intersection( startPoint, l.to_vector(), result, residuals ) ;

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
intersection( const typename K::Line_3 &l, 
			  const Splat_3<K> &c,
			  const K & k )
{
	return internal::intersection( c, l ) ;
}

} // end namespace internal

/*template <class K>
inline
Object
intersection( const typename K::Line_3 &l, 
			  const Splat_3<K> &c,
			  const K & k )
{
	return internal::intersection( c, l ) ;
}*/

template <class K>
inline
Object
intersection(const Splat_3<K> &c, const Line_3<K> &l)
{
  // return typename K::Intersect_3()(c,l);
  return internal::intersection(c,l,K());
}

template <class K>
inline
Object
intersection(const Line_3<K> &l, const Splat_3<K> &c)
{
  // return typename K::Intersect_3()(c,l);
  return internal::intersection(c,l,K());
}

} // end namespace CGAL

#endif // CGAL_SPLAT_3_LINE_3_INTERSECTION_H