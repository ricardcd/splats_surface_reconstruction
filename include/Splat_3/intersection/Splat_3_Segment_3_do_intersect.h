// TODO: This class is very naive, it should not use the intersection query!

#ifndef CGAL_SPLAT_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_SPLAT_3_SEGMENT_3_DO_INTERSECT_H

#include <iostream>

#include <CGAL/kernel_basic.h>
#include "Splat_3/intersection/Splat_3_Segment_3_intersect.h"

namespace CGAL {

namespace internal {

template <class K>
bool
do_intersect( const Splat_3<K> &c, const Segment_3<K> &s, const K & k )
{
	return !intersection( c, s ).is_empty() ;
}

template <class K>
bool
do_intersect( const Segment_3<K> &s, const Splat_3<K> &c, const K & k )
{
	return !intersection( c, s ).is_empty() ;
}

} // end namespace internal



template <class K>
inline bool do_intersect(const Segment_3<K> &s, 
						 const Splat_3<K> &c)
{
	//return typename K::Do_intersect_3()(c,s);
	return internal::do_intersect( c, s, K() ) ;
}

template <class K>
inline bool do_intersect(const Splat_3<K> &c,
						 const Segment_3<K> &s)
{
	// return typename K::Do_intersect_3()(c,s);
	return internal::do_intersect( c, s, K() );
}

} // end namespace CGAL

#endif // CGAL_SPLAT_3_SEGMENT_3_DO_INTERSECT_H