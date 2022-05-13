// TODO: This class is very naive, it should not use the intersection query!

#ifndef CGAL_SPLAT_3_LINE_3_DO_INTERSECT_H
#define CGAL_SPLAT_3_LINE_3_DO_INTERSECT_H

#include <iostream>

#include <CGAL/kernel_basic.h>
#include "Splat_3/intersection/Splat_3_Line_3_do_intersect.h"

namespace CGAL {

namespace internal {

template <class K>
bool
do_intersect( const Splat_3<K> &c, const Line_3<K> &l, const K & k )
{
	return !intersection( c, l ).is_empty() ;
}

template <class K>
bool
do_intersect( const Line_3<K> &l, const Splat_3<K> &c, const K & k )
{
	return !intersection( c, l ).is_empty() ;
}

} // end namespace internal



template <class K>
inline bool do_intersect(const Line_3<K>  &l, 
						 const Splat_3<K> &c)
{
	// return typename K::Do_intersect_3()(c,l);
	return internal::do_intersect( c, l, K() ) ;
}

template <class K>
inline bool do_intersect(const Splat_3<K> &c,
						 const Line_3<K>  &l)
{
	// return typename K::Do_intersect_3()(c,l);
	return internal::do_intersect( c, l, K() ) ;
}


} // end namespace CGAL

#endif // CGAL_SPLAT_3_LINE_3_DO_INTERSECT_H