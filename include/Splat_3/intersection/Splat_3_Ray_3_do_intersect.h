// TODO: This class is very naive, it should not use the intersection query!

#ifndef CGAL_SPLAT_3_RAY_3_DO_INTERSECT_H
#define CGAL_SPLAT_3_RAY_3_DO_INTERSECT_H

#include <iostream>

#include <CGAL/kernel_basic.h>
#include "Splat_3/intersection/Splat_3_Ray_3_intersect.h"



namespace CGAL {

namespace internal {

template <class K>
bool
do_intersect( const Splat_3<K> &c, const Ray_3<K> &r, const K & k )
{	
	return !intersection( c, r ).is_empty() ;
}

template <class K>
bool
do_intersect( const Ray_3<K> &r, const Splat_3<K> &c, const K & k )
{
	return !intersection( c, r ).is_empty() ;
}

} // end namespace internal



template <class K>
inline bool do_intersect(const Ray_3<K>  &r, 
						 const Splat_3<K> &c)
{
	// return typename K::Do_intersect_3()(c,r);
	return internal::do_intersect(c,r, K());
}

template <class K>
inline bool do_intersect(const Splat_3<K> &c,
						 const Ray_3<K>  &r)
{
	// return typename K::Do_intersect_3()(c,r);
	return internal::do_intersect(c,r, K());
}

} // end namespace CGAL

#endif // CGAL_SPLAT_3_RAY_3_DO_INTERSECT_H