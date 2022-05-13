#ifndef SQUARED_DISTANCE_SEGMENT_TRIANGLE_VERTICES_3_H
#define SQUARED_DISTANCE_SEGMENT_TRIANGLE_VERTICES_3_H

namespace CGAL {

namespace internal {

// Convenience function... should be revised
template <class K>
inline
typename K::FT
squared_distance_segment_triangle_vertices_3(	const typename K::Segment_3& origin,
												const typename K::Point_3& p1,
												const typename K::Point_3& p2,
												const typename K::Point_3& p3,
												const K& k )
{
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();

  const FT dist_origin_p1 = sq_distance(origin,p1);
  const FT dist_origin_p2 = sq_distance(origin,p2);
  const FT dist_origin_p3 = sq_distance(origin,p3);

  if (   dist_origin_p2 >= dist_origin_p1
      && dist_origin_p3 >= dist_origin_p1 )
  {
    return dist_origin_p1;
  }
  if ( dist_origin_p3 >= dist_origin_p2 )
  {
    return dist_origin_p2;
  }

  return dist_origin_p3;
}

} // End namespace internal

template <class K>
inline
typename K::FT
squared_distance_segment_triangle_vertices_3( const Segment_3<K>& origin,
												const Point_3<K>& p1,
												const Point_3<K>& p2,
												const Point_3<K>& p3 ) {
	return internal::squared_distance_segment_triangle_vertices_3( origin, p1, p2, p3, K() ) ;
}



} // End namespace CGAL

#endif // SQUARED_DISTANCE_SEGMENT_TRIANGLE_VERTICES_3_H