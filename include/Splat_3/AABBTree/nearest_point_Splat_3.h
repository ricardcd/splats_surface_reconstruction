#ifndef NEAREST_POINT_SPLAT_3_H_
#define NEAREST_POINT_SPLAT_3_H_

// Code based on the file "CGAL/internal/AABB_tree/nearest_point_triangle_3.h"

#include <CGAL/kernel_basic.h>
#include <CGAL/enum.h>


namespace CGAL {


namespace internal {

/**
 * @brief returns true if p is inside circle c. If p is not inside c,
 * result is the nearest point of c from p. WARNING: it is assumed that
 * c and p are on the same plane.
 * @param p the reference point
 * @param c the circle
 * @param result if p is not inside c, the nearest point of c from p
 * @param k the kernel
 * @return true if p is inside 
 */
template <class K>
inline
bool
is_inside_circle_3(const typename K::Point_3& p,
                     const Splat_3<K>& splat,
                     typename K::Point_3& result,
                     const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector =
    k.construct_vector_3_object();
  typename K::Compare_squared_distance_3 compare_sq_distance =
    k.compare_squared_distance_3_object();
 
  // Since they are in the same plane, the inside test boils down to a distance comparison with the radius
  bool inside = ( compare_sq_distance( splat.center(), p, splat.squared_radius() ) == CGAL::SMALLER ) ;

  if ( inside ) {
	  return true ;
  }
  else {
	  // Outside, get an approximation of the closest point on the local surface
	  // WARNING: Note that it is JUST AN APPROXIMATION

	  // Construct the vector that goes from the Circle center to the point
	  Vector_3 v = vector( splat.center(), p ) ;

	  Point_3 sphereIntersection = Point_3( splat.center() + ( v * CGAL::sqrt( splat.squared_radius() ) ) ) ;

	  // Project the closest intersection with the sphere to the underlying local surface (we consider this a good approximation...)
	  result = splat.project( sphereIntersection ) ;

	  return false ;
  }
}



/**
 * @brief Computes the closest_point from origin between bound and
 * any point of circle.
 * @param origin the origin point
 * @param circle the circle
 * @param bound the farthest point
 * @param k the kernel
 * @return nearest point: bound or a point inside circle
 */
template <class K>
typename K::Point_3
nearest_point_3(const typename K::Point_3& origin,
                const Splat_3<K>& splat,
                const typename K::Point_3& bound,
                const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();
  typename K::Compare_squared_distance_3 compare_sq_distance =
    k.compare_squared_distance_3_object();
  typename K::Construct_supporting_plane_3 supporting_plane =
    k.construct_supporting_plane_3_object();
  typename K::Construct_projected_point_3 projection =
    k.construct_projected_point_3_object();

  // Distance from origin to bound
  const FT bound_sq_dist = sq_distance(origin, bound);

  // Project origin on triangle supporting plane
  const Point_3 proj = splat.project( origin ) ;

  // If point is projected outside, return bound
  if ( compare_sq_distance(origin, proj, bound_sq_dist) == CGAL::LARGER )
  {
    return bound;
  }

  // Check if it is inside the circle
  Point_3 moved_point;
  bool inside = is_inside_circle_3( proj, splat, moved_point, k ) ;
  
  // If proj is inside circle, return it
  if ( inside )
  {
    return proj;
  }

  // if it is closest to origin than bound
  if ( compare_sq_distance(origin, moved_point, bound_sq_dist)
                                                        == CGAL::LARGER )
  {
    return bound;
  }

  return moved_point;
}


}  // end namespace internal


template <class K>
inline
Point_3<K>
nearest_point_3(const Point_3<K>& origin,
                const Splat_3<K>& splat,
                const Point_3<K>& bound)
{
  return internal::nearest_point_3(origin, splat, bound, K());
}

}  // end namespace CGAL


#endif // NEAREST_POINT_SPLAT_3_H_