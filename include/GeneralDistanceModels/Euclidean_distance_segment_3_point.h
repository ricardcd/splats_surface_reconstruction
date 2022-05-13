#include <CGAL/squared_distance_3.h>
#include <limits>
#include <CGAL/Point_3.h>
#include "SquaredDistance/squared_distance_segment_triangle_3.h"
#include "SquaredDistance/squared_distance_segment_triangle_vertices_3.h"

template < class K >
class Euclidean_distance_segment_3_point {

public :
  typedef typename K::Point_3		Point_3 ;
  typedef typename K::Segment_3		Segment_3 ;
  typedef typename K::Triangle_3	Triangle_3 ;
  typedef typename K::Segment_3		Query_item;
  typedef typename K::FT			FT;

  double transformed_distance( const Segment_3& s, const Point_3& p ) const {
    return CGAL::squared_distance( s, p ) ;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const Segment_3& s,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& rectangle ) const {
    
	FT minX = rectangle.min_coord( 0 ) ;
	FT minY = rectangle.min_coord( 1 ) ;
	FT minZ = rectangle.min_coord( 2 ) ;
	FT maxX = rectangle.max_coord( 0 ) ;
	FT maxY = rectangle.max_coord( 1 ) ;
	FT maxZ = rectangle.max_coord( 2 ) ;

	Point_3 p000( minX, minY, minZ ) ;
	Point_3 p001( minX, minY, maxZ ) ;
	Point_3 p010( minX, maxY, minZ ) ;
	Point_3 p011( minX, maxY, maxZ ) ;
	Point_3 p100( maxX, minY, minZ ) ;
	Point_3 p101( maxX, minY, maxZ ) ;
	Point_3 p110( maxX, maxY, minZ ) ;
	Point_3 p111( maxX, maxY, maxZ ) ;

	std::vector < Triangle_3 > triangles ;
	triangles.push_back( Triangle_3( p000, p010, p001 ) ) ; // [0]
	triangles.push_back( Triangle_3( p010, p011, p001 ) ) ;
	triangles.push_back( Triangle_3( p010, p110, p011 ) ) ; // [2]
	triangles.push_back( Triangle_3( p110, p111, p011 ) ) ;
	triangles.push_back( Triangle_3( p110, p100, p111 ) ) ; // [4]
	triangles.push_back( Triangle_3( p100, p101, p111 ) ) ; 
	triangles.push_back( Triangle_3( p101, p100, p000 ) ) ; // [6]
	triangles.push_back( Triangle_3( p000, p001, p101 ) ) ;
	triangles.push_back( Triangle_3( p001, p011, p101 ) ) ; // [8]
	triangles.push_back( Triangle_3( p011, p111, p101 ) ) ;
	triangles.push_back( Triangle_3( p010, p100, p110 ) ) ; // [10]
	triangles.push_back( Triangle_3( p000, p100, p010 ) ) ;
		
	typename std::vector < Triangle_3 >::iterator it ;
	FT minSqDist = std::numeric_limits<double>::infinity() ;
	for ( it = triangles.begin(); it != triangles.end(); ++it ) {

		FT curSqDist = std::numeric_limits<double>::infinity() ;
		if ( !it->is_degenerate() ) {
			curSqDist = CGAL::squared_distance( s, *it ) ;			
		}
		else {
			curSqDist = CGAL::squared_distance_segment_triangle_vertices_3(	s, 
																			it->vertex(0), 
																			it->vertex(1), 
																			it->vertex(2) ) ;
		}

		if ( curSqDist < minSqDist ) {			
			minSqDist = curSqDist ;
		}

	}

    return minSqDist ;
  }


  // The function max_distance_to_rectangle not implemented yet! 
  // just Nearest Neighbor searches work for the moment...
  template <class TreeTraits>
  double max_distance_to_rectangle(const Segment_3& s,
                                   const CGAL::Kd_tree_rectangle<TreeTraits>& rectangle ) const {
		
    return 0.0 ;

  }

  double new_distance(double& dist, double old_off, double new_off,
                      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }

}; // end of struct Distance