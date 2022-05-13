#ifndef FUZZY_CAPSULE_3
#define FUZZY_CAPSULE_3

#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Triangle_3.h>
#include "SquaredDistance/squared_distance_segment_triangle_3.h"
#include "SquaredDistance/squared_distance_segment_triangle_vertices_3.h"


namespace CGAL {

template < class K >
class Fuzzy_capsule_3 {
	
public:

	typedef typename K::Point_3		Point_3 ;
	typedef typename K::Segment_3	Segment_3 ;
	typedef typename K::Triangle_3	Triangle_3 ;
	typedef typename K::Segment_3	Line_3 ;
	typedef typename K::Point_3		Point_d ; // The d-dimensional point, required by the concept
	typedef typename K::FT			FT ;	  // The float type, required by the concept



	// Default constructor
    Fuzzy_capsule_3(){ }

	// Constructor
	Fuzzy_capsule_3( const Segment_3& segment, FT radius, FT epsilon=FT(0) ) : m_segment(segment), m_radius(radius), m_eps(epsilon) 
	{ 			
		// avoid problems if eps > r
		if (m_eps > m_radius ) m_eps = m_radius ; 
		
	}

	// Test whether the capsule contains p
	bool contains( const typename K::Point_d& p ) const {
		
		// Compute the squared distance to the segment, it has to be less than the radius
		FT sqDist = CGAL::squared_distance( p, m_segment ) ;
		FT sqRadius = m_radius * m_radius ;
				
		return sqDist < sqRadius ;

	}

	// Test whether the inner approximation of the capsule intersects a rectangle associated with a node of a tree
	bool inner_range_intersects( const Kd_tree_rectangle<FT>& rectangle ) const {                      
		// std::cout << "    inner_range_intersects start" << std::endl ;

		if ( rectangle.dimension() != 3 )
			return false ; // Capsule is of dimension 3, so rectangle must be the same!

		// Build 12 triangles on the 3D box defined by rectangle (two triangles each face)
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

		// Check if any of the triangles is at a distance < (r-eps)*(r-eps) from the capsule's segment
		FT rSq = ( m_radius - m_eps ) * ( m_radius - m_eps ) ;

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
		// std::cout << "INNER rSq = " << rSq << std::endl ; 
		// std::cout << "INNER rectangle = " << p000 << p111 << std::endl ; 
		for ( it = triangles.begin(); it != triangles.end(); ++it ) {

			FT dist = 99999999.0 ;
			if ( !it->is_degenerate() ) {
				dist = CGAL::squared_distance( m_segment, *it ) ;			
				// std::cout << "INNER Distance = " << dist << std::endl ; 
			}
			else {
				dist = CGAL::squared_distance_segment_triangle_vertices_3(	m_segment, 
																			it->vertex(0), 
																			it->vertex(1), 
																			it->vertex(2) ) ;
				// std::cout << "INNER Degenerate! Distance = " << dist << std::endl ; 
			}

			if ( dist < rSq ) {
				// std::cout << "INNER RANGE INTERSECT" << std::endl ; 
				return true ;
			}

		}

		// If all triangles are within a distance > (r-eps)*(r-eps), then check if the segment endpoints lay within the cube defined by the Kd_tree_rectangle

		bool insidePlanes =  triangles[0].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE && 
							 triangles[2].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[4].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[6].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[8].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[10].supporting_plane().oriented_side( m_segment.point( 0 ) ) == CGAL::ON_POSITIVE_SIDE  &&
							 triangles[0].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[2].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[4].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[6].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[8].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  && 
							 triangles[10].supporting_plane().oriented_side( m_segment.point( 1 ) ) == CGAL::ON_POSITIVE_SIDE  ;

		return ( insidePlanes ) ;

		
		
		/*Triangle_3 tri_000_001_010( p000, p001, p010 ) ;
		if ( !tri_000_001_010.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_000_001_010 ) < rSq )
				return true	;	
		}
		Triangle_3 tri_010_011_001( p010, p011, p001 ) ;
		if ( !tri_010_011_001.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_010_011_001 ) < rSq )
				return true ;
		}
		Triangle_3 tri_010_110_011( p010, p110, p011 ) ;
		if ( !tri_010_110_011.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_010_110_011 ) < rSq )
				return true ;
		}
		Triangle_3 tri_110_111_011( p110, p111, p011 ) ;
		if ( !tri_110_111_011.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_110_111_011 ) < rSq )
				return true ;
		}
		Triangle_3 tri_110_100_111( p110, p100, p111 ) ;
		if ( !tri_110_100_111.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_110_100_111 ) < rSq )
				return true ;
		}
		Triangle_3 tri_100_101_111( p100, p101, p111 ) ;
		if ( !tri_100_101_111.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_100_101_111 ) < rSq )
				return true ;
		}
		Triangle_3 tri_101_100_000( p101, p100, p000 ) ;
		if ( !tri_101_100_000.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_101_100_000 ) < rSq )
				return true ;
		}
		Triangle_3 tri_000_001_101( p000, p001, p101 ) ;
		if ( !tri_000_001_101.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_000_001_101 ) < rSq )
				return true ;
		}
		Triangle_3 tri_001_011_101( p001, p011, p101 ) ;
		if ( !tri_001_011_101.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_001_011_101 ) < rSq )
				return true ;
		}
		Triangle_3 tri_011_111_101( p011, p111, p101 ) ;
		if ( !tri_011_111_101.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_011_111_101 ) < rSq )
				return true ;
		}
		Triangle_3 tri_010_100_110( p010, p100, p110 ) ;
		if ( !tri_010_100_110.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_010_100_110 ) < rSq )
				return true ;
		}
		Triangle_3 tri_000_100_010( p000, p100, p010 ) ;
		if ( !tri_000_100_010.is_degenerate() ) {
			if ( CGAL::squared_distance( m_segment, tri_000_100_010 ) < rSq )
				return true ;
		}*/

		// std::cout << "    inner_range_intersects ends" << std::endl ;

		// return false ;
	}

	bool outer_range_contains(const Kd_tree_rectangle<FT>& rectangle) const { 

		// std::cout << "    outer_range_contains start" << std::endl ;

		if ( rectangle.dimension() != 3 )
			return false ; // Capsule is of dimension 3, so rectangle must be the same!

		// Build 12 triangles on the 3D box defined by rectangle (two triangles each face)
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

		FT rSq = ( m_radius + m_eps ) * ( m_radius + m_eps ) ;

		// All points are inside the radius
		if ( CGAL::squared_distance( m_segment, p000 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p001 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p010 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p011 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p100 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p101 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p110 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
		if ( CGAL::squared_distance( m_segment, p111 ) >= rSq ) {
			// std::cout << "OUTER Does not contain" << std::endl ;
			return false ;
		}
					
		// std::cout << "OUTER Contains" << std::endl ;
		return true ;

	}

private:
	Segment_3 m_segment ;
	Line_3 m_segLine ;
    FT m_radius ;
    FT m_eps;
	
} ;

} // End namespace CGAL

#endif // FUZZY_CAPSULE_3