#ifndef POINTSET_SURFACE_MESHER_H
#define POINTSET_SURFACE_MESHER_H

#define DEBUG_OUTPUT
// #define PS_SURFACE_MESHER_ALL_PTS_TEST

// CGAL includes
#include <CGAL/iterator.h> 
#include <CGAL/Object.h> 
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/bounding_box.h>

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// Other related projects
#include "RobustStatistics/Ransac.h"
#include "SegmentQueryIntersectionOracle.h"

// Std includes
#include <fstream> // Debug only...
#include <iostream> // Debug only...
// #include <direct.h> // Debug only...


namespace CGAL {
	
template < class Kernel >	
class PointSetSurfaceMesherOracle {

	typedef	PointSetSurfaceMesherOracle								Self ;
	
public: 	
	typedef typename Kernel::FT										FT ;
    typedef typename Kernel::Ray_3									Ray_3 ;
    typedef typename Kernel::Line_3									Line_3 ;
	typedef typename Kernel::Segment_3								Segment_3 ;	
    typedef typename Kernel::Point_3								Point_3 ;
	typedef typename Kernel::Vector_3								Vector_3 ;    
	typedef typename Kernel::Plane_3								Plane_3 ;    	
	typedef typename Kernel::Sphere_3								Sphere_3 ;    	
	
	typedef Self													Surface_3 ;
	typedef CGAL::Object											Object ;
	typedef Point_3													Intersection_point ; 

	typedef SegmentQueryIntersectionOracle<K>				SQIO ;
	typedef boost::shared_ptr<SQIO>									SQIOPtr ;

private: 
	SQIOPtr															m_pSegQueryIntOracle ;
	FT																m_boundingSphereRadius ;	
	Point_3															m_boundingSphereCenter ;	
	
public:

	// Constructor
	PointSetSurfaceMesherOracle( SegmentQueryIntersectionOracle<K>* sqio ) {

		m_pSegQueryIntOracle = SQIOPtr( sqio ) ; 
	
		// Compute the squared radius of the minimum Sphere enclosing all input points
		Sphere_3 bSphere = m_pSegQueryIntOracle->computeBoundingSphere() ;
		
		// Store bounding sphere's data
		m_boundingSphereCenter = bSphere.center() ;
		m_boundingSphereRadius = CGAL::sqrt( bSphere.squared_radius() ) ; // Bounding sphere radius, limiting working space
		m_boundingSphereRadius += m_boundingSphereRadius*0.1 ; // Enlarge it a little bit...
		
	}

	// Getters/Setters
	SQIOPtr getSegmentQueryIntersectionOracle() const { return m_pSegQueryIntOracle ; } 
	
	// Other functions (deprecated...)
	bool isInDomain( const Point_3& query ) const {
		return CGAL::squared_distance( m_boundingSphereCenter, query ) < ( m_boundingSphereRadius * m_boundingSphereRadius ) ;
	}
	


	/* (Start) Implementation of concept for SurfaceMeshTraits_3::Construct_initial_points */
	
	class Construct_initial_points;	
    friend class Construct_initial_points;
	
	class Construct_initial_points
	{
		const Self& self;
	public:
		Construct_initial_points(const Self& self) : self(self)
		{
		}
		
		template <typename OutputIteratorPoints>
		OutputIteratorPoints operator() (   const Surface_3& surface,
											OutputIteratorPoints out,
											int n = 20 ) const
		{
#ifdef PS_SURFACE_MESHER_ALL_PTS_TEST
			/* Return a set of the input points as the initial points */
			std::vector<Point_3> allPts = self.getSegmentQueryIntersectionOracle()->getPts() ;
			
			typedef std::vector< Point_3 >::size_type size_type ;
			// int nb_initial_points = (std::min)( n, (int)allPts.size() ) ; // Check if there are enough points in the input set...
			int nb_initial_points = (int)allPts.size() ; // ALL the points...
			
			// Chose a bunch of them randomly
			for( size_type i = 0; i < nb_initial_points; i++ ) {
				*out++ = allPts[ i ] ;		
			}

#else
			/* Return a set of the input points as the initial points */
			//std::vector<Point_3> allPts = self.getSegmentQueryIntersectionOracle()->getPts() ;
			//
			//typedef std::vector< Point_3 >::size_type size_type ;
			// int nb_initial_points = (std::min)( n, (int)allPts.size() ) ; // Check if there are enough points in the input set...
			//
			//// Chose a bunch of them randomly
			//for( size_type i = 0; i < nb_initial_points; i++ ) {
			//	const int pos = CGAL::default_random.get_int( 0, static_cast< int >( allPts.size() ) ) ;
			//	*out++ = allPts[ pos ] ;		
			//}


			/* Check the intersection of some random rays on the points' bounding box and use them as input points*/

			// Get the input points' bounding box
			std::vector<Point_3> allPts = self.getSegmentQueryIntersectionOracle()->getPts() ;
			CGAL::Iso_cuboid_3<K> bBox = CGAL::bounding_box( allPts.begin(), allPts.end() ) ;
			
			// Throw random rays and compute intersections
			int i = 0 ;
			while ( i < n ) {

				// Compute a random segment
				Point_3 rp1( CGAL::default_random.get_double( bBox.xmin(), bBox.xmax() ), 
							 CGAL::default_random.get_double( bBox.ymin(), bBox.ymax() ), 
							 CGAL::default_random.get_double( bBox.zmin(), bBox.zmax() ) ) ;
				
				Point_3 rp2( CGAL::default_random.get_double( bBox.xmin(), bBox.xmax() ), 
							 CGAL::default_random.get_double( bBox.ymin(), bBox.ymax() ), 
							 CGAL::default_random.get_double( bBox.zmin(), bBox.zmax() ) ) ;

				Segment_3 segment( rp1, rp2 ) ;

				Object obj = surface.getSegmentQueryIntersectionOracle()->intersection( segment ) ;
				if ( const Point_3 *ip = CGAL::object_cast< Point_3 >( &obj ) ) {
					// Add the intersection point as initial point
					*out++ = *ip ;
					// Update number of valid intersections
					i++ ;
				}

			}
#endif

			return out ;

		}
	} ;
	
    Construct_initial_points construct_initial_points_object() const
    {
		return Construct_initial_points(*this);
    }

	/* (End)   Implementation of concept for SurfaceMeshTraits_3::Construct_initial_points */



	/* (Start) Implementation of concept for SurfaceMeshTraits_3::Intersect_3 */

	class Intersect_3;
    friend class Intersect_3;
	
	class Intersect_3 { 

		const Self& self ;

	public:

		Intersect_3( const Self& self ) : self( self ) {}

		// Intersection between the surface and a segment
		Object operator()(const Surface_3& surface, const Segment_3& segment) const {
			
			return surface.getSegmentQueryIntersectionOracle()->intersection( segment ) ;

		}

		// Intersection between the surface and a ray
		Object operator()(const Surface_3& surface, const Ray_3& ray) const	{			
			// Change the query to reuse the machinery already coded for segments
			Segment_3 querySegment = Segment_3( ray.source(), 
												ray.source()+( ray.to_vector()*999999 ) ) ; // Create a huge segment following the query ray
			return (*this)( surface, querySegment ) ;
		}
		
		// Intersection between the surface and a line
		Object operator()(const Surface_3& surface, const Line_3& line) const
		{
			// Change the query to reuse the machinery already coded for segments
			Point_3 pointOnLine = line.point() ;
			Segment_3 querySegment = Segment_3( pointOnLine + ( line.to_vector()*-999999 ), 
												pointOnLine+( line.to_vector()*999999 ) ) ; // Create a huge segment following the query line
			return (*this)( surface, querySegment ) ;
		}
		
	} ;

	Intersect_3 intersect_3_object() const { 
		return Intersect_3( *this ) ; 
	}	
	
	/* (End)   Implementation of concept for SurfaceMeshTraits_3::Intersect_3 */


} ; // End class PointSetSurfaceMesherOracle

} // End namespace CGAL

#endif // POINTSET_SURFACE_MESHER_H

