#ifndef CGAL_AABB_TREE_SPLAT_3_MESHER_ORACLE_H
#define CGAL_AABB_TREE_SPLAT_3_MESHER_ORACLE_H

// CGAL includes
#include <CGAL/iterator.h> 
#include <CGAL/centroid.h>
#include <CGAL/AABB_tree.h>
#include "Splat_3/AABBTree/AABB_Splat_3_intersections.h"
#include "Splat_3/AABBTree/AABB_Splat_3_traits.h"
#include "Splat_3/AABBTree/AABB_Splat_3_primitive.h"
#include "Splat_3/Splat_3.h"
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Object.h>

// Boost includes
#include <boost/shared_ptr.hpp>

// Other Includes
#include <utility> 
#include <math.h>

// Other related projects
#include "RANSAC/ransac.h"

// Global variables (debug only)
int MesherIntersectionCounter = 0 ;
int MesherIntersectionCounterSegments = 0 ;
int MesherIntersectionCounterRays = 0 ;
int MesherIntersectionCounterLines = 0 ;
int foundMesherIntersectionCounter = 0 ;
int MesherSuccesfullyOptimizedPts = 0 ;
int numSegmentsMesherWithOneIntersection = 0 ;

// Definition of a null point... (maybe there is a better way to do that?)
#define NULL_POINT Point_3( 99999999, 99999999, 99999999 )

namespace CGAL {
	namespace Mesh_3 {
		namespace details {
			// -----------------------------------
			// Surface_patch_index_generator
			// To use patch_id enclosed in AABB_primitives or not
			// -----------------------------------
			template < typename Subdomain_index, typename Tag >
			struct Surface_patch_index_generator
			{
			  typedef std::pair<Subdomain_index,Subdomain_index>  Surface_patch_index;
			  typedef Surface_patch_index                         type;
  
			  template < typename Primitive_id >
			  Surface_patch_index operator()(const Primitive_id&)
			  { return Surface_patch_index(0,1); }
			};
  
			template < typename Subdomain_index >
			struct Surface_patch_index_generator<Subdomain_index, CGAL::Tag_true>
			{
			  typedef Subdomain_index       Surface_patch_index;
			  typedef Surface_patch_index   type;

			  template < typename Primitive_id >
			  Surface_patch_index operator()(const Primitive_id& primitive_id)
			  { 
				std::cout << "Surface_patch_index ( primitive_id ) = " << primitive_id->patch_id() << std::endl ;
				std::cout << "Patch id = " << primitive_id->patch_id() << std::endl ;
				return primitive_id->patch_id(); 
				}
			};


			// -----------------------------------
			// Index_generator
			// Don't use boost::variant if types are the same type
			// -----------------------------------
			template < typename Subdomain_index, typename Surface_patch_index >
			struct Index_generator
			{
			  typedef boost::variant<Subdomain_index,Surface_patch_index> Index;
			  typedef Index                                         type;
			};

			template < typename T >
			struct Index_generator<T, T>
			{
			  typedef T       Index;
			  typedef Index   type;
			};
			
		} // End namespace details
	} // End namespace Mesh_3
} // End namespace CGAL


namespace CGAL {
	
template < class Kernel, 
		   class SplatsIterator, 
		   class Use_patch_id_tag=Tag_false >	
class AABB_Tree_Splat_3_Mesher_Oracle {
		
	typedef	AABB_Tree_Splat_3_Mesher_Oracle Self ;
	
public: 	
	typedef typename Kernel::FT										FT ;
    typedef typename Kernel::Ray_3									Ray_3 ;
    typedef typename Kernel::Line_3									Line_3 ;
    typedef typename Kernel::Point_3								Point_3 ;
	typedef typename Kernel::Vector_3								Vector_3 ;
    typedef typename Kernel::Segment_3								Segment_3 ;	
	typedef typename CGAL::Splat_3<Kernel>							Splat_3 ;
	typedef typename Kernel::Plane_3								Plane_3 ;	
	
	typedef int														Subdomain_index ; // Type of indexes for cells of the input complex
	typedef boost::optional<Subdomain_index>						Subdomain ;
	typedef typename Mesh_3::details::
			Surface_patch_index_generator< Subdomain_index, 
										   Use_patch_id_tag >::type Surface_patch_index ; // Type of indexes for surface patch of the input complex
	typedef boost::optional<Surface_patch_index>		            Surface_patch;
	typedef typename Mesh_3::details::
			Index_generator< Subdomain_index, 
							 Surface_patch_index >::type			Index;    // Type of indexes to characterize the lowest dimensional face of the input complex on which a vertex lie
	
	// typedef Point_3													Intersection_point ; 
    typedef Self													Surface_mesher_traits_3 ;
	typedef Self													Surface_3 ;
	
	typedef CGAL::cpp0x::tuple< Point_3, Index, int >				Intersection ;

	// AABB tree
	typedef CGAL::AABB_Splat_3_primitive< Kernel, SplatsIterator >	Primitive ;
	typedef CGAL::AABB_Splat_3_traits< Kernel, Primitive >					AABB_Splat_3_traits ;
	typedef CGAL::AABB_tree< AABB_Splat_3_traits >					Tree ;
	typedef typename Tree::Object_and_primitive_id					Object_and_primitive_id;
	typedef typename Tree::Primitive_id								AABB_primitive_id ;
	typedef typename AABB_Splat_3_traits::Bounding_box				Bounding_box ;
	typedef boost::shared_ptr<Tree>									Tree_shared_ptr ; 	

	// Kernel_traits compatibility
	typedef Kernel R;

private: 
    Tree_shared_ptr													m_pTree ; // Pointer to the AABB tree
	std::vector< Point_3 >											m_initialPoints ;	
	double															m_distanceSigma  ;	
	double															m_ransacDistThres ;
	double															m_smallestRansacSegment ;	
	Intersection													m_computedIntersection ; // Due to the requirement of knowing beforehand if there exist an intersection in Do_intersect_surface, 
																							 // we store the intersection computation in this variable and return it in Compute_intersection
				
public: 

	// Default constructor
	AABB_Tree_Splat_3_Mesher_Oracle() { 
		m_pTree = boost::shared_ptr< Tree >() ; 
		m_distanceSigma  = 0 ;
		m_ransacDistThres = 0.1 ;
		m_smallestRansacSegment = 0.0 ;
	}
	
	// Constructor with the AABB tree as a parameter
	AABB_Tree_Splat_3_Mesher_Oracle( Tree *pTree, std::vector< Point_3 > initialPoints ) { 
		m_pTree = Tree_shared_ptr( pTree ) ;		
		m_initialPoints = initialPoints ;
		m_distanceSigma  = 0 ;
		m_ransacDistThres = 0.1 ;
		m_smallestRansacSegment = 0.0 ;
	}

	// Constructor with the AABB tree and the intersection type as parameter
	AABB_Tree_Splat_3_Mesher_Oracle(	Tree *pTree, 									
										std::vector< Point_3 > initialPoints,
										double distanceSigma = 0.5, 
										double ransacDistThres = 0.1,
										double smallestRansacSegment = 0.0 ) { 
		m_pTree = Tree_shared_ptr( pTree ) ;
		m_initialPoints = initialPoints ;		
		m_distanceSigma  = distanceSigma ;
		m_ransacDistThres = ransacDistThres ;
		m_smallestRansacSegment = smallestRansacSegment ;
	}

	void setParameters( Tree *pTree, 
						std::vector< Point_3 > initialPoints,
						double distanceSigma = 0.5, 
						double ransacDistThres = 0.1,
						double smallestRansacSegment = 0.0 ) { 
		m_pTree = Tree_shared_ptr( pTree ) ; 
		m_initialPoints = initialPoints ;
		m_distanceSigma  = distanceSigma ;
		m_ransacDistThres = ransacDistThres ;
		m_smallestRansacSegment = smallestRansacSegment ;
	}

	// Get a pointer to the AABBTree
	Tree* tree() const{ return m_pTree.get() ; }

	// Get the sigma (variance) describing the weight for a distance to the center
	double distanceSigma() const{ return m_distanceSigma ; }
	
	// Get the ransac distance threshold
	double ransacDistThres() const { return m_ransacDistThres ; }
		
	double smallestRansacSegment() const { return m_smallestRansacSegment ; }


	class Construct_intersection ;
    friend class Construct_intersection ;
	
	// Implementation of Construct_intersection class in another file to make this one less crowded...
	class Construct_intersection { 
				
		typedef boost::optional<typename Tree::Object_and_primitive_id> AABB_intersection;
		const Self& self ;
		
	public: 
		
		typedef typename Tree::Object_and_primitive_id	Obj_prim_id ;
		
		Construct_intersection( const Self& self ) : self( self ) {}

		// Intersection between the surface and a segment
		Intersection operator()( const Segment_3& segment ) const 
		{
			// Return the already computed intersection
			// (It should have been computed in Do_intersect_surface)
			// std::cout << "Recovering the computed intersection..." << std::endl ;
			return self.m_computedIntersection ;
			// return computeIntersection( segment ) ;
		}

		Intersection computeIntersection( const Segment_3& segment ) const {
			MesherIntersectionCounterSegments++ ;
			MesherIntersectionCounter++ ;

			// Check if there is any intersection
			if ( self.tree()->do_intersect( segment ) ) {
				// std::cout << "Some intersection" << std::endl ;
				foundMesherIntersectionCounter++ ;
				
				// Get the list of intersections
				std::vector< Obj_prim_id > obj_prim_vect ;
				std::vector< Point_3 >	pointsVector ;
				std::vector< Splat_3 >	splatsVector ;
				self.tree()->all_intersections( segment, std::back_inserter( obj_prim_vect ) ) ;
				typename std::vector< Obj_prim_id >::iterator it ;

				if ( obj_prim_vect.size() == 1 )
					numSegmentsMesherWithOneIntersection++ ;

				// Get the intersection points and its corresponding Splat_3 primitives
				for ( it = obj_prim_vect.begin(); it!=obj_prim_vect.end(); ++it ) {
					const Point_3* p = CGAL::object_cast< Point_3 >( &( it->first ) ) ;
									
					if ( p ) {
						pointsVector.push_back( *p ) ;
						Splat_3 s( *(it->second) ) ;
						splatsVector.push_back( s ) ;
					}
				}
				
				// Compute the intersection point				
				/*Object intersect = IntersectionWeightedCentroid( segment, 
																 pointsVector, 
																 splatsVector,
																 self.distanceSigma() ) ;*/
				Object intersect = IntersectionRansacWeightedCentroid(	segment, 
																		pointsVector, 
																		splatsVector,
																		self.distanceSigma(),
																		self.ransacDistThres() ) ;

				if ( const Point_3 *intPoint = CGAL::object_cast< Point_3 >(&intersect) ) {
					// std::cout << "Intersection : " << intPoint << std::endl ;
					return Intersection( *intPoint, 
										 self.index_from_surface_patch_index( self.make_surface_index( AABB_primitive_id() ) ),
										 2 ) ;
				}
				else {
					// std::cout << "Intersection outside the query segment" << std::endl ;
					return Intersection() ;
				}

			}
			else {
				// std::cout << "NO Intersection" << std::endl ;
				return Intersection() ;
			}
		}


		Object IntersectionWeightedCentroid( const Segment_3& segment,
											 const std::vector< Point_3 >& points, 
											 const std::vector< Splat_3 >& splats, 
											 const double& sigma ) const {
			typename std::vector< Point_3 >::const_iterator itP ;
			typename std::vector< Splat_3 >::const_iterator itC ;
			FT weightedDispSum = 0 ;
			FT totalWeight = 0 ;

			Point_3 ps = segment.source() ;
			Point_3 pe = segment.target() ;

			Vector_3 dir = segment.to_vector() ;
			// Unitize
			dir = dir / CGAL::sqrt( dir.squared_length() ) ;

			itC = splats.begin() ;
			for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
				// Displacement of the point along the segment
				FT disp = CGAL::sqrt( CGAL::squared_distance( ps, *itP ) ) ;

				// Compute the distance from the intersection to the center of the Circle_3
				FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

				// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
				FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
				FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;
				//FT weight = weighting( distCent, curSigma ) ;
				
				// Weight the current displacement
				weightedDispSum += disp * weight ;

				// Increase the total weight
				totalWeight += weight ;
			}

			FT meanDisp = weightedDispSum / totalWeight ;

			//if ( CGAL::sqrt( CGAL::squared_distance( ps, meanPoint ) ) < CGAL::sqrt( segment.squared_length() ) ) {
			if ( meanDisp < CGAL::sqrt( segment.squared_length() ) ) {
				Point_3 meanPoint = ps + ( meanDisp * dir ) ;
				return CGAL::make_object( meanPoint ) ;
			}
			else {							
				return Object() ;
			}
		}

		// The intersection is a centroid computed using RANSAC
		Object IntersectionRansacWeightedCentroid(  const Segment_3& segment,
													const std::vector< Point_3 >& points, 
													const std::vector< Splat_3 >& splats, 
													const double& sigma,
													const double& ransacDistThres ) const {

			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				// Compute weights
				typename std::vector< Point_3 >::const_iterator itP ;
				typename std::vector< Splat_3 >::const_iterator itC ;
				std::vector< double > weights ;
				itC = splats.begin() ;				
				for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
					// Compute the distance from the intersection to the center of the Circle_3
					FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

					// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
					//FT weight = centeredGaussian1D( 0, distCent, sigma ) ;
					//FT weight = centeredGaussian1D( distCent, 0, sigma ) ;
					FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
					FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;	
					//FT weight = weighting( distCent, curSigma ) ;	

					//if ( itC->score() > 0 )				
					//	weight = weight * itC->score() ; // Score also with weight

					//std::cout << "[IntersectionRansacWeightedCentroid] Radius = " << CGAL::sqrt( itC->squared_radius() ) << "   distCent = " << distCent << "   Sigma = " << curSigma << "   Weight = " << weight << std::endl ;

					weights.push_back( weight ) ;
				}
				
				bool exec =	ransac::ransacFitWeightedCentroid< Kernel >( points,
																		 weights,
																		 ransacDistThres,
																		 centroidRansac,
																		 inliers ) ;

				if ( !exec ) {
					//std::cerr << MesherIntersectionCounter << " Ransac problem..." << std::endl ;
					//return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
					return Object() ;
				}

				FT xSum = 0, ySum = 0, zSum = 0 ;
				std::vector< int >::iterator inIt ;
				std::vector< Point_3 > inlierPoints ;
				std::vector< Splat_3 > inlierCircles ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					inlierPoints.push_back( points[ *inIt ] ) ;					
					inlierCircles.push_back( splats[ *inIt ] ) ;
				}

				return IntersectionWeightedCentroid( segment, inlierPoints, inlierCircles, sigma ) ;
				// Point_3 intPoint = IntersectionWeightedCentroid( segment, inlierPoints, inlierCircles, sigma ) ;
				// return CGAL::make_object( intPoint ) ;

				// return centroidRansac ;

			}
			else {
				//std::cerr << MesherIntersectionCounter << " Not enough points (" << points.size() << ")" << std::endl ;
				if (  points.size() == 2 && 
					  CGAL::squared_distance( points[0], points[1] ) < ransacDistThres ) {
					return IntersectionWeightedCentroid( segment, points, splats, sigma ) ;
					// Point_3 intPoint = IntersectionWeightedCentroid( segment, points, splats, sigma ) ;
					// return CGAL::make_object( intPoint ) ;
				}
				else	
					//return points[0];
					return Object() ;
					//return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
			}
		}

		// Intersection between the surface and a ray
		Intersection operator()(const Ray_3& ray) const
		{
			MesherIntersectionCounterRays++ ;
			MesherIntersectionCounter++ ;
			//surface.increment_number_of_intersection_tests() ;
			
			// Change the query to reuse the machinery already coded for segments
			Segment_3 querySegment = Segment_3( ray.source(), 
												ray.source()+( ray.to_vector()*999999 ) ) ; // Create a huge segment following the query ray
			return (*this)( querySegment ) ;
		}
		
		// Intersection between the surface and a line
		Intersection operator()(const Line_3& line) const
		{
			MesherIntersectionCounterLines++ ;
			MesherIntersectionCounter++ ;
			////surface.increment_number_of_intersection_tests() ;
			
			// Change the query to reuse the machinery already coded for segments
			Point_3 pointOnLine = line.point() ;
			Segment_3 querySegment = Segment_3( pointOnLine + ( line.to_vector()*-999999 ), 
												pointOnLine+( line.to_vector()*999999 ) ) ; // Create a huge segment following the query line
			return (*this)( querySegment ) ;
		}


	} ;
	
	
	Construct_intersection construct_intersection_object() const { 
		return Construct_intersection( *this ) ; 
	} 

	

	class Construct_initial_points;	
    friend class Construct_initial_points;
	
	class Construct_initial_points
	{
		const Self& self;
	public:
		Construct_initial_points(const Self& self) : self(self)
		{
		}
		
		template <typename OutputIterator>
		OutputIterator operator() ( OutputIterator out ) const
		{
			// std::cout << "AABB_polyhedral_oracle: empty initial point set" << std::endl;		  
			for ( int i = 0; i<self.m_initialPoints.size(); i++ ) {
				*out++ = std::make_pair( self.m_initialPoints[ i ], Index(0) ) ;
			}
			return out ;
		}

	} ;
	
    Construct_initial_points construct_initial_points_object() const
    {
		return Construct_initial_points(*this);
    }
	
	
    template <class P>
    bool is_in_volume(const Surface_3& surface, const P& p)
    {
		const Bounding_box bbox = surface.tree()->root_bbox();
		if(p.x() < bbox.xmin() || p.x() > bbox.xmax())
			return false;
		if(p.y() < bbox.ymin() || p.y() > bbox.ymax())
			return false;
		if(p.z() < bbox.zmin() || p.z() > bbox.zmax())
			return false;
		
		const double diameter = max_bbox_length(bbox) * 2;
		
		typename CGAL::Random_points_on_sphere_3<Point_3> random_point(FT(1));
		typename Kernel::Construct_vector_3 vector =
        Kernel().construct_vector_3_object();
		typename Kernel::Construct_segment_3 segment =
        Kernel().construct_segment_3_object();
		typename Kernel::Construct_translated_point_3 translate =
        Kernel().construct_translated_point_3_object();
		typename Kernel::Construct_scaled_vector_3 scale =
        Kernel().construct_scaled_vector_3_object();
		
		const Segment_3 seg =
        segment(p, translate(p,
                             scale(vector(ORIGIN,
                                          *random_point),
                                   diameter)));
		return (surface.tree()->number_of_intersections(seg) % 2) == 1;
    }
		
	/* Just in Mesh_3 */

	struct Is_in_domain {
		Is_in_domain(const AABB_Tree_Splat_3_Mesher_Oracle& domain) : r_domain_(domain) {}

		Subdomain operator()(const Point_3& p) const {
			return Subdomain( Subdomain_index( 1 ) ) ; // IS THIS CORRECT???
		}
		private:
		const AABB_Tree_Splat_3_Mesher_Oracle& r_domain_;
	};

	Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

	struct Do_intersect_surface {
		Do_intersect_surface(const AABB_Tree_Splat_3_Mesher_Oracle& domain) : r_domain_(domain) {}

		Surface_patch operator()(const Segment_3& s) const {
			// std::cout << "HERE!!!!!!!" << std::endl ;

			// Check first the bounding_tree			
			boost::optional<AABB_primitive_id> primitive_id = r_domain_.m_pTree == 0 ? boost::none : r_domain_.m_pTree->any_intersected_primitive(s) ;
      
			if ( primitive_id ) { 				
				// Some planar primitives are intersected, compute if there really is an intersection
				const_cast<AABB_Tree_Splat_3_Mesher_Oracle&>( r_domain_ ).m_computedIntersection = r_domain_.construct_intersection_object().computeIntersection( s ) ;
				if ( r_domain_.m_computedIntersection != Intersection() ) {
				//Intersection inter = r_domain_.construct_intersection_object().computeIntersection( s ) ;
				//if ( inter != Intersection() ) {				
					// std::cout << "INTERSECTION" << std::endl ;
					return Surface_patch( r_domain_.make_surface_index( *primitive_id ) ) ;				
				}
				else {
					// std::cout << "NO INTERSECTION" << std::endl ;
					return Surface_patch() ;
				}
			}
			else {				
				// std::cout << "NO INTERSECTED PRIMITIVE" << std::endl ;
				return Surface_patch() ;
			}
		}

		Surface_patch operator()(const Ray_3& r) const {
			// Change the query to reuse the machinery already coded for segments
			Segment_3 querySegment = Segment_3( r.source(), 
												r.source()+( r.to_vector()*999999 ) ) ; // Create a huge segment following the query ray
			return (*this)( querySegment ) ;
		}

		Surface_patch operator()(const Line_3& l) const {
			// Change the query to reuse the machinery already coded for segments
			Point_3 pointOnLine = l.point() ;
			Segment_3 querySegment = Segment_3( pointOnLine + ( l.to_vector()*-999999 ), 
												pointOnLine+( l.to_vector()*999999 ) ) ; // Create a huge segment following the query line
			return (*this)( querySegment ) ;
		}

		private:
		const AABB_Tree_Splat_3_Mesher_Oracle& r_domain_;
	};

	Do_intersect_surface do_intersect_surface_object() const
	{
		return Do_intersect_surface(*this);
	}

	/**
	* Returns the index to be stored in a vertex lying on the surface identified
	* by \c index.
	*/
	Index index_from_surface_patch_index(const Surface_patch_index& index) const { 		
		return Index(index); 		
	}

	/**
	* Returns the index to be stored in a vertex lying in the subdomain
	* identified by \c index.
	*/
	Index index_from_subdomain_index(const Subdomain_index& index) const
	{ return Index(index); }

	/**
	* Returns the \c Surface_patch_index of the surface patch
	* where lies a vertex with dimension 2 and index \c index.
	*/
	Surface_patch_index surface_patch_index(const Index& index) const
	{ return boost::get<Surface_patch_index>(index); }

	/**
	* Returns the index of the subdomain containing a vertex
	*  with dimension 3 and index \c index.
	*/
	Subdomain_index subdomain_index(const Index& index) const
	{ return boost::get<Subdomain_index>(index); }

	/* Definition of private functions */
	private:
    double max_bbox_length(const Bounding_box& bbox) const
    {
		return (std::max)(bbox.xmax()-bbox.xmin(),
						  (std::max)(bbox.ymax()-bbox.ymin(),
									 bbox.zmax()-bbox.zmin()));
    }

public:
	Surface_patch_index make_surface_index(
    const AABB_primitive_id& primitive_id = AABB_primitive_id() ) const
  {
	Mesh_3::details::Surface_patch_index_generator<Subdomain_index,Use_patch_id_tag>
      generator;
    
    return generator(primitive_id);
  }



private:
	/* Auxiliar functions regarding gaussians */
	
	// Function to compute a gaussian centered at center in one dimension
	static FT centeredGaussian1D( const FT& x, const FT& center = 0, const FT& sigma = 0.5 ) {
		FT a = 1 / ( sigma * sqrt( 2.0 * CGAL_PI ) ) ;

		// --- Debug (Start) ---
		//std::cout << "x = " << x << std::endl ;
		//std::cout << "center = " << center << std::endl ;
		//std::cout << "sigma = " << sigma << std::endl ;
		//std::cout << "a = " << a << std::endl ;
		//std::cout << "exp part = " << exp( - ( ( ( x - center ) * ( x - center ) ) / ( 2.0*sigma*sigma) ) ) << std::endl ;
		//
		//std::cout << "( ( x - center ) * ( x - center ) ) = " << ( ( x - center ) * ( x - center ) ) << std::endl ;
		//std::cout << "( 2.0*sigma*sigma) = " << ( 2.0*sigma*sigma) << std::endl ;
		//std::cout << "Gaussian = " <<  a * exp( - ( ( ( x - center ) * ( x - center ) ) / ( 2.0*sigma*sigma) ) ) << std::endl ;
		// --- Debug  (End)  ---

		// General formula
		return a * exp( - ( ( ( x - center ) * ( x - center ) ) / ( 2.0*sigma*sigma) ) ) ;
	}

	// Another weighting function
	static inline FT weighting( FT x, FT radius ) {
		//spline patched Gaussian
		if ( radius < x )
			return 0 ;
		else {
			x /= radius ;
			if( x < 0.5f )
				return exp( -8.0*x*x ) ;
			else {
				x = 1.0 - x ;
				x = x*x ;
				return 2.1653645317858030703*x*x ;
			}
		}
	}
	
	/*static FT gaussianKernel3D( Vector_3 x, FT sigma ) {		
		return ( 1 / pow( CGAL::sqrt( 2 * CGAL_PI ) * sigma, 3 ) ) * exp( - x.squared_length() / ( 2*sigma*sigma ) ) ;
	}*/


  }; // end class AABB_Tree_Splats_Mesher_Oracle

} // End namespace CGAL

#endif // CGAL_AABB_TREE_SPLAT_3_MESHER_ORACLE_H


