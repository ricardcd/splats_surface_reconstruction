#ifndef SPLATSUNSIGNEDDISTANCEFUNCTIONNCUT_H
#define SPLATSUNSIGNEDDISTANCEFUNCTIONNCUT_H

/* Includes */
// CGAL
#include <CGAL/trace.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
	// Spatial searching
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include "SplatsCreation/KDTreeSearch/My_Search_traits_3.h"
#include "SplatsCreation/KDTreeSearch/Fuzzy_capsule_3.h"
	// Triangulation (refinement)
#include "SplatsIsosurfaceMesher/Splats_distance_function_triangulation_3.h"
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Robust_circumcenter_filtered_traits_3.h>
	// AABB_tree
#include <CGAL/AABB_tree.h>
#include "Splat_3/AABBTree/AABB_Splat_3_primitive.h"
#include "Splat_3/AABBTree/AABB_Splat_3_intersections.h"
#include "Splat_3/AABBTree/AABB_Splat_3_traits.h"
	// Random number generator
 #include <CGAL/Random.h>
	// Other
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/bounding_box.h>
// Related to this project
#include "Splat_3/Splat_3.h"
#include "RobustStatistics/Ransac.h"
#include "RobustStatistics/Models/FitCentroid/CentroidModelEstimator.h"
#include "UnsignedDistanceFunctionRefiner/Unsigned_distance_refine_triangulation.h"
#include "Octree/CGAL_Octree.h"
// Boost
#include <boost/shared_ptr.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/filesystem.hpp>
// Eigen
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <Eigen/Eigenvalues> 
#include <Eigen/Sparse> 
// Std
#include <algorithm>
#include <cassert>
// NCut
#include "SplatsIsosurfaceMesher/NCut.h"
// Undefined value define
#ifndef UNDEFINED_VALUE
	#define UNDEFINED_VALUE 999999.0
#endif
// Test
// #include "SplatsIsosurfaceMesher/GenericFunction.h"

//namespace boost {
//  enum vertex_component_t { vertex_component = 111 } ;
//  enum vertex_triangulation_handle_t { vertex_triangulation_handle = 112 } ;
//  BOOST_INSTALL_PROPERTY( vertex, component ) ;
//  BOOST_INSTALL_PROPERTY( vertex, triangulation_handle ) ;
//}
//
//using namespace boost ;

template < class K >
class SplatsDistanceFunctionNCut { // : public GenericFunction<K> {

public:
	/* Typedefs */
	typedef	SplatsDistanceFunctionNCut		Self ;

	typedef typename K::FT												FT ;
	typedef typename K::Point_3											Point_3 ;
	typedef typename K::Sphere_3										Sphere_3 ;
	typedef typename K::Ray_3											Ray_3 ;
	typedef typename K::Vector_3										Vector_3 ;
	typedef typename CGAL::Splat_3< K >									Splat_3 ;
	typedef typename K::Segment_3										Segment_3 ;
	typedef typename K::Plane_3											Plane_3 ;
		
	typedef typename K::Iso_cuboid_3									BBox ; // Needed for marching cubes...

	typedef typename CGAL::Monge_via_jet_fitting<K, K>					Monge_via_jet_fitting ;
	typedef typename Monge_via_jet_fitting::Monge_form					Monge_form ;

	/// Internal 3D triangulation, of type Reconstruction_triangulation_3.
	// Note: unsigned_distance_refine_triangulation() requires a robust circumcenter computation.
	typedef CGAL::Splats_distance_function_triangulation_3<
				CGAL::Robust_circumcenter_filtered_traits_3<K> >		Triangulation;

	// Repeat Triangulation types
	typedef typename Triangulation::Triangulation_data_structure		Triangulation_data_structure;
	typedef typename K::Triangle_3										Triangle ;
	typedef typename K::Tetrahedron_3									Tetrahedron ;
	typedef typename Triangulation::Cell_handle							Cell_handle ;
	typedef typename Triangulation::Vertex_handle						Vertex_handle ;
	typedef typename Triangulation::Cell								Cell ;
	typedef typename Triangulation::Vertex								Vertex ;
	typedef typename Triangulation::Facet								Facet ;
	typedef typename Triangulation::Edge								Edge ;
	typedef typename Triangulation::Cell_circulator						Cell_circulator ;
	typedef typename Triangulation::Facet_circulator					Facet_circulator ;
	typedef typename Triangulation::Cell_iterator						Cell_iterator ;
	typedef typename Triangulation::Facet_iterator						Facet_iterator ;
	typedef typename Triangulation::Edge_iterator						Edge_iterator ;
	typedef typename Triangulation::Vertex_iterator						Vertex_iterator ;
	typedef typename Triangulation::Point_iterator						Point_iterator ;
	typedef typename Triangulation::Finite_vertices_iterator			Finite_vertices_iterator ;
	typedef typename Triangulation::Finite_cells_iterator				Finite_cells_iterator ;
	typedef typename Triangulation::Finite_facets_iterator				Finite_facets_iterator ;
	typedef typename Triangulation::Finite_edges_iterator				Finite_edges_iterator ;
	typedef typename Triangulation::All_cells_iterator					All_cells_iterator ;
	typedef typename Triangulation::Locate_type							Locate_type ;

	// AABB Tree related definitions
	typedef typename std::vector< Splat_3 >::const_iterator						SplatsIterator ;
	typedef typename CGAL::AABB_Splat_3_primitive< K, SplatsIterator >			Primitive ;
	typedef typename CGAL::AABB_Splat_3_traits< K, Primitive >							AABB_Splat_3_traits ;
	typedef typename CGAL::AABB_tree< AABB_Splat_3_traits >						SplatsAABBTree ;
	typedef typename SplatsAABBTree::Object_and_primitive_id					Object_and_primitive_id;
	typedef typename boost::shared_ptr< SplatsAABBTree >						pSplatsAABBTree ;

    // Octree related
    typedef typename CGAL_Octree::CGAL_Octree< K >                      Octree ;
    typedef typename CGAL_Octree::CubicBounds< K >                      CubicBounds ;

	// Eigen related
	typedef Eigen::Matrix< FT, 
						   Eigen::Dynamic, 
						   1 >						EVector ;
	typedef Eigen::Matrix< FT,
						   ::Eigen::Dynamic,
						   ::Eigen::Dynamic>		EMatrix ;
	typedef Eigen::SparseMatrix< FT >				SparseMatrix ;

	// --- Nearest neighbors search definitions ---

	// Property maps to store point (center) and its corresponding splat
	typedef boost::tuple< Point_3, Splat_3 >							Point_and_splat;

	class Point_and_splat_property_map{
		typedef Point_3 value_type ;
		typedef const value_type& reference ;
		typedef const Point_and_splat& key_type ;
		typedef boost::readable_property_map_tag category;

		// Get function for the property map
		friend reference get( Point_and_splat_property_map mppm,
							  typename Point_and_splat_property_map::key_type p ) {
			return boost::get<0>(p);
		}

	} ;
		
	typedef typename CGAL::Search_traits_3< K >							TreeTraitsBase ;
	typedef CGAL::Search_traits_adapter< Point_and_splat,
										 Point_and_splat_property_map,
										 TreeTraitsBase >				TreeTraits ;
	typedef typename CGAL::Kd_tree<TreeTraits>							NNTree ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< TreeTraits >	Neighbor_search ;
	typedef boost::shared_ptr<NNTree>									pNNTree ;		
	typedef typename CGAL::Fuzzy_sphere< TreeTraits >					Sphere_range_query ;
	
	typedef boost::tuple< Point_3, Vertex_handle >						Point_and_vertex ;
	class Point_and_vertex_property_map{
		typedef Point_3 value_type ;
		typedef const value_type& reference ;
		typedef const Point_and_vertex& key_type ;
		typedef boost::readable_property_map_tag category;

		// Get function for the property map
		friend reference get( Point_and_vertex_property_map mppm,
							  typename Point_and_vertex_property_map::key_type p ) {
			return boost::get<0>(p);
		}

	} ;

	typedef typename CGAL::Search_traits_3< K >							PVTreeTraitsBase ;
	typedef CGAL::Search_traits_adapter< Point_and_vertex,
										 Point_and_vertex_property_map,
										 PVTreeTraitsBase >				PVTreeTraits ;
	typedef typename CGAL::Kd_tree<PVTreeTraits>						PVNNTree ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< TreeTraits >	PV_Neighbor_search ;
	typedef boost::shared_ptr<PVNNTree>									pPVNNTree ;		
	typedef typename CGAL::Fuzzy_sphere< PVTreeTraits >					PV_sphere_range_query ;

	typedef CGAL::Search_traits_3< K >										PtsTreeTraits ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< PtsTreeTraits >	PtsNeighborSearch ;
	typedef typename PtsNeighborSearch::Tree											PtsNNTree ;
	typedef boost::shared_ptr<PtsNNTree>									pPtsNNTree ;
	typedef typename CGAL::Fuzzy_sphere< PtsTreeTraits >					Pts_sphere_range_query ;

	// --- Boost graphs ---

	struct VertexProperties {
		int					component ;
		Vertex_handle		triVertHandle ;
		int					ind ;
	} ;

	template <typename ComponentMap>
	struct vertexComponent {

		vertexComponent() {}

		vertexComponent(ComponentMap component, int f_component) : m_component(component), m_f_component(f_component) {}

		template <typename Vertex>
		bool operator()(const Vertex& v) const {
			return (get(m_component, v) == m_f_component);
		}

		ComponentMap m_component;
		int m_f_component;
	} ;

	// Graph type (custom properties for vertices, float weight as properties for edges)
	typedef boost::adjacency_list< boost::vecS,
								   boost::vecS,
								   boost::undirectedS,
								   VertexProperties,
								   boost::property< boost::edge_weight_t, FT > > Graph ;

	typedef typename boost::property_map< Graph, int VertexProperties::* >::type ComponentMap ;
	typedef boost::filtered_graph< Graph,
								   boost::keep_all,
								   vertexComponent<ComponentMap> > FilteredGraph ;



	/* Functions */

	SplatsDistanceFunctionNCut( const std::vector< Splat_3 >& splats,
							const FT& gaussianRadiusFactor = 0.2, 
							const FT smoothPow = 4.0, 
							const FT areaConst = 1e-5,
							const FT functionErrorBound = 1e-3,
							const FT radiusEdgeBound = 2.5, 
							const int fixKnn = -1,
							const double bandRadius = 0.0,
							const double defMax = 1000,
							const bool contourAtMedian = false,
							const int bandKnn = -1, 
                            const double epsilon = 0.001,
                            const int odepth = -1,
                            const bool debugOutput = false)
							: m_tr(new Triangulation) // , GenericFunction<K>()
	{
		m_inputType = SPLATS ;
		m_smoothPow = smoothPow ;
		m_areaConst = areaConst ;
		m_correct = true ; // By default, we assume no error...
		m_radiusEdgeBound = radiusEdgeBound ;
		m_k = fixKnn ;
		m_defMax = 1000 ;
		m_functionErrorBound = functionErrorBound ;
		m_bandKnn = bandKnn ;
		m_timeImplicit = 0.0 ;
		m_timeSigning = 0.0 ;
		m_epsilon = epsilon ;
		m_debugOutput = debugOutput;

		// --- Debug (Start) ---
		if (m_debugOutput) {
            // Initialize output directory
            boost::filesystem::path dir( "./_OUT" ) ;
            if (boost::filesystem::create_directory( dir ) )
                std::cout << "[DEBUG] Created Debug Output Directory (./_OUT)" << "\n";
		}
		// --- Debug  (End)  ---

		// Extract centers
		typename std::vector< Splat_3 >::const_iterator it ;
		std::vector< Point_3 > points ;
		for( it = splats.begin(); it != splats.end(); ++it ) {
			points.push_back( it->center() ) ;
		}
		
		// Build initial triangulation
        if ( odepth > 0 ) {
            // Use an octree to decrease the number of base points in the implicit function
            std::cout << "  - Decimating points in an octree of depth " << odepth << "..." << std::flush ;
            Octree octree ;
            CubicBounds obounds( points ) ;
            octree.build( points, odepth, obounds ) ;
            std::vector<Point_3> meanPts ;
            octree.getMeanPoints( meanPts ) ;
            std::cout << "done" << std::endl ;

            std::cout << "  - Building triangulation (decimated " << points.size() << " to " << meanPts.size() << ")..." << std::flush ;
            m_tr->insert( meanPts.begin(), meanPts.end() ) ;
            std::cout << "done" << std::endl ;
        }
        else {
            // Use all points to build the implicit function
            std::cout << "  - Building triangulation with all the points..." << std::flush ;
            m_tr->insert( points.begin(), points.end() ) ;
            std::cout << "done" << std::endl ;
        }

        this->addBoundingBoxSteinerPoints(points);
		
		// Construct the KD-Tree with these center points
		m_pTree = pNNTree( new NNTree( boost::make_zip_iterator(boost::make_tuple( points.begin(), splats.begin() ) ),
									   boost::make_zip_iterator(boost::make_tuple( points.end(), splats.end() ) ) ) ) ;
		
		FT average_spacing = 0 ;
		if ( gaussianRadiusFactor < 1 ) {
			m_gaussianH = CGAL::sqrt( gaussianRadiusFactor * gaussianRadiusFactor * this->bounding_sphere().squared_radius() ) ;
			m_bandRadius = CGAL::sqrt( bandRadius * bandRadius * this->bounding_sphere().squared_radius() ) ;
		}
		else {
			// Computes average spacing
			average_spacing = CGAL::compute_average_spacing( points.begin(), points.end(), 6 /* knn = 1 ring */ ) ;
			// Compute gaussian weight H variable
			m_gaussianH = gaussianRadiusFactor * average_spacing ;
			m_bandRadius = bandRadius * average_spacing ;
		}
		
		// --- Debug (Start) ---
		std::cout << "  - Some Debug Info:" << std::endl ;
		std::cout << "    Function error bound = " << m_functionErrorBound << std::endl ;
		if ( gaussianRadiusFactor > 1 ) {
			std::cout << "    average_spacing = " << average_spacing << std::endl ;
			std::cout << "    m_gaussianH (w.r.t. AS) = " << m_gaussianH << std::endl ;
			std::cout << "    m_bandRadius (w.r.t. AS) = " << m_bandRadius << std::endl ;
		}
		else {
			std::cout << "    m_gaussianH (w.r.t. BSR) = " << m_gaussianH << std::endl ;
			std::cout << "    m_bandRadius (w.r.t. BSR) = " << m_bandRadius << std::endl ;
		}
		// --- Debug  (End)  ---

		CGAL::Timer timer ;

		/* Banded version (should work better)*/
		// Compute the (unsigned) implicit function
		timer.start() ;
		m_correct = compute_implicit_function_band() ;
		m_timeImplicit = timer.time() ;
		if ( !m_correct ) return ;

		// --> TEST!!!
		/*std::cout << "  - Unitizing" << std::endl ;
		unitize_f() ;*/

		// Sign the function
		timer.reset() ;
		m_correct = sign_implicit_function_band_ncut() ;
		// m_correct = sign_implicit_function_band_ncut_single_cc() ;
		m_timeSigning = timer.time() ;
		if ( !m_correct ) return ;
		// sign_implicit_function_extended_band_ncut() ;

		/* Global version */
		// Compute the (unsigned) implicit function
		// compute_implicit_function() ;
		// Sign the function
		// sign_implicit_function_ncut() ;
		
		// --- Debug (Start) ---
		//Finite_vertices_iterator v, e;
		//int i = 0 ;
		//std::cout << "Checking values at vertices:" << std::endl ;
		//for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
		//	 v != e;
		//	 ++v, i++)
		//{
		//	std::cout<< v->f() << std::endl ;
		//}
		// --- Debug  (End)  ---

		if ( contourAtMedian ) {
			shift_f( median_value_at_input_vertices() ) ;
		}

	}

	SplatsDistanceFunctionNCut( const std::vector< Point_3 >& points,
								const FT& gaussianRadiusFactor = 0.2, 
								const FT smoothPow = 4.0, 
								const FT areaConst = 1e-5,
								const FT functionErrorBound = 1e-3,
								const FT radiusEdgeBound = 2.5, 
								const int fixKnn = -1,
								const double bandRadius = 0.0,
								const double defMax = 1000, 
								const bool contourAtMedian = false,
								const int bandKnn = -1,
                                const double epsilon = 0.001,
                                const int odepth = -1,
                                const bool debugOutput = false)
							: m_tr(new Triangulation) // , GenericFunction<K>()
	{
		m_inputType = POINTS ;
		if ( fixKnn <= 0 ) {
			std::cout << "[ERROR] When using raw point sets, the k-NN must be greater than 0!" << std::endl ;
			m_correct = false ;
			return ;
		}
		m_smoothPow = smoothPow ;
		m_areaConst = areaConst ;
		m_correct = true ; // By default, we assume no error...
		m_radiusEdgeBound = radiusEdgeBound ;
		m_k = fixKnn ;
		m_defMax = defMax ;
		m_functionErrorBound = functionErrorBound ;		// Error bound for linear approximation
		m_bandKnn = bandKnn ;
		m_timeImplicit = 0.0 ;
		m_timeSigning = 0.0 ;
		m_epsilon = epsilon ;
		m_debugOutput = debugOutput;
		
		// --- Debug (Start) ---
		if (m_debugOutput) {
            // Initialize output directory
            boost::filesystem::path dir("./_OUT");
            if (boost::filesystem::create_directory(dir))
                std::cout << "[DEBUG] Created Debug Output Directory (./_OUT)" << "\n";
        }
		// --- Debug  (End)  ---

		// Build initial triangulation
        if ( odepth > 0 ) {
            // Use an octree to decrease the number of base points in the implicit function
            std::cout << "  - Decimating points in an octree of depth " << odepth << "..." << std::flush ;
            Octree octree ;
            CubicBounds obounds( points ) ;
            octree.build( points, odepth, obounds ) ;
            std::vector<Point_3> meanPts ;
            octree.getMeanPoints( meanPts ) ;
            std::cout << "done" << std::endl ;

            std::cout << "  - Building triangulation (decimated " << points.size() << " to " << meanPts.size() << ")..." << std::flush ;
            m_tr->insert( meanPts.begin(), meanPts.end() ) ;
            std::cout << "done" << std::endl ;
        }
        else {
            // Use all points to build the implicit function
            std::cout << "  - Building triangulation with all the points..." << std::flush ;
            m_tr->insert( points.begin(), points.end() ) ;
            std::cout << "done" << std::endl ;
        }

        this->addBoundingBoxSteinerPoints(points);
		
		// Construct the KD-Tree with these center points
		m_ptsTree = pPtsNNTree( new PtsNNTree( points.begin(), points.end() ) ) ;
		
		FT average_spacing = 0 ;
		if ( gaussianRadiusFactor < 1 ) {
			m_gaussianH = CGAL::sqrt( gaussianRadiusFactor * gaussianRadiusFactor * this->bounding_sphere().squared_radius() ) ;
			m_bandRadius = CGAL::sqrt( bandRadius * bandRadius * this->bounding_sphere().squared_radius() ) ;
		}
		else {
			// Computes average spacing
			average_spacing = CGAL::compute_average_spacing( points.begin(), points.end(), 6 /* knn = 1 ring */ ) ;
			// Compute gaussian weight H variable
			m_gaussianH = gaussianRadiusFactor * average_spacing ;
			m_bandRadius = bandRadius * average_spacing ;
		}
		
		// --- Debug (Start) ---
		std::cout << "  - Some Debug Info:" << std::endl ;
		std::cout << "    Function error bound = " << m_functionErrorBound << std::endl ;
		if ( gaussianRadiusFactor > 1 ) {
			std::cout << "    average_spacing = " << average_spacing << std::endl ;
			std::cout << "    m_gaussianH (w.r.t. AS) = " << m_gaussianH << std::endl ;
			std::cout << "    m_bandRadius (w.r.t. AS) = " << m_bandRadius << std::endl ;
		}
		else {
			std::cout << "    m_gaussianH (w.r.t. BSR) = " << m_gaussianH << std::endl ;
			std::cout << "    m_bandRadius (w.r.t. BSR) = " << m_bandRadius << std::endl ;
		}
		// --- Debug  (End)  ---
		
		CGAL::Timer timer ;

		/* Banded version (should work better)*/
		// Compute the (unsigned) implicit function
		timer.start() ;
		m_correct = compute_implicit_function_band() ;
		m_timeImplicit = timer.time() ;
		if ( !m_correct ) return ;

		// Shift the zero! (should improve results...)
		shift_f( -minimum_value_at_input_vertices() ) ;

		// Sign the function
		timer.reset() ;
		m_correct = sign_implicit_function_band_ncut() ;
		// m_correct = sign_implicit_function_band_ncut_single_cc() ;
		m_timeSigning = timer.time() ;
		if ( !m_correct ) return ;


		// sign_implicit_function_extended_band_ncut() ;

		/* Global version */
		// Compute the (unsigned) implicit function
		// compute_implicit_function() ;
		// Sign the function
		// sign_implicit_function_ncut() ;
		
		// --- Debug (Start) ---
		//Finite_vertices_iterator v, e;
		//int i = 0 ;
		//std::cout << "Checking values at vertices:" << std::endl ;
		//for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
		//	 v != e;
		//	 ++v, i++)
		//{
		//	std::cout<< v->f() << std::endl ;
		//}
		// --- Debug  (End)  ---

		if ( contourAtMedian ) {
			shift_f( median_value_at_input_vertices() ) ;
		}

	}

    void addBoundingBoxSteinerPoints(const std::vector< Point_3 >& points) {
        // Determine min/max of the given set of points
        FT min_x = points[0].x() ;
        FT min_y = points[0].y() ;
        FT min_z = points[0].z() ;

        FT max_x = points[0].x() ;
        FT max_y = points[0].y() ;
        FT max_z = points[0].z() ;


        typename std::vector<Point_3>::const_iterator it ;
        for ( it = points.begin();
              it != points.end();
              ++it )
        {
            if ( it->x() < min_x ) min_x = it->x() ;
            if ( it->y() < min_y ) min_y = it->y() ;
            if ( it->z() < min_z ) min_z = it->z() ;
            if ( it->x() > max_x ) max_x = it->x() ;
            if ( it->y() > max_y ) max_y = it->y() ;
            if ( it->z() > max_z ) max_z = it->z() ;
        }

        // Steiner points
        std::vector<Point_3> steinerBB;
        steinerBB.push_back(Point_3(min_x, min_y, min_z));
        steinerBB.push_back(Point_3(min_x, max_y, min_z));
        steinerBB.push_back(Point_3(min_x, max_y, max_z));
        steinerBB.push_back(Point_3(min_x, min_y, max_z));
        steinerBB.push_back(Point_3(max_x, min_y, min_z));
        steinerBB.push_back(Point_3(max_x, max_y, min_z));
        steinerBB.push_back(Point_3(max_x, max_y, max_z));
        steinerBB.push_back(Point_3(max_x, min_y, max_z));

        m_tr->insert( steinerBB.begin(), steinerBB.end() ) ;
    }

	FT getTimeImplicitFunctionCreation() const {
		return m_timeImplicit ;
	}

	FT getTimeSigningFunction() const {
		return m_timeSigning ;
	}

	/// Get a guess for the bounding sphere containing the input points
	Sphere_3 bounding_sphere() const {
		return m_tr->input_points_bounding_sphere();
	}
	

	/// Computes enlarged geometric bounding sphere of the embedded triangulation.
	Sphere_3 enlarged_bounding_sphere(FT ratio) const
	{
		Sphere_3 bsphere = this->bounding_sphere(); // triangulation's bounding sphere

		// --- Debug (Start) ---
		//std::cout << "BSphere:" << std::endl ;
		//std::cout << bsphere << std::endl ;
		//std::cout << "Enlarged BSphere:" << std::endl ;
		//std::cout << Sphere_3(bsphere.center(), bsphere.squared_radius() * ratio*ratio) << std::endl ;
		// --- Debug  (End)  ---

		return Sphere_3(bsphere.center(), bsphere.squared_radius() * ratio*ratio);
	}


	BBox bounding_box() const {
		std::vector< Point_3 > points ;
		Finite_vertices_iterator v, e;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			points.push_back( v->point() ) ;
		}

		return CGAL::bounding_box( points.begin(), points.end() ) ;
	}


	/// Compute the implicit function at the vertices of a refined triangulation (similar to an octree decomposition)
	bool compute_implicit_function_band() {

		// Initialize implicit function values at triangulation vertices...
		/*std::cout << "  - Initializing implicit function..."  ;
		Finite_vertices_iterator v, e;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			v->f() = eval( v->point() ) ;
		}
		std::cout << "done" << std::endl ;*/

		CGAL::Timer task_timer; task_timer.start() ;
		std::cout << "  - Delaunay refinement:" << std::endl ;

		// Delaunay refinement
		const FT radius_edge_ratio_bound = m_radiusEdgeBound ; // Modified --> Original = 2.5
		const unsigned int max_vertices = (unsigned int)1e7 ; // max 10M vertices
		const FT enlarge_ratio = 1.5 ;
		const FT radius = sqrt( bounding_sphere().squared_radius() ) ; // get triangulation's radius
		const FT cell_radius_bound = radius/5. ; // large
		std::cout << "    - Parameters:" << std::endl ;
		std::cout << "        radius_edge_ratio_bound = " << m_radiusEdgeBound << std::endl ;
		std::cout << "        enlarge_ratio = " << enlarge_ratio << std::endl ;
		std::cout << "        cell_radius_bound = " << cell_radius_bound << std::endl ;
		std::cout << "    - Refining..." ;

		unsigned int nb_vertices_added = delaunay_refinement( radius_edge_ratio_bound,
															  cell_radius_bound,
															  m_gaussianH,
															  max_vertices,
															  enlarge_ratio ) ;

		// Prints status
		std::cout << "done, added " << nb_vertices_added << " Steiner points in "
						  << task_timer.time() << " seconds, "
						  << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
						  << std::endl ;
		task_timer.reset() ;

		std::cout << "  - Evaluating implicit function..." ;
		int verbosePercent = 0 ; int verboseLastPercent = 0 ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;

		// Computes the Implicit function eval()
		// at each vertex of the triangulation.
		// Finite_vertices_iterator v, e;
		int i = 0 ;
		m_maxVal = 0 ;
		m_numVertOnBand = 0 ;
		int numVert = m_tr->number_of_vertices() ;
		Finite_vertices_iterator v, e;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v)
		{
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numVert ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}

			// Eval function
			if ( v->f() < 0 ) {
				if ( v->type() == Triangulation::INPUT && m_k > 0 ) m_k++ ; // Workaround for k-NN distance based when used on input points
				v->f() = eval( v->point() ) ; 
				if ( v->type() == Triangulation::INPUT && m_k > 0 ) m_k-- ; // Workaround for k-NN distance based when used on input points
			}
			
			// Assign also a linear index ONLY IF IT IS ON THE BAND
			// if ( v->f() != 999999 ) {
			if ( !this->invalidValue( v->f() ) ) {
				v->ind() = i++ ;
				m_numVertOnBand++ ;
			}
			else {
				if ( m_bandRadius > m_gaussianH ) {
					// Add the points on a secondary band with a default value
					Point_3 p = v->point() ;
					Neighbor_search search( *m_pTree, p, 1 ) ;
				
					/*Splat_3 curSplat = boost::get<1>( search.begin()->first ) ;
					Point_3 projectedPoint = curSplat.project( p ) ;
					FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;*/

					Splat_3 curSplat = boost::get<1>( search.begin()->first ) ;
					FT dist = CGAL::sqrt( CGAL::squared_distance( curSplat.center(), p ) ) ;

					// Inside second band
					if ( dist < m_bandRadius ) {
						// std::cout << "Inside second band" << std::endl ;
						v->ind() = i++ ;
						m_numVertOnBand++ ;

						if ( m_bandKnn > 0 ) { // TODO: FACTORIZE THIS PART!!!!!!!!!!!!!
							// Also replace the function with the approximation provided by kNN
							Neighbor_search search( *m_pTree, p, m_bandKnn ) ;
							bool any = false ;
							FT sumWeightedDist = 0.0 ;
							FT sumWeights = 0.0 ;
							for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
								// Distance measure
								Splat_3 curSplat = boost::get<1>( it->first ) ;
								Point_3 projectedPoint = curSplat.project( p ) ;
								FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;
								
								//if ( dist > m_bandRadius ) {
								//	std::cout << "Projection outside the band!!!!!" << std::endl ;
								//	dist = 0.0 ;
								//}

								// Weighting
								if ( dist < m_bandRadius ) {
									FT weight = this->weightGaussian( dist, m_bandRadius ) / m_bandKnn ;

									sumWeightedDist += weight * dist ;
									sumWeights += weight ;
									if ( !any )
										any = true ;
								}
							}
							if ( any ) 
								v->f() = sumWeightedDist / sumWeights ;
							else {
								v->ind() = -1 ;
								i-- ;
								m_numVertOnBand-- ;
							}

							/*if ( v->f() > m_bandRadius ) {
								std::cout << "[ERROR] f() = " << v->f() << " > " << m_bandRadius << std::endl ;
							}

							if ( v->f() != v->f() ) {
								std::cout << "[ERROR] f() = NAN = " << v->f() << " | sumWeightedDist = " << sumWeightedDist << " | sumWeights = " << sumWeights << search.size() << std::endl ;
							}*/

							// std::cout << "secondary band knn = " << v->f() << std::endl ;
						}
					}
					else {
						v->ind() = -1 ; 
					}
				}
				else {
					v->ind() = -1 ; 
				}
			}

			if ( v->ind() >= 0 && m_maxVal < v->f() && !this->invalidValue( v->f() ) ) {
				m_maxVal = v->f() ;
			}
		}
		
		// Prints status
		CGAL_TRACE_STREAM << "\b\b\b\b100%, done (" << task_timer.time() << " s)" << std::endl ;
		task_timer.reset() ;

		// --- Debug (Start) ---
		std::cout << "[DEBUG] Implicit function max val = " << m_maxVal << std::endl ;
		// Compute the value of the function at each of its vertices
		FT mv = median_value_at_input_vertices() ;
		std::cout << "[DEBUG] Median value at input points: " << mv << std::endl ;
		// --- Debug  (End)  ---

		return true;

	}

	bool sign_implicit_function_band_ncut_single_cc() {

		std::cout << "  - Sign the implicit function:" << std::endl ;
		
		// SparseMatrix W = buildLaplacian() ;

		int numVert = m_tr->number_of_vertices() ;
		int numEdges = m_tr->number_of_edges() ;

		std::cout << "    - Building graph..." ;
		int verbosePercent = 0 ; int verboseLastPercent = 0 ; int i = 0 ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;

		std::vector< Eigen::Triplet<double> > triplets ;
		int k = 0 ;
		Finite_edges_iterator it;
		CGAL::Timer task_timer; task_timer.start() ;
		
		for ( it = m_tr->finite_edges_begin();
			  it !=	m_tr->finite_edges_end();
			  ++it, i++ )
		{
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numEdges ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}

			Cell_handle c = (it->first) ;
			Vertex_handle v0 = c->vertex( it->second ) ;
			Vertex_handle v1 = c->vertex( it->third ) ;

			// std::cout << v0->ind() << ", " << v1->ind() << std::endl ;

			if ( v0->ind() >= 0 && v1->ind() >= 0 ) {
				// if ( v0->f() == 999999 || v1->f() == 999999 ) {
				if ( invalidValue( v0->f() ) || invalidValue( v1->f() ) ) {
					triplets.push_back( Eigen::Triplet<double>( v0->ind(), v1->ind(), m_defMax ) ) ;
					triplets.push_back( Eigen::Triplet<double>( v1->ind(), v0->ind(), m_defMax ) ) ;
				}
				else {
					FT val = ( v0->f() + v1->f() ) / 2 ;
					FT len = CGAL::sqrt( CGAL::squared_distance( v0->point(), v1->point() ) ) ;

					// w = edge_weight( v0->point(), v1->point() ) ;
					// w = pow( ( ( v0->f() + v1->f() ) / 2 ), m_smoothPow ) ;
					// w = this->weightInverseGaussian( ( v0->f() + v1->f() ) / 2, m_gaussianH ) ;
					// w = this->weightInverseWendland( val, m_gaussianH ) ;

					/*val = val*val ;
					len = len*len ;*/
					// FT wf = this->weightInverseWendland( val, m_maxVal ) ;
					FT wf = val ;
					// FT we = this->weightInverseWendland( len, maxEdgeLength ) ;
					wf = pow( wf, m_smoothPow ) ;
					FT w = wf ;
					// std::cout << "w = " << w << std::endl ;
					// w = wf*we ;
					//std::cout << "wf = " << wf << " / we = " << we << " / w = " << w << std::endl ;
					//std::cout << "val = " << val << " / len = " << len << std::endl ;

					triplets.push_back( Eigen::Triplet<double>( v0->ind(), v1->ind(), w ) ) ;
					triplets.push_back( Eigen::Triplet<double>( v1->ind(), v0->ind(), w ) ) ;
				}
			}

			// Full matrices --> untractable
			/*W( v0->ind(), v1->ind() ) = w ;
			W( v1->ind(), v0->ind() ) = w ;*/

			// Sparse matrices --> incremental filling is slow even if we reserve the space beforehand
			// W.insert( v0->ind(), v1->ind() ) = w ;
			// W.insert( v1->ind(), v0->ind() ) = w ;

			// Sparse matrices --> use triplets to initialize the matrix (fastest!)
			
		}
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << " s)" << std::endl ;

		// Compute connected components
		typedef boost::adjacency_list< boost::vecS, 
									   boost::vecS, 
									   boost::undirectedS,
									   boost::no_property, 
									   boost::property< boost::edge_weight_t, FT > > Graph ;
		Graph g ;
		std::vector< Eigen::Triplet<double> >::iterator itt ;
		for ( itt = triplets.begin(); itt != triplets.end(); ++itt ) {
			// std::cout << "( " << itt->row() << ", " << itt->col() << " ) = " << itt->value() << std::endl ;
			if ( itt->value() != itt->value() ) {
				std::cout << "NAN!!!!!!!!!!!!" << std::endl ;
				std::cout << "( " << itt->row() << ", " << itt->col() << " ) = " << itt->value() << std::endl ;
			}
			boost::add_edge( static_cast<int>( itt->row() ), 
							 static_cast<int>( itt->col() ), 
											   itt->value(), g ) ;
		}
		std::vector<int> comp( num_vertices(g) ) ;
		int numCC = boost::connected_components( g, &comp[0] ) ;
		std::cout << "    - Number of Connected Components = " << numCC << std::endl ;

		SparseMatrix W( m_numVertOnBand, m_numVertOnBand ) ;
		// W.reserve( triplets.size() + m_numVertOnBand ) ;
		W.setFromTriplets( triplets.begin(), triplets.end() ) ;

		EMatrix eigenVectors;
		EVector eigenValues;
		//std::tie( eigenvectors, eigenvalues ) = std::move( ncut( W, 2 ) ) ;
		
		std::cout << "    - Normalized cut..." << std::flush ; 
		task_timer.reset() ;
		ncut<FT>( W, 2, eigenVectors, eigenValues ) ;
		std::cout << "done (" << task_timer.time() << " s)" << std::endl ;
		std::cout << "[DEBUG] EigenValues = < " << eigenValues(0) << ", " << eigenValues(1) << " >" << std::endl ;
		
		// --- Debug (Start) ---
        std::ofstream ofio, ofoo, ofe, oft;
		if (m_debugOutput) {
            ofio.open("./_OUT/InsideVertices_OPT.xyz", std::ios_base::out);
            ofoo.open("./_OUT/OutsideVertices_OPT.xyz", std::ios_base::out);
            ofe.open("./_OUT/2ndSmallestEigenVector.txt", std::ios_base::out);
            oft.open("./_OUT/TriVert.xyz", std::ios_base::out);
        }
		// --- Debug  (End)  ---

		// Sign function according to second smallest eigenvector
		std::cout << "    - Signing function..." << std::flush ; 
		task_timer.reset() ;
		EVector secondEigenVector = eigenVectors.col(0) ;
		// std::cout << "[DEBUG] 2nd largest EigenVector = < " << secondEigenVector << " >" << std::endl ;
		Finite_vertices_iterator v, e;
		bool any = false ; 
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			oft << v->point() << std::endl ;
			if ( v->ind() >= 0 ) {
				ofe << secondEigenVector( v->ind() ) << std::endl ;
				if ( secondEigenVector( v->ind() ) < 0 ) {
					v->f() = -v->f() ;
					ofio << v->point() << std::endl ;
					if ( !any ) any = true ;
				}
				else {
					ofoo << v->point() << std::endl ;
				}
			}
		}
		ofio.close() ;
		ofoo.close() ;
		oft.close() ;
		ofe.close() ;
		std::cout << "done (" << task_timer.time() << ")" << std::endl ;

		return any ;

	}


	//struct edge_in_largest_component {
	//	edge_in_largest_component() { }
	//	edge_in_largest_component( std::vector<int>& components ) : m_weight(weight) { }
	//	
	//	template <typename Edge>
	//	bool operator()(const Edge& e) const {
	//		int s = boost::source( *edgePair.first, g ) ;
	//		int t = boost::target( *edgePair.first, g ) ;

	//		return ( m_components[s] == largestCCInd && m_components[t] == largestCCInd )
	//	}

	//	// Attributes
	//	std::vector<int> m_components ;
	//};

	


	bool sign_implicit_function_extended_band_ncut() {

		std::cout << "  - Sign the implicit function:" << std::endl ;

		int numVert = m_tr->number_of_vertices() ;
		int numEdges = m_tr->number_of_edges() ;

		Finite_vertices_iterator v, e;
		std::vector< Point_3 > points ;
		std::vector< Vertex_handle > vertices ;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			points.push_back( v->point() ) ;
			vertices.push_back( v ) ;
		}

		PVNNTree tree( boost::make_zip_iterator(boost::make_tuple( points.begin(), vertices.begin() ) ),
					   boost::make_zip_iterator(boost::make_tuple( points.end(), vertices.end() ) ) ) ; 

		std::cout << "    - Building graph..." ;
		int verbosePercent = 0 ; int verboseLastPercent = 0 ; int i = 0 ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;

		std::vector< Eigen::Triplet<double> > triplets ;
		int k = 0 ;
		CGAL::Timer task_timer; task_timer.start() ;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++ )
		{
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numVert ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}
			
			if ( v->ind() >= 0 ) {
				/* Radial NN */
				std::vector< Point_and_vertex > neighs ;
				PV_sphere_range_query searchSphere( v->point(), m_gaussianH/2 ) ;
				tree.search( std::back_inserter( neighs ), searchSphere ) ;

				// std::cout << "Neighs = " << neighs.size() << std::endl ;

				typename std::vector< Point_and_vertex >::iterator it ;
				for ( it = neighs.begin(); it!= neighs.end(); ++it ) {
					// std::cout << "here0" << std::endl ;
					Vertex_handle v1 = boost::get<1>( *it ) ;
					// std::cout << "here1" << std::endl ;
					if ( v1->ind() >= 0 && v->ind() != v1->ind() ) {
						// FT w = pow( ( ( v->f() + v1->f() ) / 2 ), m_smoothPow ) ;
						FT val = min_along_edge( v->point(), v1->point() ) ;
						// FT w = pow( val, m_smoothPow ) ;
						FT w = this->weightInverseGaussian( val, m_gaussianH ) ;
						triplets.push_back( Eigen::Triplet<double>( v->ind(), v1->ind(), w ) ) ;
						triplets.push_back( Eigen::Triplet<double>( v1->ind(), v->ind(), w ) ) ;
					}
				}
			}

		}
		SparseMatrix W( m_numVertOnBand, m_numVertOnBand ) ;
		W.reserve( triplets.size() + m_numVertOnBand ) ;
		W.setFromTriplets( triplets.begin(), triplets.end() ) ;
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << " s)" << std::endl ;

		std::cout << "    - Filling diagonal..." << std::flush ; task_timer.reset() ;
		for ( i = 0; i < m_numVertOnBand; i++ ) {
			FT val = 0 ;
			for( typename SparseMatrix::InnerIterator ii( W, i ); ii; ++ii) {
				val += fabs( ii.value() ) ;
			}
			W.insert( i, i ) = val ; 
		}
		std::cout << "done (" << task_timer.time() << " s)" << std::endl ;

		EMatrix eigenVectors;
		EVector eigenValues;
		//std::tie( eigenvectors, eigenvalues ) = std::move( ncut( W, 2 ) ) ;

		std::cout << "    - Normalized cut..." << std::flush ; task_timer.reset() ;
		ncut<FT>( W, 2, eigenVectors, eigenValues ) ;
		std::cout << "done (" << task_timer.time() << " s)" << std::endl ;

		std::cout << "[DEBUG] EigenValues = < " << eigenValues(0) << ", " << eigenValues(1) << " >" << std::endl ;
		
		// --- Debug (Start) ---
        std::ofstream ofio, ofoo, oft;
        if (m_debugOutput) {
            ofio.open("./_OUT/InsideVertices_OPT.xyz", std::ios_base::out);
            ofoo.open("./_OUT/OutsideVertices_OPT.xyz", std::ios_base::out);
            oft.open("./_OUT/TriVert.xyz", std::ios_base::out);
        }
		// --- Debug  (End)  ---

		// Sign function according to second smallest eigenvector
		std::cout << "    - Signing function..." << std::flush ; task_timer.reset() ;
		EVector secondEigenVector = eigenVectors.col(0) ;
		// std::cout << "[DEBUG] 2nd largest EigenVector = < " << secondEigenVector << " >" << std::endl ;
		
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			oft << v->point() << std::endl ;
			if ( v->ind() >= 0 ) {
				if ( secondEigenVector( v->ind() ) < 0 ) {
					v->f() = -v->f() ;
					ofio << v->point() << std::endl ;
				}
				else {
					ofoo << v->point() << std::endl ;
				}
			}
		}
		ofio.close() ;
		ofoo.close() ;
		oft.close() ;
		std::cout << "done (" << task_timer.time() << ")" << std::endl ;
		
		return true ;

	}


	// We need to do things like this because of the inheritance...
	virtual FT operator()(const Point_3& p) const {
		return f(p) ;
	}

	/// 'ImplicitFunction' interface: evaluates the implicit function at a given 3D query point.
	virtual FT f(const Point_3& p) const {

		m_hint = m_tr->locate(p,m_hint);

		if(m_tr->is_infinite(m_hint)) {
			int i = m_hint->index(m_tr->infinite_vertex());
			FT f = m_hint->vertex((i+1)&3)->f() ;
			// if ( f == 999999 )
			if ( invalidValue( f ) )
				return m_maxVal ;
			else
				return f ;
		}

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);

		// Compute the function value at the vertices of the tetrahedron
		// If one of them is a negative value, default to -1 (i.e., mark it as invalid)
		FT f0 = m_hint->vertex(0)->f() ;
		// if ( f0 == 999999 ) {
		if ( invalidValue( f0 ) ) {
			f0 = m_maxVal ;
		}
		// std::cout << "f(0) = " << f0 << std::endl ;
		FT f1 = m_hint->vertex(1)->f() ; 
		// if ( f1 == 999999 ) {
		if ( invalidValue( f1 ) ) {
			f1 = m_maxVal ;
		}
		// std::cout << "f(1) = " << f1 << std::endl ;
		FT f2 = m_hint->vertex(2)->f() ;
		// if ( f2 == 999999 ) {
		if ( invalidValue( f2 ) ) {
			f2 = m_maxVal ;
		}
		// std::cout << "f(2) = " << f2 << std::endl ;
		FT f3 = m_hint->vertex(3)->f() ;
		// if ( f3 == 999999 ) {
		if ( invalidValue( f3 ) ) {
			f3 = m_maxVal ;
		}
		// std::cout << "f(3) = " << f3 << std::endl ;

		// --- Debug (Start) ---
		//if ( f0 == m_maxVal || f1 == m_maxVal || f2 == m_maxVal || f3 == m_maxVal ) {
		//	return m_maxVal ;
		//}
		// --- Debug  (End)  ---

		// --- Debug (Start) ---
		//if ( f0 < 0 || f1 < 0 || f2 < 0 || f3 < 0 ) {
		//	std::cout << "Eval = " << a*f0 + b*f1 + c*f2 + d*f3 << std::endl ;
		//}
		// std::cout << "f = " << a*f0 + b*f1 + c*f2 + d*f3 << std::endl ;
		// --- Debug  (End)  ---

		return a*f0 + b*f1 + c*f2 + d*f3 ;

	}

	// Marching cubes version
	FT operator()(const Point_3& p, bool& ok) const {

		ok = true ;

		m_hint = m_tr->locate(p,m_hint);

		if(m_tr->is_infinite(m_hint)) {
			int i = m_hint->index(m_tr->infinite_vertex());
			FT f = m_hint->vertex((i+1)&3)->f() ;
			// if ( f == 999999 ) {
			if ( this->invalidValue( f ) ) {
				ok = false ;
				return m_maxVal ;
			}
			else
				return f ;
		}

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);

		// Compute the function value at the vertices of the tetrahedron
		// If one of them is a negative value, default to -1 (i.e., mark it as invalid)
		FT f0 = m_hint->vertex(0)->f() ;
		// std::cout << "f(0) = " << f0 << std::endl ;
		// if ( f0 == 999999 ) {
		if ( this->invalidValue( f0 ) ) {
			f0 = m_maxVal ;
		}
		FT f1 = m_hint->vertex(1)->f() ; 
		// std::cout << "f(1) = " << f1 << std::endl ;
		// if ( f1 == 999999 ) {
		if ( this->invalidValue( f1 ) ) {
			f1 = m_maxVal ;
		}
		FT f2 = m_hint->vertex(2)->f() ; 
		// std::cout << "f(2) = " << f2 << std::endl ;
		// if ( f2 == 999999 ) {
		if ( this->invalidValue( f2 ) ) {
			f2 = m_maxVal ;
		}
		FT f3 = m_hint->vertex(3)->f() ; 
		// std::cout << "f(3) = " << f3 << std::endl ;
		// if ( f3 == 999999 ) {
		if ( this->invalidValue( f3 ) ) {
			f3 = m_maxVal ;
		}

		// --- Debug (Start) ---
		//if ( f0 == m_maxVal || f1 == m_maxVal || f2 == m_maxVal || f3 == m_maxVal ) {
		//	return m_maxVal ;
		//}
		// --- Debug  (End)  ---

		// --- Debug (Start) ---
		//if ( f0 < 0 || f1 < 0 || f2 < 0 || f3 < 0 ) {
		//	std::cout << "Eval = " << a*f0 + b*f1 + c*f2 + d*f3 << std::endl ;
		//}
		// --- Debug  (End)  ---

		return a*f0 + b*f1 + c*f2 + d*f3 ;

	}

	FT unsigned_func( const Point_3& p ) const {
		m_hint = m_tr->locate(p,m_hint);

		if(m_tr->is_infinite(m_hint)) {
			int i = m_hint->index(m_tr->infinite_vertex());
			FT f = std::abs( m_hint->vertex((i+1)&3)->f() ) ;
			if ( invalidValue( f ) )
				return m_maxVal ;
			else
				return f ;
		}

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);

		// Compute the function value at the vertices of the tetrahedron
		FT f0 = std::abs( m_hint->vertex(0)->f() ) ;
		if ( invalidValue( f0 ) ) {
			f0 = m_maxVal ;
		}
		FT f1 = std::abs( m_hint->vertex(1)->f() ) ;
		if ( invalidValue( f1 ) ) {
			f1 = m_maxVal ;
		}
		FT f2 = std::abs( m_hint->vertex(2)->f() ) ;
		if ( invalidValue( f2 ) ) {
			f2 = m_maxVal ;
		}
		FT f3 = std::abs( m_hint->vertex(3)->f() ) ;
		if ( invalidValue( f3 ) ) {
			f3 = m_maxVal ;
		}

		return a*f0 + b*f1 + c*f2 + d*f3 ;
	}

	// Note: this version does not take into account the maximum
	// This is done on purpose so that the scale computation detects points outside the band as invalid!
	FT unsigned_func_for_scale( const Point_3& p ) const {
		m_hint = m_tr->locate(p,m_hint);

		if(m_tr->is_infinite(m_hint)) {
			int i = m_hint->index(m_tr->infinite_vertex());
			FT f = std::abs( m_hint->vertex((i+1)&3)->f() ) ;
			return f ;
		}

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);

		// Compute the function value at the vertices of the tetrahedron
		FT f0 = std::abs( m_hint->vertex(0)->f() ) ;
		FT f1 = std::abs( m_hint->vertex(1)->f() ) ;
		FT f2 = std::abs( m_hint->vertex(2)->f() ) ;
		FT f3 = std::abs( m_hint->vertex(3)->f() ) ;
		
		return a*f0 + b*f1 + c*f2 + d*f3 ;
	}


	// --> Weighting functions <--

	/// Gaussian weight function
	static FT weightGaussian( const FT& dist, const FT& h = 1.0 ) {
		return exp( -( ( dist*dist ) / ( h*h ) ) ) ;
	}
	
	// Inverse Gaussian weight function (for edges)
	static FT weightInverseGaussian( const FT& dist, const FT& h = 1.0 ) {
		return 1 - exp( -( ( dist*dist ) / ( h*h ) ) ) ;
	}

	/// Wendland weight function
	static FT weightWendland( const FT& dist, const FT& h = 1.0 ) {
		return ( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*(4*dist / h + 1 ) ;
	}

	/// Inverse Wendland weight function
	static FT weightInverseWendland( const FT& dist, const FT& h = 1.0 ) {
		return 1 - ( ( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*(4*dist / h + 1 ) ) ;
	}

	// Get the triangulation
	boost::shared_ptr<Triangulation> getTriangulationPtr() {
		return m_tr ;
	}

	// Is the computed implicit function correct?
	bool correct() const {
		return m_correct ;
	}

	FT distanceToClosestSplatCenter( const Point_3 &query ) const {
		Neighbor_search search( *m_pTree, query, 1 ) ;
		typename Neighbor_search::iterator it = search.begin() ;
		Splat_3 curSplat = boost::get<1>( it->first ) ;
		return CGAL::sqrt( CGAL::squared_distance( query, curSplat.center() ) ) ;
		// return CGAL::squared_distance( query, curSplat.center() ) ;
	}

	bool allPointsOnRestrictedBand( const Triangle &t, const FT &rad ) const
	{
		// FT val = m_gaussianH * rad ;
		FT val = m_maxVal * rad ;

		for ( int i = 0; i < 3; i++ ) {
			// std::cout << abs( (*this)(t.vertex(i)) ) << " / " << val << std::endl ;
			if ( std::abs( (*this)(t.vertex(i)) ) > val )
				return false ;
		}

		return true ;
	}

	/// Evaluate the function at a given point
	FT eval( const Point_3 &p ) const {
		if( m_inputType == SPLATS )
			return evalSplats(p) ; 
		else 
			return evalPoints(p) ;
	}

	/// Evaluate the function using Splats
	FT evalSplats( const Point_3 &p ) const {
		FT sumWeightedDist = 0.0 ;
		FT sumWeights = 0.0 ;

		// std::cout << "query p = " << p << std::endl ;

		/* K-NN */
		if ( m_k > 0 ) {
			Neighbor_search search( *m_pTree, p, m_k ) ; 
			bool any = false ;
			
			for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
				// Distance measure
				Splat_3 curSplat = boost::get<1>( it->first ) ;
				Point_3 projectedPoint = curSplat.project( p ) ;
				FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;

				// Weighting
				if ( dist < m_gaussianH ) {
					FT weight = this->weightGaussian( dist, m_gaussianH ) / m_k ;

					sumWeightedDist += weight * dist ;
					sumWeights += weight ;
					if ( !any )
						any = true ;
				}
				// Else, weight is almost zero, do not take this point into account...
			
			}

			if ( !any ) {
				// return std::numeric_limits<double>::infinity() ;
				return 999999 ;
			}

		}
		else {
			/* Radial NN */
			std::vector< Point_and_splat > neighs ;
			Sphere_range_query searchSphere( p, m_gaussianH ) ;
			m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;

			if ( neighs.empty() ) {
				// return std::numeric_limits<double>::infinity() ;
				return 999999 ;
			}

			int numNeigh = neighs.size() ;
			/*if ( neighs.size() > 0 )
				std::cout << "Num neighs = " << neighs.size() << std::endl ;*/

			typename std::vector< Point_and_splat >::iterator it ;
			for ( it = neighs.begin(); it!= neighs.end(); ++it ) {
				// Distance measure
				Splat_3 curSplat = boost::get<1>( *it ) ;
				Point_3 projectedPoint = curSplat.project( p ) ;
				FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;

				// Weighting
				// FT weight = this->weightGaussian( dist, m_gaussianH ) ;
				FT weight = this->weightGaussian( dist, m_gaussianH ) / numNeigh ;

			/*	if ( dist > m_gaussianH ) {
					std::cout << "Dist = " << dist << " / " << m_gaussianH << std::endl ;
					std::cout << "Weight = " << weight << std::endl ;
				}*/

				sumWeightedDist += weight * dist ;
				sumWeights += weight ;
			}
		}

		//std::cout << "Eval = " << sumWeightedDist / sumWeights << std::endl ;

		/*if ( ( sumWeightedDist / sumWeights ) > m_gaussianH ) {
			std::cout << "Outside the band! = " << ( sumWeightedDist / sumWeights ) << std::endl ;
		}*/

		return sumWeightedDist / sumWeights ;
	}


	/// Evaluate the function using Points
	FT evalPoints( const Point_3 &p ) const {
		FT sumDist = 0.0 ;
		
		// std::cout << "query p = " << p << std::endl ;

		/* K-NN */
		PtsNeighborSearch search( *m_ptsTree, p, m_k ) ; 

		int k = m_k ;

		FT gaussSq = m_gaussianH * m_gaussianH ;

		for( typename PtsNeighborSearch::iterator it = search.begin(); it != search.end(); ++it) {
			// Distance measure
			if ( it->first != p ) {
				FT sqDist = CGAL::squared_distance( it->first, p ) ;
				sumDist += sqDist ;
			}
			else {
				k-- ;
			}
		}
		
		FT dist = CGAL::sqrt( sumDist / k ) ;

		if ( dist < m_gaussianH ) {	
			return dist ;
		}
		else {
			return 999999 ;
		}

	}



	class EvaluateExactImplicitFunction ;
	friend class EvaluateExactImplicitFunction ;

	class EvaluateExactImplicitFunction {
		// Attributes
		const Self& self ;

		public: 
			EvaluateExactImplicitFunction( const Self& self ) : self( self ) {}

			/// 'ImplicitFunction' interface: evaluates the implicit function at a given 3D query point.
			FT operator()(const Point_3& p) const {
				return self.eval( p ) ;
			}

	} ;
	


private:
	/* Attributes */
	boost::shared_ptr<Triangulation> m_tr ; // Main triangulation holding the precomputed implicit function
	pNNTree m_pTree ;
	pPtsNNTree m_ptsTree ;
	Sphere_3 m_boundingSphere ;
	FT m_gaussianH ;
	mutable Cell_handle m_hint; // last cell found = hint for next search
	FT m_maxVal ;
	FT m_smoothPow ;
	FT m_areaConst ;
	bool m_correct ;
	FT m_stWeight ;
	FT m_functionErrorBound ;
	FT m_radiusEdgeBound ;
	Sphere_3 m_enlargedBoundingSphere ;
	int m_numVertOnBand ;
	int m_k ;
	FT m_bandRadius ;
	FT m_defMax ;
	int m_inputType ;
	enum InputTypes{ SPLATS, POINTS } ;
	int m_bandKnn ;
	FT m_timeImplicit ;
	FT m_timeSigning ;
	FT m_epsilon ;
    bool m_debugOutput;

	/* Functions */

	void barycentric_coordinates(	const Point_3& p,
									Cell_handle cell,
									FT& a,
									FT& b,
									FT& c,
									FT& d) const
	{
		const Point_3& pa = cell->vertex(0)->point();
		const Point_3& pb = cell->vertex(1)->point();
		const Point_3& pc = cell->vertex(2)->point();
		const Point_3& pd = cell->vertex(3)->point();

		FT v = volume(pa,pb,pc,pd);
		a = CGAL::abs(volume(pb,pc,pd,p) / v);
		b = CGAL::abs(volume(pa,pc,pd,p) / v);
		c = CGAL::abs(volume(pb,pa,pd,p) / v);
		d = CGAL::abs(volume(pb,pc,pa,p) / v);
	}


	/// Delaunay refinement (break bad tetrahedra, where
	/// bad means badly shaped or too big). The normal of
	/// Steiner points is set to zero.
	/// Returns the number of vertices inserted.
	unsigned int delaunay_refinement( FT radius_edge_ratio_bound, ///< radius edge ratio bound (ignored if zero)
									  FT cell_radius_bound, ///< cell radius bound (ignored if zero)
									  FT functionErrorBound,
									  unsigned int max_vertices, ///< number of vertices bound
									  FT enlarge_ratio ) ///< bounding box enlarge ratio
	{
		// CGAL_TRACE("Calls delaunay_refinement(radius_edge_ratio_bound=%lf, cell_radius_bound=%lf, max_vertices=%u, enlarge_ratio=%lf)\n",
					//radius_edge_ratio_bound, cell_radius_bound, max_vertices, enlarge_ratio);

		EvaluateExactImplicitFunction funcExact( *this ) ; 
		m_enlargedBoundingSphere = enlarged_bounding_sphere(enlarge_ratio);
		unsigned int nb_vertices_added = unsigned_distance_refine_triangulation( *m_tr,
																				 radius_edge_ratio_bound,
																				 cell_radius_bound,
																				 max_vertices,
																				 m_enlargedBoundingSphere ) ;

		// CGAL_TRACE("End of delaunay_refinement()\n");

		return nb_vertices_added;
	}


	/// Gets median value of the implicit function over input vertices.
	FT median_value_at_input_vertices() const
	{
		std::deque<FT> values;
		Finite_vertices_iterator v, e;
		for(v = m_tr->finite_vertices_begin(),
			e= m_tr->finite_vertices_end();
			v != e; 
			v++)
			if(v->type() == Triangulation::INPUT)
			values.push_back(v->f());

		int size = values.size();
		if(size == 0)
		{
			std::cerr << "Contouring: no input points\n";
			return 0.0;
		}

		std::sort(values.begin(),values.end());
		int index = size/2;
		// return values[size/2];
		return 0.5 * (values[index] + values[index+1]); // avoids singular cases
	}



	/// Gets minimum value of the implicit function over input vertices.
	FT minimum_value_at_input_vertices() const
	{
		std::deque<FT> values;
		Finite_vertices_iterator v, e;
		FT minVal = std::numeric_limits<double>::infinity() ;
		for(v = m_tr->finite_vertices_begin(),
			e= m_tr->finite_vertices_end();
			v != e; 
			v++)
		{
			if ( minVal > v->f() ) {
				minVal = v->f() ;
			}
		}

		return minVal ;
	}



	FT edge_weight( Point_3 p1, Point_3 p2, int steps = 5 ) {
		// Edge joining the two points
		Segment_3 s( p1, p2 ) ;
		FT stepSize = CGAL::sqrt( s.squared_length() ) / ( steps-1 ) ;
		Vector_3 dirNorm = s.to_vector() / CGAL::sqrt( s.to_vector().squared_length() ) ;

		FT val = 0 ;
		for ( int i = 0; i < steps; i++ ) {
			val += (*this)( p1 + ( dirNorm * stepSize * i ) ) ;
		}
		
		return pow( ( val / steps ), m_smoothPow ) + m_areaConst ;
	}



	FT min_along_edge( Point_3 p1, Point_3 p2, int steps = 5 ) {
		// Edge joining the two points
		Segment_3 s( p1, p2 ) ;
		FT stepSize = CGAL::sqrt( s.squared_length() ) / ( steps-1 ) ;
		Vector_3 dirNorm = s.to_vector() / CGAL::sqrt( s.to_vector().squared_length() ) ;

		FT val = 999999999999999 ;
		for ( int i = 0; i < steps; i++ ) {
			FT t = (*this)( p1 + ( dirNorm * stepSize * i ) ) ;
			if ( t < val ) 
				val = t ;
		}
		
		return val ;
	}
	


	/// Shift the surface contouring value
	void shift_f(const FT shift)
	{
		Finite_vertices_iterator v, e;
		for(v = m_tr->finite_vertices_begin(),
			e = m_tr->finite_vertices_end();
			v!= e;
			v++) 
		{
			if ( !invalidValue( v->f() ) )
				v->f() += shift;
		}
	}

	/// Unitize the f value (all values will be between 0 and 1)
	void unitize_f()
	{
		// Get minimum (maximum is supposed to be already computed)
		FT min = minimum_value_at_input_vertices() ;
		
		Finite_vertices_iterator v, e;
		for(v = m_tr->finite_vertices_begin(),
			e = m_tr->finite_vertices_end();
			v!= e;
			v++) 
		{
			if ( !invalidValue( v->f() ) )
				v->f() = ( v->f() - min ) / ( m_maxVal - min ) ;
		}
	}



	/// Builds a graph from the
	void build_graph( Graph& g, std::map< Vertex_handle, typename Graph::vertex_descriptor >& tri2graphVert ) {
		
		std::cout << "    - Building graph..." ;
		int verbosePercent = 0 ; int verboseLastPercent = 0 ; int i = 0 ;

		CGAL::Timer task_timer; task_timer.start() ;
		int k = 0 ;

		Finite_vertices_iterator itv ;
		int numAddedVert = 0 ;
		for ( itv = m_tr->finite_vertices_begin();
			  itv != m_tr->finite_vertices_end();
			  ++itv )
		{
			if ( itv->ind() >= 0 ) {
				typename Graph::vertex_descriptor v = boost::add_vertex( g ) ;
				g[v].triVertHandle = itv ;
				tri2graphVert[ itv ] = v ;
				g[v].ind = -1 ;
				numAddedVert++ ;
			}
		}
		if ( numAddedVert != m_numVertOnBand ) 
			std::cout << "[ERROR] The number of added vertices must match the vertices on the band --> " << numAddedVert << " - " << m_numVertOnBand << std::endl ;


		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
		Finite_edges_iterator it ;
		int numEdges = m_tr->number_of_edges() ;
		for ( it = m_tr->finite_edges_begin();
			  it !=	m_tr->finite_edges_end();
			  ++it, i++ )
		{
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numEdges ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}

			Cell_handle c = (it->first) ;
			Vertex_handle v0 = c->vertex( it->second ) ;
			Vertex_handle v1 = c->vertex( it->third ) ;

			if ( v0->ind() >= 0 && v1->ind() >= 0 ) {
				// if ( v0->f() == 999999 || v1->f() == 999999 ) {
				if ( this->invalidValue( v0->f() ) || this->invalidValue( v1->f() ) ) {
					if ( m_bandKnn > 0 )
						std::cout << "[ERROR] Should never enter here!!" << std::endl ;
					// triplets.push_back( Eigen::Triplet<double>( v0->ind(), v1->ind(), m_defMax ) ) ;
					// triplets.push_back( Eigen::Triplet<double>( v1->ind(), v0->ind(), m_defMax ) ) ;
					boost::add_edge( tri2graphVert[ v0 ], tri2graphVert[ v1 ], m_defMax, g ) ;
					boost::add_edge( tri2graphVert[ v1 ], tri2graphVert[ v0 ], m_defMax, g ) ;
					// std::cout << "MAX = " << m_defMax << std::endl ;
				}
				else {
					FT val = ( v0->f() + v1->f() ) / 2 ;
					// val = val / m_maxVal ; // TEST!!!
					// val = this->weightInverseWendland( val, m_maxVal ) ;
					// FT w = pow( wf, m_smoothPow ) ;
					FT w = pow( val, m_smoothPow ) ;
										
					boost::add_edge( tri2graphVert[ v0 ], tri2graphVert[ v1 ], w, g ) ;
					boost::add_edge( tri2graphVert[ v1 ], tri2graphVert[ v0 ], w, g ) ;

					// std::cout << "val = " << val << std::endl ;
				}
			}
		}
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << " s, " << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated)" << std::endl ;
	}



	/// Filter the graph representation so it just contains the largest connected component
	FilteredGraph keep_largest_connected_component( Graph& g ) {

		ComponentMap comp = boost::get( &VertexProperties::component, g ) ;

		// int numCC = connected_components( g, boost::get( &VertexProperties::component, g ) ) ;
		int numCC = boost::connected_components( g, comp ) ;
		std::cout << "    - Number of Connected Components = " << numCC << std::endl ;

		// Search for the largest CC
		std::vector< int > numElemsCC ;
		numElemsCC.reserve( numCC ) ;
		for ( int k = 0; k < numCC; k++ ) {
			numElemsCC.push_back(0) ;
		}

		typedef typename boost::graph_traits< Graph >::vertex_iterator Bgl_vertex_iterator ;
		Bgl_vertex_iterator vi, vi_end ;
		for ( boost::tie( vi, vi_end ) = boost::vertices( g ); 
			  vi != vi_end; 
			  ++vi ) 
		{
			numElemsCC.at( g[*vi].component )++;
		}
		std::vector<int>::const_iterator itMax;
		itMax = std::max_element( numElemsCC.begin(), numElemsCC.end() ) ;
		int largestCCInd = itMax - numElemsCC.begin() ;

		// Filter the part of the map corresponding to the largest connected component
		vertexComponent<ComponentMap> vfilter( boost::get( &VertexProperties::component, g ), largestCCInd ) ;
		FilteredGraph fg( g, boost::keep_all(), vfilter ) ;

		return fg ;
	}


	SparseMatrix build_laplacian( FilteredGraph& fg, Graph& g, std::map< Vertex_handle, typename Graph::vertex_descriptor >& tri2graphVert ) const 
	{
		// Build the Laplacian matrix
		std::cout << "    - Building Laplacian matrix..." << std::flush ;
		CGAL::Timer task_timer; task_timer.start() ;

		typedef typename boost::graph_traits < FilteredGraph >::vertex_iterator Bgl_filtered_vertex_iterator ;
		Bgl_filtered_vertex_iterator fvit, fvit_end ;
		int i = 0 ;
		for ( boost::tie( fvit, fvit_end ) = boost::vertices( fg );
			  fvit != fvit_end; 
			  ++fvit, i++ ) 
		{
			fg[*fvit].ind = i ;
		}
		int numVert = i ;

		std::vector< Eigen::Triplet<double> > triplets ;
		typename boost::property_map< FilteredGraph, boost::edge_weight_t >::type EdgeWeightMap = boost::get( boost::edge_weight, fg ) ;
		typedef typename boost::graph_traits< FilteredGraph >::edge_iterator Bgl_filtered_edge_iterator ;
		Bgl_filtered_edge_iterator feit, feit_end ;
		for ( boost::tie( feit, feit_end ) = boost::edges( fg );
			  feit != feit_end; 
			  ++feit ) 
		{
			typename Graph::vertex_descriptor s = boost::source( *feit, fg ) ;
			typename Graph::vertex_descriptor t = boost::target( *feit, fg ) ;
			FT w = boost::get( EdgeWeightMap, *feit ) ;

			// --- Debug (Start) ---
			// Is the edge part of the original triangulation?
			/*Cell_handle ch ;
			int da, db ;
			if ( !m_tr->is_edge( g[s].triVertHandle, g[t].triVertHandle, ch, da, db ) ) {
				std::cout << "[ERROR] This edge does not exist in the triangulation!" << std::endl ;
			}*/
			// std::cout << g[s].ind << ", " << g[t].ind << "=" << w << std::endl ;
			// --- Debug (Start) ---
			
			triplets.push_back( Eigen::Triplet<double>( g[s].ind, g[t].ind, w ) ) ;
			triplets.push_back( Eigen::Triplet<double>( g[t].ind, g[s].ind, w ) ) ;
		}

		SparseMatrix W( numVert, numVert ) ;
		W.setFromTriplets( triplets.begin(), triplets.end() ) ;

		std::cout << "done (" << task_timer.time() << " s, " << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated)" << std::endl ;

		return W ;
	}


	/// Apply normalized cut to the computed graph Laplacian
	bool normalized_cut( SparseMatrix& W, Graph& g, std::map< Vertex_handle, typename Graph::vertex_descriptor >& tri2graphVert ) const
	{
		std::cout << "    - Normalized cut..." << std::flush ; 
		CGAL::Timer task_timer; task_timer.start() ;
		
		EMatrix eigenVectors;
		EVector eigenValues;
		task_timer.reset() ;
		ncut<FT>( W, 2, eigenVectors, eigenValues ) ;
		std::cout << "done (" << task_timer.time() << " s)" << std::endl ;
		std::cout << "[DEBUG] EigenValues = < " ;
		for ( int k = 0 ; k < 2; k++ ) 
			std::cout << eigenValues(k) << ", " ; 
		std::cout << " >" << std::endl ;
		
		// --- Debug (Start) ---
        std::ofstream ofio, ofoo, ofe, oft, ofb;
		if (m_debugOutput) {
            ofio.open("./_OUT/InsideVertices_OPT.xyz", std::ios_base::out);
            ofoo.open("./_OUT/OutsideVertices_OPT.xyz", std::ios_base::out);
            ofe.open("./_OUT/2ndSmallestEigenVector.txt", std::ios_base::out);
            oft.open("./_OUT/TriVert.xyz", std::ios_base::out);
            ofb.open("./_OUT/VertInBand.xyz", std::ios_base::out);
            // std::ofstream off( "./_OUT/DistanceFunction.xyz", std::ios_base::out ) ;
            /*std::ofstream oftt( "./_OUT/Triangulation.tri", std::ios_base::out ) ;
            oftt << *m_tr ;
            oftt.close() ;*/
        }
		// --- Debug  (End)  ---

		// Sign function according to second smallest eigenvector
		std::cout << "    - Signing function..." << std::flush ; 
		task_timer.reset() ;
		EVector secondEigenVector = eigenVectors.col(0) ;
		// std::cout << "[DEBUG] 2nd largest EigenVector = < " << secondEigenVector << " >" << std::endl ;
		Finite_vertices_iterator v, e;
		bool any = false ; 
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v )
		{
			oft << v->point() << std::endl ;
			if ( v->ind() >= 0 ) {
				typename Graph::vertex_descriptor vd = tri2graphVert[ v ] ;
				if ( g[vd].ind >= 0 ) {
					ofb << v->point() << " " << v->f() << std::endl ;
					ofe << secondEigenVector( g[vd].ind ) << std::endl ;
					if ( secondEigenVector( g[vd].ind ) < 0 ) {
						// Vertex_handle vh = fg[vd].triVertHandle ;
						v->f() = -v->f() ;
						ofio << v->point() << std::endl ;
						if ( !any ) any = true ;
					}
					else {
						ofoo << v->point() << std::endl ;
					}
				}
			}
		}
		// --- Debug (Start) ---
		ofio.close() ;
		ofoo.close() ;
		oft.close() ;
		ofe.close() ;
		ofb.close() ;
		// --- Debug  (End)  ---
		std::cout << "done (" << task_timer.time() << ")" << std::endl ;
		
		return any ;
	}


	/// Sign the implicit function inside the band
	bool sign_implicit_function_band_ncut() 
	{
		std::map< Vertex_handle, typename Graph::vertex_descriptor > tri2graphVert ;
		Graph g ;
		build_graph( g, tri2graphVert ) ;
		FilteredGraph fg = keep_largest_connected_component( g ) ;
		SparseMatrix W = build_laplacian( fg, g, tri2graphVert ) ;

		return normalized_cut( W, g, tri2graphVert ) ;
	}

	bool invalidValue( const FT& val ) const {
		// const FT epsilon = 0.01 ;
		// std::cout << "std::abs( val - 999999.0 ) = " << std::abs( val - 999999.0 ) << std::endl ;
		return std::abs( val - UNDEFINED_VALUE ) < m_epsilon ;
	}

} ;

#endif // SPLATSUNSIGNEDDISTANCEFUNCTIONNCUT_H
