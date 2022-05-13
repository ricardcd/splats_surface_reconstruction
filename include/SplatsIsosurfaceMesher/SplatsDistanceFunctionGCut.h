#ifndef SPLATSUNSIGNEDDISTANCEFUNCTION_H
#define SPLATSUNSIGNEDDISTANCEFUNCTION_H

/* Includes */
// CGAL
#include <CGAL/trace.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
	// AABB_tree
#include "Splat_3/AABBTree/AABB_Splat_3_primitive.h"
#include "Splat_3/AABBTree/AABB_Splat_3_intersections.h"
#include "Splat_3/AABBTree/AABB_Splat_3_traits.h"
#include <CGAL/AABB_tree.h>
	// Spatial searching
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include "SplatsCreation/KDTreeSearch/My_Search_traits_3.h"
#include "SplatsCreation/KDTreeSearch/Fuzzy_capsule_3.h"
	// Triangulation (refinement)
#include <SplatsIsosurfaceMesher/Splats_distance_function_triangulation_3.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Robust_circumcenter_filtered_traits_3.h>
	// Random number generator
 #include <CGAL/Random.h>
	// Other
#include <CGAL/compute_average_spacing.h>
#include <CGAL/linear_least_squares_fitting_3.h>
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
#include <boost/filesystem.hpp>
// Graph Cuts
#include <graph.h>
// Std
#include <algorithm>
// Undefined value define
#ifndef UNDEFINED_VALUE
	#define UNDEFINED_VALUE 999999.0
#endif


template < class K >
class SplatsDistanceFunctionGCut {

public:
	/* Typedefs */
	typedef	SplatsDistanceFunctionGCut		Self ;

	typedef typename K::FT												FT ;
	typedef typename K::Point_3											Point_3 ;
	typedef typename K::Sphere_3										Sphere_3 ;
	typedef typename K::Ray_3											Ray_3 ;
	typedef typename K::Vector_3										Vector_3 ;
	typedef typename CGAL::Splat_3< K >									Splat_3 ;
	typedef typename K::Segment_3										Segment_3 ;
	typedef typename K::Plane_3											Plane_3 ;
		
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
    typedef typename std::vector< Splat_3 >::const_iterator				SplatsIterator ;
    typedef typename CGAL::AABB_Splat_3_primitive< K, SplatsIterator >  Primitive ;
    typedef typename CGAL::AABB_Splat_3_traits< K, Primitive >			AABB_Splat_3_traits ;
    typedef typename CGAL::AABB_tree< AABB_Splat_3_traits >				SplatsAABBTree ;
    typedef typename SplatsAABBTree::Object_and_primitive_id			Object_and_primitive_id;
    typedef typename boost::shared_ptr< SplatsAABBTree >				pSplatsAABBTree ;

    // Octree type
    typedef typename CGAL_Octree::CGAL_Octree< K >                      Octree ;
    typedef typename CGAL_Octree::CubicBounds< K >                      CubicBounds ;

	/* Nearest neighbors search definitions */
	// Property map, each point takes as 
	typedef boost::tuple< Point_3, Splat_3 >							Point_and_splat;

	class My_point_property_map{
    public:
		typedef Point_3 value_type ;
		typedef const typename My_point_property_map::value_type& reference ;
		typedef const Point_and_splat& key_type ;
		typedef boost::readable_property_map_tag category;

		// Get function for the property map
		friend reference get( My_point_property_map mppm,
							  typename My_point_property_map::key_type p ) {
			return boost::get<0>(p);
		}

	} ;

	typedef typename CGAL::Search_traits_3< K >							TreeTraitsBase ;
	typedef CGAL::Search_traits_adapter< Point_and_splat,
										 My_point_property_map,
										 TreeTraitsBase>				TreeTraits ;
	typedef typename CGAL::Kd_tree<TreeTraits>							NNTree ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< TreeTraits >	Neighbor_search ;
	typedef boost::shared_ptr<NNTree>									pNNTree ;		
	typedef typename CGAL::Fuzzy_sphere< TreeTraits >					Sphere_range_query ;
	
	// S-T Cut Graph type
	typedef Graph< double, double, double >								GraphType ;



	/* Functions */
	/// Constructor, computes the (unsigned) implicit function in a refined Delaunay triangulation build with the centers of the splats
	/// Also signs the function using Graph Cuts
	SplatsDistanceFunctionGCut( const std::vector< Splat_3 >& splats,
							const int& k = 25,
							const FT& gaussianRadiusFactor = 0.2, 
							const FT& ransacThres = 0.1,
							const int rayTries = 10,
							const FT confidenceThreshold = 0.7,
							const FT uncertainBand = 1.0,
							const FT smoothPrior = 4.0, 
							const FT areaConst = 1e-5,
							const FT stWeight = 1e5, 
							const FT dataWeight = 1.0, 
							const FT functionErrorBound = 1e-3,
							const FT radiusEdgeBound = 2.5, 
							const bool boundedSurface = false,
							const double bandRadius = 0.0,
							const int bandKnn = -1,
							const double epsilon = 0.001,
                            const bool signPointsOutsideBand = false,
                            const int odepth = -1 )
	: m_tr(new Triangulation) 
	{
		m_k = k ;
		m_ransacThres = ransacThres ;
		m_rayTries = rayTries ;
		m_confidenceThres = confidenceThreshold ;
		m_uncertainBand = uncertainBand ;
		m_smoothPrior = smoothPrior ;
		m_areaConst = areaConst ;
		m_correct = true ; // By default, we assume no error...
		m_stWeight = stWeight ;
		m_radiusEdgeBound = radiusEdgeBound ;
		m_dataWeight = dataWeight ;
		m_boundedSurface = boundedSurface ;
		m_timeImplicit = 0.0 ;
		m_timeSigning = 0.0 ;
		m_bandKnn = bandKnn ;
		m_epsilon = epsilon ;
		m_signPointsOutsideBand = signPointsOutsideBand ;
		
		FT factor = gaussianRadiusFactor ;
		if ( m_k > 0 && ( gaussianRadiusFactor < 0 || gaussianRadiusFactor > 1 ) ) {
			std::cerr << "[WARNING] gaussianRadiusFactor must be between 0 and 1! Defaulting to gaussianRadiusFactor=0.2." << std::endl ;
			factor = 0.2 ;
		}
		
		// --- Debug (Start) ---
		// Initialize output directory
		boost::filesystem::path dir( "./_OUT" ) ;
		if (boost::filesystem::create_directory( dir ) )
			std::cout << "[DEBUG] Created Debug Output Directory (./_OUT)" << "\n";
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

		
		// Construct the KD-Tree with these center points
		m_pTree = pNNTree( new NNTree( boost::make_zip_iterator(boost::make_tuple( points.begin(), splats.begin() ) ),
									   boost::make_zip_iterator(boost::make_tuple( points.end(), splats.end() ) ) ) ) ;

		// Construct the AABB tree for intersection queries
		m_splatsTree = pSplatsAABBTree( new SplatsAABBTree( splats.begin(), splats.end() ) ) ;
		
		// Computes average spacing
		// FT average_spacing = CGAL::compute_average_spacing(points.begin(), points.end(), 6 /* knn = 1 ring */ ) ;

		// Compute gaussian weight H variable
		/*if ( m_k > 0 ) {
			m_gaussianH = CGAL::sqrt( factor * factor * this->bounding_sphere().squared_radius() ) ;
			m_bandRadius = CGAL::sqrt( bandRadius * bandRadius * this->bounding_sphere().squared_radius() ) ;
		}
		else {
			m_gaussianH = factor * average_spacing ;
			m_bandRadius = bandRadius * average_spacing ;
		}*/
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

		
		// Error bound for linear approximation (w.r.t. the boundig sphere radius)
		//if ( functionErrorBound > 0 ) 
		//	m_functionErrorBound = CGAL::sqrt( functionErrorBound * functionErrorBound * this->bounding_sphere().squared_radius() ) ;
		//else
			m_functionErrorBound = functionErrorBound ;
		
		// m_functionErrorBound = functionErrorBound ;

		m_splats = splats ;

		// --- Debug (Start) ---
		std::cout << "  - Some Debug Info:" << std::endl ;
		std::cout << "    Function error bound = " << m_functionErrorBound << std::endl ;
		std::cout << "    average_spacing = " << average_spacing << std::endl ;
		std::cout << "    m_gaussianH     = " << m_gaussianH << std::endl ;
		// --- Debug  (End)  ---

		CGAL::Timer timer ;

		// Compute the (unsigned) implicit function
		timer.start() ;
		m_correct = compute_implicit_function() ;
		m_timeImplicit = timer.time() ;
		if ( !m_correct ) return ;

		// Prepare the trick, if needed
		if ( m_boundedSurface ) {
			// prepare_open_surface_trick( points ) ;
			FT q = linear_least_squares_fitting_3( points.begin(),
												   points.end(),
												   m_sphericalCapPlane,
												   CGAL::Dimension_tag<0>() ) ;
			m_sphericalCapPlane = m_sphericalCapPlane.opposite() ;
		}

		// Sign the function
		timer.reset() ;
		m_correct = sign_implicit_function() ;
		m_timeSigning = timer.time() ;
		if ( !m_correct ) return ;

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

	}


	
	/// Compute the implicit function at the vertices of a refined triangulation (similar to an octree decomposition)
	bool compute_implicit_function() {
		
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
		std::cout << "    - Refining..."  ;

		unsigned int nb_vertices_added = delaunay_refinement( radius_edge_ratio_bound,
															  cell_radius_bound,
															  max_vertices,
															  enlarge_ratio ) ;

		// Prints status
		std::cout << "done, added " << nb_vertices_added << " Steiner points in "
						  << task_timer.time() << " seconds, "
						  << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
						  << std::endl ;
		task_timer.reset() ;

		std::cout << "  - Compute implicit function..." ;
		int verbosePercent = 0 ; int verboseLastPercent = 0 ; int i = 0 ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;

		// Computes the Implicit function eval()
		// at each vertex of the triangulation.
		m_maxVal = 0 ;
		int numVert = m_tr->number_of_vertices() ;
		Finite_vertices_iterator v, e;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++)
		{
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numVert ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}

			// Most function values are already computed in the refinement process...
			// However, some of them may remain unknown
			if ( v->f() < 0 )
				v->f() = eval( v->point() ) ;

			// Secondary band using knn
			if ( this->invalidValue( v->f() ) )
				v->f() = secondaryBandKnn( v->point() ) ;

			// Assign also a linear index
			v->ind() = i ; 

			if ( m_maxVal < v->f() && !this->invalidValue( v->f() ) ) {
				m_maxVal = v->f() ;
			}
		}

		// Prints status
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << " s)" << std::endl ;
		task_timer.reset() ;

		// --- Debug (Start) ---
		std::cout << "[DEBUG] Implicit function max val = " << m_maxVal << std::endl ;
		// Compute the value of the function at each of its vertices
		// FT mv = median_value_at_input_vertices() ;
		// std::cout << "[DEBUG] Median value at input points: " << mv << std::endl ;
		// --- Debug  (End)  ---

		return true;

	}


	bool sign_implicit_function() {

		// --- Debug (Start) ---
		if ( m_maxVal < m_uncertainBand ) {
			std::cout << "[WARNING] Maximum value computed on the implicit function is lower than the uncertain band..." << std::endl ;
			m_uncertainBand = m_maxVal ;
		}
		// --- Debug  (End)  ---

		// Create graph structure
		
		GraphType g( /*estimated # of nodes*/ m_tr->number_of_vertices(), /*estimated # of edges*/ m_tr->number_of_edges() ) ; 
		g.add_node( m_tr->number_of_vertices() ) ;

		// Fill graph
		if ( m_signPointsOutsideBand )
			computeSTWeights( g ) ;
		else
			computeSTWeightsRestricted( g ) ;
		computeSmoothWeights( g ) ;
		
		std::cout << "    - Graph cut..." << std::flush ;
		CGAL::Timer task_timer; task_timer.start() ;
		int flow = g.maxflow();
		std::cout << "done (" << task_timer.time() << "s )" << std::endl ;

		// --- Debug (Start) ---
		std::ofstream ofio( "./_OUT/InsideVertices_OPT.xyz", std::ios_base::out ) ;
		std::ofstream ofoo( "./_OUT/OutsideVertices_OPT.xyz", std::ios_base::out ) ;
		// --- Debug  (End)  ---

		Finite_vertices_iterator v, e;
		bool any1 = false ;
		bool any2 = false ;
		int i ;
		for( i = 0, v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++)
		{
			if ( g.what_segment(i) == GraphType::SINK ) {
				v->f() = -v->f() ;
				if (!any1) any1 = true ;
				// --- Debug (Start) ---
				ofio << v->point() << std::endl ;
				// --- Debug  (End)  ---
			}
			// --- Debug (Start) ---
			else {
				if (!any2) any2 = true ;
				ofoo << v->point() << std::endl ;
			}
			// --- Debug  (End)  ---
		}
		// --- Debug (Start) ---
		ofio.close() ;
		ofoo.close() ;
		// --- Debug  (End)  ---

		// Check correctness of signing
		m_correct = any1 && any2 ;
		
		return m_correct ;
	}


	/// 'ImplicitFunction' interface: evaluates the implicit function at a given 3D query point.
	FT operator()(const Point_3& p) const {

		m_hint = m_tr->locate(p,m_hint);

		if(m_tr->is_infinite(m_hint)) {
			int i = m_hint->index(m_tr->infinite_vertex());
			FT f = m_hint->vertex((i+1)&3)->f() ;
			if ( this->invalidValue( f ) )
				return m_maxVal ;
			else
				return f ;
		}

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);

		// Compute the function value at the vertices of the tetrahedron
		// If one of them is a negative value, default to -1 (i.e., mark it as invalid)
		FT f0 = m_hint->vertex(0)->f() ;
		// std::cout << "f(0) = " << f0 << std::endl ;
		if ( this->invalidValue( f0 ) ) {
			f0 = m_maxVal ;
		}
		FT f1 = m_hint->vertex(1)->f() ; 
		// std::cout << "f(1) = " << f1 << std::endl ;
		if ( this->invalidValue( f1 ) ) {
			f1 = m_maxVal ;
		}
		FT f2 = m_hint->vertex(2)->f() ; 
		// std::cout << "f(2) = " << f2 << std::endl ;
		if ( this->invalidValue( f2 ) ) {
			f2 = m_maxVal ;
		}
		FT f3 = m_hint->vertex(3)->f() ; 
		// std::cout << "f(3) = " << f3 << std::endl ;
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
	
	/// Inverse Gaussian weight function (for edges)
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
	bool correct() {
		return m_correct ;
	}

	bool allPointsOnRestrictedBand( const Triangle &t, const FT &rad ) const
	{
		// FT val = m_gaussianH * rad ;
		FT val = m_maxVal * rad ;

		for ( int i = 0; i < 3; i++ ) {
			//std::cout << abs( (*this)(t.vertex(i)) ) << " / " << val << std::endl ;
			if ( std::abs( (*this)(t.vertex(i)) ) > val )
				return false ;
			/*else if ( distanceToClosestSplatCenter( t.vertex(i) ) > val )
				return false ;*/
		}

		return true ;
	}

	FT distanceToClosestSplatCenter( const Point_3 &query ) const {
		Neighbor_search search( *m_pTree, query, 1 ) ;
		typename Neighbor_search::iterator it = search.begin() ;
		Splat_3 curSplat = boost::get<1>( it->first ) ;
		return CGAL::sqrt( CGAL::squared_distance( query, curSplat.center() ) ) ;
		// return CGAL::squared_distance( query, curSplat.center() ) ;
	}

	/// Evaluate the function at a given point
	FT eval( const Point_3 &p ) const {

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
					FT weight = this->weightGaussian( dist, m_gaussianH ) ;

					sumWeightedDist += weight * dist ;
					sumWeights += weight ;
					if ( !any )
						any = true ;
				}
				// Else, weight is almost zero, do not take this point into account...
			
			}

			if ( !any ) {
				// return std::numeric_limits<double>::infinity() ;
				return UNDEFINED_VALUE ;
			}
		}
		else {
			///* Radial NN */
			//std::vector< Point_and_splat > neighs ;
			//Sphere_range_query searchSphere( p, m_gaussianH ) ;
			//m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;

			//if ( neighs.empty() ) {
			//	// return std::numeric_limits<double>::infinity() ;
			//	return UNDEFINED_VALUE ;
			//}

			////if ( neighs.size() > 1 )
			////	std::cout << "Num neighs = " << neighs.size() << std::endl ;

			//typename std::vector< Point_and_splat >::iterator it ;
			//for ( it = neighs.begin(); it!= neighs.end(); ++it ) {
			//	// Distance measure
			//	Splat_3 curSplat = boost::get<1>( *it ) ;
			//	Point_3 projectedPoint = curSplat.project( p ) ;
			//	FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;

			//	// Weighting
			//	FT weight = this->weightGaussian( dist, m_gaussianH ) ;

			//	sumWeightedDist += weight * dist ;
			//	sumWeights += weight ;
			//}

			/* Radial NN */
			std::vector< Point_and_splat > neighs ;
			Sphere_range_query searchSphere( p, m_gaussianH ) ;
			m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;

			if ( neighs.empty() ) {
				// return std::numeric_limits<double>::infinity() ;
				return UNDEFINED_VALUE ;
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

				//if ( dist > m_gaussianH ) {
				//	std::cout << "Dist = " << dist << " / " << m_gaussianH << std::endl ;
				//	std::cout << "Weight = " << weight << std::endl ;
				//}

				sumWeightedDist += weight * dist ;
				sumWeights += weight ;
			}
		}

		////if ( neighs.size() > 1 )
		//std::cout << "Eval = " << sumWeightedDist / sumWeights << std::endl ;

		return sumWeightedDist / sumWeights ;
	}

	/// Get the time spent creating the unsigned distance function
	FT getTimeImplicitFunctionCreation() {
		return m_timeImplicit ;
	}

	/// Get the time spent signing the distance function
	FT getTimeSigningFunction() {
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
		return Sphere_3(bsphere.center(), bsphere.squared_radius() * ratio*ratio);
	}
	
	std::vector< Point_3 > getIntersections( Segment_3 & seg ) const {
		std::vector< Point_3 > intersections ;
		if ( m_splatsTree->do_intersect( seg ) ) {
			std::vector< Object_and_primitive_id > obj_prim_vect ;
			m_splatsTree->all_intersections( seg, std::back_inserter( obj_prim_vect ) ) ;
			typename std::vector< Object_and_primitive_id >::iterator it ;
			// Get the intersection points and its corresponding Splat_3 primitives
			for ( it = obj_prim_vect.begin(); it != obj_prim_vect.end(); ++it ) {
				Point_3 p ;
				if ( CGAL::assign( p, it->first ) ) {
					intersections.push_back( p ) ;
				}
			}
		}
		return intersections ;
	}

	// Compute the intersections between splats and a ray using RANSAC
	std::vector< Point_3 > getIntersectionsRANSAC( Ray_3& ray ) const {
		std::vector< Point_3 > allIntersections ;
				
		if ( m_splatsTree->do_intersect( ray ) ) {
			std::vector< Object_and_primitive_id > obj_prim_vect ;
			m_splatsTree->all_intersections( ray, std::back_inserter( obj_prim_vect ) ) ;
			typename std::vector< Object_and_primitive_id >::iterator it ;

			// Get the intersection points and its corresponding Splat_3 primitives
			for ( it = obj_prim_vect.begin(); it != obj_prim_vect.end(); ++it ) {
				Point_3 p ;
				if ( CGAL::assign( p, it->first ) ) {

					//if ( (*this)(p) > 0 ) 
					allIntersections.push_back( p ) ;
				}
			}
		}
		
		// std::cout << "----- All Intersect = " << allIntersections.size() << std::endl ;

		// Apply RANSAC-based intersection detection
		CentroidModelEstimator *estimator = new CentroidModelEstimator( m_ransacThres ) ;

		/* Compute RANSAC Local Smooth Surface */
		std::vector< Point_3 > intersections ;
		bool good = true ;
		do {
			std::vector< typename std::vector< Point_3 >::const_pointer > inliersPtr ;	
			std::vector< typename std::vector< Point_3 >::const_pointer > outliersPtr ;	
			std::vector< double > modelParam ;
			int usedData = Ransac< Point_3, double >::run(	estimator, 
											 				allIntersections,									  
															modelParam,
															inliersPtr,
															outliersPtr,
															0.99,
															100,
															100 ) ;

			if ( inliersPtr.size() >= 2 ) {
				// Store new point
				Point_3 pint( modelParam[0], modelParam[1], modelParam[2] ) ;
				intersections.push_back( pint ) ;
				
				std::vector< Point_3 > outliers ;
				for ( int j = 0; j < outliersPtr.size(); j++ ) {
					outliers.push_back( *outliersPtr[ j ] ) ;
				}

				// std::cout << "inliers.size() = " << inliersPtr.size() << std::endl ;
				// std::cout << "outliers.size() = " << outliers.size() << std::endl ;

				allIntersections= outliers ;
				
				// std::cout << "allIntersections.size() = " << allIntersections.size() << std::endl ;
			}
			else {
				good = false ;
			}

		} while( good ) ;

		// std::cout << "----- Num. Intersect = " << intersections.size() << std::endl ;

		return intersections ;

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
	std::vector< Splat_3 > m_splats ;
	Sphere_3 m_boundingSphere ;
	int m_k ;
	FT m_gaussianH ;
	mutable Cell_handle m_hint; // last cell found = hint for next search
	boost::shared_ptr< SplatsAABBTree > m_splatsTree ;
	int m_rayTries ;
	FT m_ransacThres ;
	FT m_maxVal ;
	FT m_confidenceThres ;
	FT m_uncertainBand ;
	FT m_smoothPrior ;
	FT m_areaConst ;
	bool m_correct ;
	FT m_stWeight ;
	FT m_functionErrorBound ;
	FT m_radiusEdgeBound ;
	FT m_dataWeight ;
	bool m_boundedSurface ;
	Sphere_3 m_enlargedBoundingSphere ;
	Plane_3 m_sphericalCapPlane ;
	FT m_timeImplicit ;
	FT m_timeSigning ;
	FT m_bandRadius ;
	int m_bandKnn ;
	FT m_epsilon ;
	bool m_signPointsOutsideBand ;



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
	unsigned int delaunay_refinement(FT radius_edge_ratio_bound, ///< radius edge ratio bound (ignored if zero)
									FT cell_radius_bound, ///< cell radius bound (ignored if zero)
									unsigned int max_vertices, ///< number of vertices bound
									FT enlarge_ratio) ///< bounding box enlarge ratio
	{
		// CGAL_TRACE("Calls delaunay_refinement(radius_edge_ratio_bound=%lf, cell_radius_bound=%lf, max_vertices=%u, enlarge_ratio=%lf)\n",
					//radius_edge_ratio_bound, cell_radius_bound, max_vertices, enlarge_ratio);

		EvaluateExactImplicitFunction funcExact( *this ) ; 
		m_enlargedBoundingSphere = enlarged_bounding_sphere(enlarge_ratio);
		//unsigned int nb_vertices_added = unsigned_distance_refine_triangulation( *m_tr,
		//																		 funcExact, 
		//																		 radius_edge_ratio_bound,
		//																		 cell_radius_bound,
		//																		 m_functionErrorBound,
		//																		 max_vertices,
		//																		 0,
		//																		 m_enlargedBoundingSphere ) ;
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

	FT random_in( const double a,
                  const double b)
	{
		double r = rand() / (double)RAND_MAX;
		return (FT)(a + (b - a) * r);
	}

	Vector_3 random_vector()
	{
		/* Using random_in function (appears to have some problem returning the same number for several iterations...) */
		// FT x = random_in(0.0,1.0);
		// FT y = random_in(0.0,1.0);
		// FT z = random_in(0.0,1.0);
		
		/* Using CGAL's random number generator */
		FT x = CGAL::default_random.get_double() ;
		FT y = CGAL::default_random.get_double() ;
		FT z = CGAL::default_random.get_double() ;
		
		return Vector_3( x, y, z ) ;
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
		
		return pow( ( val / steps ), m_smoothPrior ) + m_areaConst ;
	}



	/* Sphere intersection */
	// Modified functions extracted from "Sphere_oracle_3.h"

	/* Derived from boost::tuple< int, FT, FT >
					intersection_line_sphere_lambda( const Surface_3& sphere,
													 const Point& a,
													 const Point& b ) const */
	boost::tuple< int, FT, FT >
	intersection_line_sphere( const Sphere_3& sphere,
							  const Point_3& a,
							  const Point_3& b ) const
	{
		/*
			Let the vectorial line equation:
			m = a + lambda * ( b - a )
			(a, b, and m are points, and lambda if a real.)
          
			Let c be the center of the sphere, of radius r.
			The intersection of the line and the sphere is given by:
			(c-m)^2 = r^2
			That is:
						((c-a)^2 - r^2)
			- 2 lambda (c-a)*(b-a)
			+ lambda^2 (b-a)^2 == 0

			(second degre equation)

			deltaprime = delta/4 = ((c-a)(b-a))^2 - (b-a)^2 * ( (c-a)^2 -r^2 )

			if delta > 0, root_1 = ((c-a)(b-a) - \sqrt(delta/4)) / (b-a)^2
						root_2 = ((c-a)(b-a) + \sqrt(delta/4)) / (b-a)^2
					(root_1 < root_2)
		*/

		typedef typename K::Vector_3 Vector_3;

		typename K::Construct_vector_3 vector =
			K().construct_vector_3_object();
		typename K::Compute_scalar_product_3 scalar_product = 
			K().compute_scalar_product_3_object();
		typename K::Compute_squared_distance_3 squared_distance = 
			K().compute_squared_distance_3_object();
		typename K::Construct_center_3 center =
			K().construct_center_3_object();
		typename K::Compute_squared_radius_3 squared_radius =
			K().compute_squared_radius_3_object();

		const Point_3 c = center(sphere);
		const Vector_3 ab = vector(a, b);
		const Vector_3 ac = vector(a, c);
		const FT ab_ac = scalar_product(ab, ac);
		const FT ab2 = squared_distance(a, b);
		const FT ac2 = squared_distance(a, c);
		const FT r2 = squared_radius(sphere);
		const FT deltaprime = ab_ac * ab_ac - ab2 * ( ac2 - r2 );

		switch( CGAL::sign(deltaprime) )
		{
		case CGAL::ZERO:
			return boost::make_tuple(1, ab_ac / ab2, 0);
		case CGAL::POSITIVE:
		{
			const FT sqrt_deltaprime = CGAL::sqrt(deltaprime);
			return boost::make_tuple(2,
						(ab_ac - sqrt_deltaprime) / ab2,
						(ab_ac + sqrt_deltaprime) / ab2);
		}
		case CGAL::NEGATIVE:
			break;
		}

		return boost::make_tuple(0, 0, 0);
	} //end intersection_line_sphere_lambda



	/** 
		Derived from clip_ray( const Surface_3& sphere,
							   const Ray_3& r,
							   Point_3& a,
							   Point_3& b) 
		In "CGAL/Sphere_oracle_3.h" 
	*/
	bool ray_sphere_intersection( const Sphere_3& sphere,
								  const Ray_3& r,
								  Point_3& int1,
								  Point_3& int2) const
	{
		typedef typename K::Vector_3 Vector;
		typename K::Construct_point_on_3 point_on =
		K().construct_point_on_3_object();
		typename K::Construct_vector_3 vector =
		K().construct_vector_3_object();
		typename K::Construct_scaled_vector_3 scaled_vector = 
		K().construct_scaled_vector_3_object();
		typename K::Construct_translated_point_3 translated_point = 
		K().construct_translated_point_3_object();

		Point_3 a = point_on( r, 0 ) ;
		Point_3 b = point_on( r, 1 ) ;

		int number_of_roots ;
		FT root_1, root_2 ;
        
		boost::tie( number_of_roots, root_1, root_2 ) = intersection_line_sphere( sphere, a, b ) ;

		if( number_of_roots == 2 && root_2 > FT(0) )
		{
			const Vector ab = vector( a, b ) ;
			int1 = translated_point( a, scaled_vector( ab, root_2 ) ) ;
			if( root_1 > FT(0))
				int2 = translated_point( a, scaled_vector( ab, root_1 ) ) ;
			// if root_1 <= 0, a is in the ball
			return true;
		}
		// else r does not intersect the sphere
		return false;
	} // end clip_ray



	bool spherical_cap_intersection( const Sphere_3& sphere,
									 const Ray_3& r,
									 const Plane_3& p ) const
	{

		Point_3 int1, int2 ;
		ray_sphere_intersection( sphere, r, int1, int2 ) ;

		// std::cout << "Ray = " << r << std::endl ;
		// std::cout << "Ray-Sphere intersection 1 = " << int1 << std::endl ;
		// std::cout << "Ray-Sphere intersection 2 = " << int2 << std::endl ;

		return ( p.oriented_side( int1 ) == CGAL::ON_POSITIVE_SIDE ) ;
	}



	// ---> Computations of weights on the graph <---
	
	// S-T weights computation
	void computeSTWeights( GraphType& g ) {


		CGAL::Timer task_timer; task_timer.start() ;
		std::cout << "  - Sign confidence..." << std::flush ;
		
		// Run over each vertex
		Finite_vertices_iterator v, e;
		std::vector< int > odds ;
		std::vector< int > evens ;
		int i = 0 ;
		int numVert = m_tr->number_of_vertices() ;
		int verbosePercent = 0 ;
		int verboseLastPercent = 0 ;	
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
		
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

			// Just apply this part to the points outside the uncertainty band...
			// if ( (*this)( v->point() ) >= m_uncertainBand ) {
			if ( this->invalidValue( v->f() ) ) {

				// Throw a series of random rays and compute parity
				int even = 0 ;
				int odd = 0 ;
				for ( int i = 0; i < m_rayTries; i++ ) {
					std::vector< Point_3 > intersections ;
					// std::cout << "- Test rays:" << std::endl ;
					// do {
						// Random ray originating from the vertex
						Ray_3 ray( v->point(), random_vector() ) ;
						// std::cout << ray << std::endl ;

						// Compute intersections using RANSAC
						intersections = getIntersectionsRANSAC( ray ) ;
					// } while ( intersections.size() > 0 ) ; // TEST: Do not trust empty intersections (they are the majority, and may biase the results)

					int numIntersect = intersections.size() ;

					if ( m_boundedSurface &&
						 spherical_cap_intersection( m_enlargedBoundingSphere, 
													 ray,
													 m_sphericalCapPlane ) )
						numIntersect++ ;

					if ( numIntersect % 2 == 0 ) 
						even++ ;
					else 
						odd++ ;
				}

				// Simple sign
				//if ( even < odd ) {
				//	// std::cout << "Sign changed!" << std::endl ;
				//	v->f() = -v->f() ;
				//}

				odds.push_back( odd ) ;
				evens.push_back( even ) ;

			}
			else {
				odds.push_back( 0 ) ;
				evens.push_back( 0 ) ;
			}

		}
		CGAL_TRACE_STREAM << "\b\b\b\b100%, " << task_timer.time() << " seconds, "
						  << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
						  << std::endl ;
		

		// --- Debug (Start) ---		
		i = 0 ;
		std::ofstream ofi( "./_OUT/InsideVertices_Init.xyz", std::ios_base::out ) ;
		std::ofstream ofo( "./_OUT/OutsideVertices_Init.xyz", std::ios_base::out ) ;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++ )
		{
			if ( odds[i] != 0 | evens[i] != 0 ) {
				if ( odds[i] > evens[i] )
					ofi << v->point() << std::endl ;
				else
					ofo << v->point() << std::endl ;
			}
		}
		ofi.close() ;
		ofo.close() ;
		// --- Debug (End) ---

		task_timer.reset() ;
		std::cout << "  - Sign optimization: " << std::endl ;
		std::cout << "    - Building graph (S-T weights)..." ;

		// --- Debug (Start) ---
		std::ofstream ofis( "./_OUT/InsideVertices_SOURCE.xyz", std::ios_base::out ) ;
		std::ofstream ofos( "./_OUT/OutsideVertices_SINK.xyz", std::ios_base::out ) ;
		// --- Debug  (End)  ---

		
		i = 0 ;
		verbosePercent = 0 ;
		verboseLastPercent = 0 ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
		
		// ---> Normal in/out confidence computation <---
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++)
		{
		
			// Show status on screen
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( numVert ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				CGAL_TRACE_STREAM << "\b\b\b\b" << std::setfill(' ') << std::setw(3) << verbosePercent << "%" ;
				verboseLastPercent = verbosePercent ;
			}
			
			//if ( v->f() >= m_uncertainBand ) {
			if ( this->invalidValue( v->f() ) ) {
				// FT confidence = 2.0* max( evens[i], odds[i] ) / m_rayTries - 1 ;

				// FT confidenceE = 2.0 * static_cast<double>( evens[i] ) / m_rayTries - 1 ;
				// FT confidenceO = 2.0 * static_cast<double>( odds[i] ) / m_rayTries - 1 ;

				FT confidenceE = static_cast<double>( evens[i] ) / m_rayTries ;
				FT confidenceO = static_cast<double>( odds[i] ) / m_rayTries ;

				FT confidence = std::max<double>( confidenceE, confidenceO ) ;
			
				/*std::cout << "Evens = " << evens[i] << std::endl ;
				std::cout << "Odds = " << odds[i] << std::endl ;
				std::cout << "ConfidenceE = " << confidenceE << std::endl ;
				std::cout << "ConfidenceO = " << confidenceO << std::endl ;
				std::cout << "Confidence = " << confidence << std::endl ;*/
				// g -> add_tweights( i,   /* capacities */ confidence*evens[i], confidence*odds[i] );*/

				if ( confidence > m_confidenceThres ) {
					// Just add links to terminals for highly confident in/out
					// g -> add_tweights( i,   /* capacities */ confidence*evens[i], confidence*odds[i] );
					// g -> add_tweights( i,   /* capacities */ confidence*( evens[i] / m_rayTries )*m_stWeight, confidence*( odds[i] / m_rayTries )*m_stWeight ) ;
					g.add_tweights(	v->ind(), 
										confidenceE*m_stWeight, /* capacities for Outside label */ 
										confidenceO*m_stWeight  /* capacities for Inside label */ ) ;

					// --- Debug (Start) ---
					if ( odds[i] > evens[i] ) {
						ofis << v->point() << std::endl ;
					}
					else {
						ofos << v->point() << std::endl ;
					}
					// --- Debug  (End)  ---
				}
				else {
					// Ambiguous...
					g.add_tweights( v->ind(),   /* capacities */ 0.5, 0.5 ) ;
				}
			}
			else {
				// Ambiguous...
				g.add_tweights( v->ind(),   /* capacities */ 0.5, 0.5 ) ;
			}
		}

		//// --- Debug (Start) ---
		ofis.close() ;
		ofos.close() ;
		// --- Debug  (End)  ---
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << ")" << std::endl ;
	}



	void computeSTWeightsRestricted( GraphType& g ) {

		// --- Debug (Start) ---
		std::ofstream oi( "./_OUT/PointsInInterphase.xyz", std::ios_base::out ) ;
		std::ofstream o_in( "./_OUT/PointsInInterphase_IN.xyz", std::ios_base::out ) ;
		std::ofstream o_out( "./_OUT/PointsInInterphase_OUT.xyz", std::ios_base::out ) ;
		// --- Debug  (End)  ---

		std::cout << "    - Building graph (s-t weights)..." ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << 0 << "%" ;
		
		CGAL::Timer task_timer; task_timer.start() ;

		Finite_vertices_iterator v, e ;
		int i = 0, verbosePercent = 0, verboseLastPercent = 0 ;
		int numVert = m_tr->number_of_vertices() ;
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

			if ( this->invalidValue( v->f() ) ) {
				// The vertex is outside the band, check that at least one adjacent vertex is inside...
				std::vector< Vertex_handle > adjVert ;
				m_tr->finite_adjacent_vertices( v, std::back_inserter( adjVert ) ) ;
				typename std::vector< Vertex_handle >::iterator itav = adjVert.begin() ;
				bool someInBand = false ;
				while ( !someInBand && itav != adjVert.end() ) {
					someInBand = ( !this->invalidValue( (*itav)->f() ) ) ;
					++itav ;
				}
				if ( someInBand ) {
					// Vertex is on the interphase between the band and the rest of in/out volume
					// Compute confidence for being in/out using the random ray-splats intersection counting mechanism
					FT confidenceIn , confidenceOut ;
					stochasticInOutLabelling( v->point(), confidenceIn , confidenceOut ) ;
				
					g.add_tweights(	v->ind(),
									confidenceOut*m_stWeight, /* capacities for Outside label */ 
									confidenceIn*m_stWeight   /* capacities for Inside label */ ) ;

					// --- Debug (Start) ---
					oi << v->point() << std::endl ;
					if ( confidenceIn > confidenceOut ) 
						o_in << v->point() << std::endl ;
					else
						o_out << v->point() << std::endl ;
					// --- Debug  (End)  ---
				}
				/*else {
					g.add_tweights(	v->ind(), 0.5, 0.5 ) ;
				}*/
			}
			/*else {
				g.add_tweights(	v->ind(), 0.5, 0.5 ) ;
			}*/
		}
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << ")" << std::endl ;

		// --- Debug (Start) ---
		oi.close() ;
		// --- Debug  (End)  ---
	}



	// Smooth weights computation
	void computeSmoothWeights( GraphType& g ) 
	{
		std::cout << "    - Building graph (edge weights)..." ;
		CGAL_TRACE_STREAM << std::setfill(' ') << std::setw(3) << 0 << "%" ;
		
		CGAL::Timer task_timer; task_timer.start() ;

		int verbosePercent = 0 ;
		int verboseLastPercent = 0 ;
		Finite_edges_iterator it;
		int numEdges = m_tr->number_of_edges() ;
		task_timer.reset() ;
		int i = 0 ;
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

			if ( !this->invalidValue( v0->f() ) && !this->invalidValue( v1->f() ) ) {

				// Original
				//FT ew = pow( ( ( v0->f() + v1->f() ) / 2 ), m_smoothPrior ) + m_areaConst ;
				//g->add_edge( v0->ind(), v1->ind(), ( ew / m_maxVal )*m_dataWeight, ( ew / m_maxVal )*m_dataWeight ) ;
				//std::cout << "Weight = " << ( ew / m_maxVal )*m_dataWeight << std::endl ;

				// N-Cut like
				FT val = ( v0->f() + v1->f() ) / 2 ;
				// FT w = this->weightInverseWendland( val, m_maxVal ) + m_areaConst ;
				// FT w = pow( this->weightInverseWendland( val, m_maxVal ), m_smoothPrior ) + m_areaConst ;
				FT w = pow( val/m_maxVal, m_smoothPrior ) + m_areaConst ;
				g.add_edge( v0->ind(), v1->ind(), w*m_dataWeight, w*m_dataWeight ) ;
				// std::cout << "Weight = " << w << std::endl ;
				
				///*std::cout << "val = " << val << std::endl ;
				//std::cout << "w = " << w << std::endl ;
				//FT ew = pow( ( ( v0->f() + v1->f() ) / 2 ), m_smoothPrior ) + m_areaConst ;
				//std::cout << "ew = " << ew << std::endl ;*/
			}
			else {
				
				if ( this->invalidValue( v0->f() ) && this->invalidValue( v1->f() ) ) {
					// Use NO WEIGHT when using non-restricted s-t labelling!
					g.add_edge( v0->ind(), v1->ind(), 1, 1 ) ;
				}
				else {
					FT ew = 0 ;
					if ( this->invalidValue( v0->f() ) ) {
						// Original
						//ew = pow( v1->f(), m_smoothPrior ) + m_areaConst ;
						// N-Cut like
						// ew = this->weightInverseWendland( v1->f(), m_maxVal ) + m_areaConst ;
						// ew = pow( this->weightInverseWendland( v1->f(), m_maxVal ), m_smoothPrior ) + m_areaConst ;
						ew = pow( v1->f()/m_maxVal, m_smoothPrior ) + m_areaConst ;
						// std::cout << "v1->f() = " << v1->f() << std::endl ;
					}
					else {
						// Original
						// ew = pow( v0->f(), m_smoothPrior ) + m_areaConst ;
						// N-Cut like
						// ew = this->weightInverseWendland( v0->f(), m_maxVal ) + m_areaConst ;
						// ew = pow( this->weightInverseWendland( v0->f(), m_maxVal ), m_smoothPrior ) + m_areaConst ;
						ew = pow( v0->f()/m_maxVal, m_smoothPrior ) + m_areaConst ;
						// std::cout << "v0->f() = " << v0->f() << std::endl ;
					}
					//std::cout << "Weight(!) = " << ( ew / m_maxVal )*m_dataWeight << std::endl ;
					// std::cout << "Weight(!) = " << ew << std::endl ;

					// Original
					// g->add_edge( v0->ind(), v1->ind(), ( ew / m_maxVal )*m_dataWeight, ( ew / m_maxVal )*m_dataWeight ) ;
					// N-Cut like
					g.add_edge( v0->ind(), v1->ind(), ew*m_dataWeight, ew*m_dataWeight ) ;
				}
				
			}
		}
		std::cout << "\b\b\b\b100%, done (" << task_timer.time() << ")" << std::endl ;
	}


	void stochasticInOutLabelling( const Point_3& p, FT& confidenceIn, FT& confidenceOut ) 
	{
		// Throw a series of random rays and compute parity
		int even = 0 ;
		int odd = 0 ;
		for ( int i = 0; i < m_rayTries; i++ ) {
			std::vector< Point_3 > intersections ;
			// std::cout << "- Test rays:" << std::endl ;
			// do {
				// Random ray originating from the vertex
				Ray_3 ray( p, random_vector() ) ;
				// std::cout << ray << std::endl ;

				// Compute intersections using RANSAC
				intersections = getIntersectionsRANSAC( ray ) ;
			// } while ( intersections.size() > 0 ) ; // TEST: Do not trust empty intersections (they are the majority, and may biase the results)

			int numIntersect = intersections.size() ;

			if ( m_boundedSurface &&
					spherical_cap_intersection( m_enlargedBoundingSphere, 
												ray,
												m_sphericalCapPlane ) ) {
				numIntersect++ ;
				// std::cout << "Intersects the spherical cap" << std::endl ;
			}

			if ( numIntersect % 2 == 0 ) 
				even++ ;
			else 
				odd++ ;
		}

		confidenceIn = static_cast<double>( even ) / m_rayTries ;
		confidenceOut = static_cast<double>( odd ) / m_rayTries ;
	}


	FT secondaryBandKnn( const Point_3& p ) const
	{
		if ( m_bandKnn <= 0 ) 
			return UNDEFINED_VALUE ;

		if ( m_bandRadius > m_gaussianH ) {

			// Add the points on a secondary band with a default value
			Neighbor_search search( *m_pTree, p, 1 ) ;
				
			Splat_3 curSplat = boost::get<1>( search.begin()->first ) ;
			FT dist = CGAL::sqrt( CGAL::squared_distance( curSplat.center(), p ) ) ;

			// Inside second band
			if ( dist < m_bandRadius ) {
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
					return sumWeightedDist / sumWeights ;
				else {
					return UNDEFINED_VALUE ;
				}				
			}
			else {
				return UNDEFINED_VALUE ;
			}
		}
		else {
			return UNDEFINED_VALUE ;
		}
	}

	
	bool invalidValue( const FT& val ) const {
		// const FT epsilon = 0.001 ;
		// std::cout << "std::abs( val - 999999.0 ) = " << std::abs( val - 999999.0 ) << std::endl ;
		return std::abs( val - UNDEFINED_VALUE ) < m_epsilon ;
	}

} ;

#endif // SPLATSUNSIGNEDDISTANCEFUNCTION_H
