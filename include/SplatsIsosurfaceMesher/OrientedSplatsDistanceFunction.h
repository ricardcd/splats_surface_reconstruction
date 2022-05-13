#ifndef ORIENTEDSPLATSUNSIGNEDDISTANCEFUNCTION_H
#define ORIENTEDSPLATSUNSIGNEDDISTANCEFUNCTION_H

/* Includes */
// CGAL
#include <CGAL/trace.h>
	// Spatial searching
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
//#include "KDTreeSearch/My_Search_Traits_3.h"
//#include "KDTreeSearch/Fuzzy_capsule_3.h"
	// Triangulation (refinement)
#include <SplatsIsosurfaceMesher/Splats_distance_function_triangulation_3.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/poisson_refine_triangulation.h>
#include <CGAL/Robust_circumcenter_filtered_traits_3.h>
	// Other
#include <CGAL/compute_average_spacing.h>
	// Related to this project
#include "Splat_3/Splat_3.h"
// Boost
#include <boost/shared_ptr.hpp>
#include <boost/iterator/zip_iterator.hpp>
// Graph Cuts
#include <graph.h>
// Std
#include <algorithm>
// Undefined value define
#ifndef UNDEFINED_VALUE
	#define UNDEFINED_VALUE 999999.0
#endif



template < class K >
class OrientedSplatsDistanceFunction {

public:
	/* Typedefs */
	typedef typename K::FT												FT ;
	typedef typename K::Point_3											Point_3 ;
	typedef typename K::Sphere_3										Sphere_3 ;
	typedef typename K::Ray_3											Ray_3 ;
	typedef typename K::Vector_3										Vector_3 ;
	typedef typename CGAL::Splat_3< K >									Splat_3 ;
	typedef typename K::Segment_3										Segment_3 ;
		
	typedef typename CGAL::Monge_via_jet_fitting<K, K>					Monge_via_jet_fitting ;
	typedef typename Monge_via_jet_fitting::Monge_form					Monge_form ;

	/// Internal 3D triangulation, of type Reconstruction_triangulation_3.
	// Note: poisson_refine_triangulation() requires a robust circumcenter computation.
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

	/* Nearest neighbors search definitions */
	// Property map, each point takes as 
	typedef boost::tuple< Point_3, Splat_3 >							Point_and_splat;

	class My_point_property_map{
	public:
		typedef Point_3 value_type ;
		typedef const value_type& reference ;
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
	


	/* Functions */

	OrientedSplatsDistanceFunction( const std::vector< Splat_3 >& splats,
							const int& k = 25,
							const FT& gaussianRadiusFactor = 0.2,
							const double epsilon = 0.001 ) 
							: m_tr(new Triangulation) {
		m_k = k ;
		m_epsilon = epsilon ;
		
		FT factor = gaussianRadiusFactor ;
		//if ( gaussianRadiusFactor < 0 || gaussianRadiusFactor > 1 ) {
		//	std::cerr << "[SplatsUnsignedDistanceFunction] WARNING gaussianRadiusFactor must be between 0 and 1! Defaulting to gaussianRadiusFactor=0.2." << std::endl ;
		//	factor = 0.2 ;
		//}
		
		// Extract centers
		typename std::vector< Splat_3 >::const_iterator it ;
		typename std::vector< Point_3 > points ;
		for( it = splats.begin(); it != splats.end(); ++it ) {
			points.push_back( it->center() ) ;
		}
		
		// Build initial triangulation
		m_tr->insert( points.begin(), points.end() ) ;
		
		// Construct the KD-Tree with these center points
		m_pTree = pNNTree( new NNTree( boost::make_zip_iterator(boost::make_tuple( points.begin(), splats.begin() ) ),
									   boost::make_zip_iterator(boost::make_tuple( points.end(), splats.end() ) ) ) ) ;

		
		// Compute gaussian weight H variable
		if ( gaussianRadiusFactor < 1 ) {
			m_gaussianH = CGAL::sqrt( factor * factor * this->bounding_sphere().squared_radius() ) ;
		}
		else {
			// Computes average spacing
			FT average_spacing = CGAL::compute_average_spacing( points.begin(), points.end(), 6 /* knn = 1 ring */ ) ;
			// Compute gaussian weight H variable
			m_gaussianH = gaussianRadiusFactor * average_spacing ;
		}

		m_splats = splats ;

		// Compute the (unsigned) implicit function
		compute_implicit_function() ;

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

	/// Compute the implicit function at the vertices of a refined triangulation (similar to an octree decomposition)
	bool compute_implicit_function() {

		CGAL::Timer task_timer; task_timer.start() ;
		CGAL_TRACE_STREAM << "Delaunay refinement...\n" ;

		// Delaunay refinement
		const FT radius_edge_ratio_bound = 2.5 ; // Modified --> Original = 2.5
		const unsigned int max_vertices = (unsigned int)1e7 ; // max 10M vertices
		const FT enlarge_ratio = 1.5 ;
		const FT radius = sqrt( bounding_sphere().squared_radius() ) ; // get triangulation's radius
		const FT cell_radius_bound = radius/5. ; // large
		unsigned int nb_vertices_added = delaunay_refinement( radius_edge_ratio_bound,
															  cell_radius_bound,
															  max_vertices,
															  enlarge_ratio ) ;

		// Prints status
		CGAL_TRACE_STREAM << "Delaunay refinement: " << "added " << nb_vertices_added << " Steiner points, "
						  << task_timer.time() << " seconds, "
						  << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
						  << std::endl ;
		task_timer.reset() ;

		CGAL_TRACE_STREAM << "Compute implicit function...\n" ;

		// Computes the Implicit function eval()
		// at each vertex of the triangulation.
		Finite_vertices_iterator v, e;
		int i = 0 ;
		m_maxVal = 0 ;
		for( v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end();
			 v != e;
			 ++v, i++)
		{
			v->f() = eval( v->point() ) ;
			v->ind() = i ; // Assign also a linear index

			if ( m_maxVal < v->f() && invalidValue( v->f() ) ) {
				m_maxVal = v->f() ;
			}
		}

		// Prints status
		CGAL_TRACE_STREAM << "Compute implicit function: " << task_timer.time() << " seconds, "
						  << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
						  << std::endl ;
		task_timer.reset() ;

		// Compute the value of the function at each of its vertices
		FT mv = median_value_at_input_vertices() ;
		std::cout << "Median value at input points: " << mv << std::endl ;

		return true;

	}


	/// 'ImplicitFunction' interface: evaluates the implicit function at a given 3D query point.
	FT operator()(const Point_3& p) const {

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
		// std::cout << "f(0) = " << f0 << std::endl ;
		if ( invalidValue( f0 ) ) {
			f0 = m_maxVal ;
		}
		FT f1 = m_hint->vertex(1)->f() ; 
		// std::cout << "f(1) = " << f1 << std::endl ;
		if ( invalidValue( f1 ) ) {
			f1 = m_maxVal ;
		}
		FT f2 = m_hint->vertex(2)->f() ; 
		// std::cout << "f(2) = " << f2 << std::endl ;
		if ( invalidValue( f2 ) ) {
			f2 = m_maxVal ;
		}
		FT f3 = m_hint->vertex(3)->f() ; 
		// std::cout << "f(3) = " << f3 << std::endl ;
		if ( invalidValue( f3 ) ) {
			f3 = m_maxVal ;
		}

		// --- Debug (Start) ---
		//if ( f0 < 0 || f1 < 0 || f2 < 0 || f3 < 0 ) {
		//	std::cout << "Eval = " << a*f0 + b*f1 + c*f2 + d*f3 << std::endl ;
		//}
		// --- Debug  (End)  ---

		return a*f0 + b*f1 + c*f2 + d*f3 ;

	}


	// --> Weighting functions <--

	/// Gaussian weight function
	static FT weightGaussian( const FT& dist, const FT& h = 1.0 ) {
		return exp( -( ( dist*dist ) / ( h*h ) ) ) ;
	}
	
	/// Wendland weight function
	static FT weightWendland( const FT& dist, const FT& h = 1.0 ) {
		return ( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*( 1 - dist/h )*(4*dist / h + 1 ) ;
	}

	// Get the triangulation
	boost::shared_ptr<Triangulation> getTriangulationPtr() {
		return m_tr ;
	}

private:
	/* Attributes */
	boost::shared_ptr<Triangulation> m_tr ; // Main triangulation holding the precomputed implicit function
	pNNTree m_pTree ;
	std::vector< Splat_3 > m_splats ;
	Sphere_3 m_boundingSphere ;
	int m_k ;
	FT m_gaussianH ;
	mutable Cell_handle m_hint; // last cell found = hint for next search
	FT m_maxVal ;
	FT m_epsilon ;

	
	/* Functions */
	
	/// Evaluate the function at a given point
	FT eval( const Point_3 &p ) const {

		FT sumWeightedDist = 0.0 ;
		FT sumWeights = 0.0 ;

		/* K-NN */
		if ( m_k > 0 ) {
			Neighbor_search search( *m_pTree, p, m_k ) ; 
			bool any = false ;
			for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
				// Distance measure
				Splat_3 curSplat = boost::get<1>( it->first ) ;
				FT dist = curSplat.algebraicDistance( p ) ;

				// FT dist2 = curSplat.algebraicDistance( p ) ;
				//Point_3 projectedPoint = curSplat.project( p ) ;
				//FT dist = CGAL::sqrt( CGAL::squared_distance( projectedPoint, p ) ) ;

				//if ( dist2 < 0 )
				//	dist = -dist ;

				//std::cout << dist << "/" << dist2 << std::endl ;

				// Weighting
				if ( abs( dist ) < m_gaussianH ) {
					FT weight = this->weightGaussian( dist, m_gaussianH ) ;

					// std::cout << "weight = " << weight << std::endl ;

					sumWeightedDist += weight * dist ;
					sumWeights += weight ;
					if ( !any )
						any = true ;
				}
				// Else, weight is almost zero, do not take this point into account...
				/*else {
					std::cout << "m_gaussianH = " << m_gaussianH << " / abs( dist ) = " << abs( dist ) << std::endl ;
				}*/
			}

			if ( !any ) {
				// return std::numeric_limits<double>::infinity() ;
				// std::cout << "OUT OF DOMAIN" << std::endl ;				
				return UNDEFINED_VALUE ;
			}
		}
		else {
			/* Radial NN */
			std::vector< Point_and_splat > neighs ;
			Sphere_range_query searchSphere( p, m_gaussianH ) ;
			m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;

			if ( neighs.empty() ) {
				// return std::numeric_limits<double>::infinity() ;
				return UNDEFINED_VALUE ;
			}

			//if ( neighs.size() > 1 )
			//	std::cout << "Num neighs = " << neighs.size() << std::endl ;

			typename std::vector< Point_and_splat >::iterator it ;
			for ( it = neighs.begin(); it!= neighs.end(); ++it ) {
				// Distance measure
				Splat_3 curSplat = boost::get<1>( *it ) ;
				FT dist = curSplat.algebraicDistance( p ) ;

				// Weighting
				FT weight = this->weightGaussian( dist, m_gaussianH ) ;

				sumWeightedDist += weight * dist ;
				sumWeights += weight ;
			}
		}

		////if ( neighs.size() > 1 )
				// std::cout << "Eval = " << sumWeightedDist / sumWeights << std::endl ;

		return sumWeightedDist / sumWeights ;
	}


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
		CGAL_TRACE("Calls delaunay_refinement(radius_edge_ratio_bound=%lf, cell_radius_bound=%lf, max_vertices=%u, enlarge_ratio=%lf)\n",
					radius_edge_ratio_bound, cell_radius_bound, max_vertices, enlarge_ratio);

		Sphere_3 elarged_bsphere = enlarged_bounding_sphere(enlarge_ratio);
		unsigned int nb_vertices_added = poisson_refine_triangulation(*m_tr,radius_edge_ratio_bound,cell_radius_bound,max_vertices,elarged_bsphere);

		CGAL_TRACE("End of delaunay_refinement()\n");

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

	bool invalidValue( const FT& val ) const {
		// const FT epsilon = 0.01 ;
		// std::cout << "std::abs( val - 999999.0 ) = " << std::abs( val - 999999.0 ) << std::endl ;
		return std::abs( val - UNDEFINED_VALUE ) < m_epsilon ;
	}


} ;

#endif // ORIENTEDSPLATSUNSIGNEDDISTANCEFUNCTION_H
