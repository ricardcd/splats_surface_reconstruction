//
// Based on: CGAL/poisson_refine_triangulation.h
//
// Author: Ricard Campos

#ifndef CGAL_UNSIGNED_DISTANCE_REFINE_TRIANGULATION_H
#define CGAL_UNSIGNED_DISTANCE_REFINE_TRIANGULATION_H

// CGAL
#include <CGAL/trace.h>
#include <CGAL/Mesher_level.h>
#include "Unsigned_distance_refine_cells_3.h"
#include "Unsigned_distance_mesh_cell_criteria_3.h"
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>

namespace CGAL {


/// Utility class for unsigned_distance_refine_triangulation():
/// implements Delaunay refinement in a loose bounding
/// box of point set (break bad tetrahedra, where
/// bad means badly shaped or too big).
///
/// This class must be derived to inherit from Mesher_level.
template <class Tr,
          class Criteria,
          class Surface,
          class Oracle,
          class Container = Meshes::Double_map_container<typename Tr::Cell_handle,
                                                         typename Criteria::Cell_quality>
>
class Unsigned_distance_mesher_level_impl_base :
    public Mesh_3::Unsigned_distance_refine_tets_with_oracle_base<Tr, Criteria, Surface, Oracle, Container>
{
  typedef Mesh_3::Unsigned_distance_refine_tets_with_oracle_base<Tr, Criteria, Surface, Oracle, Container> Base;

public:
  // Inherited methods and fields used below
  using Base::triangulation_ref_impl;
  using Base::oracle;
  using Base::surface;
  using Base::should_be_refined;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Point Point;
  typedef typename Base::Cell_quality Cell_quality;

public:
  /** \name CONSTRUCTORS */

  Unsigned_distance_mesher_level_impl_base(Tr& t, Criteria crit, unsigned int max_vert, Surface& surface, Oracle& oracle)
    : Base(t, crit, surface, oracle),
      max_vertices(max_vert) ///< number of vertices bound (ignored if zero)

  {}

protected:
  /* --- protected functions --- */

  bool test_if_cell_is_bad(const Cell_handle c)
  {
    Cell_quality q;
    if( is_in_domain(c) && should_be_refined(c, q) )
    {
      this->add_bad_element(c, q);
      return true;
    }
    return false;
  }

  bool is_in_domain(const Cell_handle c)
  {
    return oracle.is_in_volume(surface, triangulation_ref_impl().dual(c));
  }

public:
  /* Overriden functions of this level: */
  /* we override all methods that call test_if_cell_is_bad() */

  void scan_triangulation_impl()
  {
    for(typename Tr::Finite_cells_iterator cit = triangulation_ref_impl().finite_cells_begin(),
        eit = triangulation_ref_impl().finite_cells_end();
        cit != eit;
        ++cit)
      test_if_cell_is_bad(cit);
  }

  void after_insertion_impl(const Vertex_handle& v)
  {
    update_star(v);
  }

  void update_star(const Vertex_handle& v)
  {
    // scan refiner
    typedef std::vector<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;
    incident_cells.reserve(32);
    triangulation_ref_impl().incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator cit = incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
    {
      if( ! triangulation_ref_impl().is_infinite(*cit) )
        test_if_cell_is_bad(*cit);
    }
  }

  /// Tells if the algorithm is done.
  bool no_longer_element_to_refine_impl() const
  {
    return Base::no_longer_element_to_refine_impl() ||
           (max_vertices > 0 && triangulation_ref_impl().number_of_vertices() >= max_vertices);
  }


private:
  /* --- private datas --- */
  unsigned int max_vertices; ///< number of vertices bound (ignored if zero)

}; // end Unsigned_distance_mesher_level_impl_base


/// Utility class for unsigned_distance_refine_triangulation():
/// glue class that inherits from both Mesher_level
/// and Unsigned_distance_mesher_level_impl_base.
template <typename Tr,
          typename Criteria,
          typename Surface,
          typename Oracle = typename CGAL::Surface_mesh_traits_generator_3<Surface>::type,
          typename PreviousLevel = Null_mesher_level
 >
class Unsigned_distance_mesher_level :
  public Unsigned_distance_mesher_level_impl_base<Tr, Criteria, Surface, Oracle>,
  public Mesher_level <
    Tr,
    Unsigned_distance_mesher_level<Tr, Criteria, Surface, Oracle, PreviousLevel>,
    typename Tr::Cell_handle,
    PreviousLevel,
    Triangulation_mesher_level_traits_3<Tr>
  >
{
  typedef Unsigned_distance_mesher_level_impl_base<Tr, Criteria, Surface, Oracle> Base;

public:
  typedef Mesher_level <
    Tr,
    Unsigned_distance_mesher_level<Tr, Criteria, Surface, Oracle, PreviousLevel>,
    typename Tr::Cell_handle,
    PreviousLevel,
    Triangulation_mesher_level_traits_3<Tr>
  > Mesher;

  Unsigned_distance_mesher_level(Tr& t, Criteria criteria, unsigned int max_vertices, Surface& surface, Oracle& oracle, PreviousLevel& previous_level)
    : Base(t, criteria, max_vertices, surface, oracle),
      Mesher(previous_level)
  {
  }

}; // end class Unsigned_distance_mesher_level


/// Delaunay refinement in a loose bounding box
/// of input point set (break bad tetrahedra, where
/// bad means badly shaped or too big).
/// @return the number of vertices inserted.
///
/// @commentheading Preconditions:
/// - Tr must use a geometric traits with robust circumcenter computation.
/// - convergence is guaranteed if radius_edge_ratio_bound >= 1.0.
///
/// @commentheading Template Parameters:
/// @param Tr 3D Delaunay triangulation.
/// @param Surface Sphere_3 or Iso_cuboid_3.
template <	typename Tr,
			typename Surface>
unsigned int unsigned_distance_refine_triangulation( Tr& tr,
													 double radius_edge_ratio_bound, ///< radius edge ratio bound (>= 1.0)
													 double cell_radius_bound, ///< cell radius bound (ignored if zero)
													 unsigned int max_vertices, ///< number of vertices bound (ignored if zero)
													 Surface& enlarged_bbox) ///< new bounding sphere or box
{
	// CGAL_TRACE("Calls unsigned_distance_refine_triangulation()\n");

	// Convergence is guaranteed if radius_edge_ratio_bound >= 1.0
	CGAL_surface_reconstruction_points_precondition( radius_edge_ratio_bound >= 1.0 ) ;

	// Basic geometric types
	typedef typename Tr::Geom_traits Gt;
	typedef typename Gt::FT FT;
	typedef typename Gt::Point_3 Point;

	// Mesher_level types
	typedef Unsigned_distance_mesh_cell_criteria_3< Tr > Tets_criteria;
	typedef typename CGAL::Surface_mesh_traits_generator_3<Surface>::type Oracle;
	typedef Unsigned_distance_mesher_level<Tr, Tets_criteria, Surface, Oracle, Null_mesher_level> Refiner;

	// CGAL_TRACE("  Creates queue\n");

	int nb_vertices = tr.number_of_vertices(); // get former #vertices

	// Delaunay refinement
	Tets_criteria tets_criteria( radius_edge_ratio_bound*radius_edge_ratio_bound, 
								 cell_radius_bound ) ;

	Oracle oracle;
	Null_mesher_level null_mesher_level;
	Refiner refiner(tr, tets_criteria, max_vertices, enlarged_bbox, oracle, null_mesher_level);
	
	// std::cout << "    - Scanning triangulation..." ;
	refiner.scan_triangulation(); // Push bad cells to the queue
	// std::cout << "done." << std::endl ;
	
	refiner.refine(Null_mesh_visitor()); // Refine triangulation until queue is empty

	int nb_vertices_added = tr.number_of_vertices() - nb_vertices;

	// CGAL_TRACE("End of unsigned_distance_refine_triangulation()\n");

	return (unsigned int) nb_vertices_added ;
}

} //namespace CGAL

#endif // CGAL_UNSIGNED_DISTANCE_REFINE_TRIANGULATION_H
