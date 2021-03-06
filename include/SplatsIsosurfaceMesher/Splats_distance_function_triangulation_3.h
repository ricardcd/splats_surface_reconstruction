#ifndef CGAL_SPLATS_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
#define CGAL_SPLATS_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H

#include <CGAL/property_map.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/centroid.h>

#include <vector>

namespace CGAL {


/// The Splats_distance_function_vertex_base_3 class is the default
/// vertex class of the Reconstruction_triangulation_3 class.
///
/// It provides the interface requested by the Poisson_reconstruction_function class:
/// - Each vertex stores a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
///
/// @heading Parameters:
/// @param Gt   Geometric traits class
/// @param Cb   Vertex base class, model of TriangulationVertexBase_3.

template < typename Gt,
           typename Vb = Triangulation_vertex_base_3<Gt> >
class Splats_distance_function_vertex_base_3 : public Vb
{
// Public types
public:

  /// Geometric traits class
  typedef Gt Geom_traits;

  // Repeat Triangulation_vertex_base_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef typename Vb::Cell_handle Cell_handle;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other                       Vb2;
    typedef Splats_distance_function_vertex_base_3<Geom_traits, Vb2> Other;
  };
  /// @endcond

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;

// data members
private:

  // TODO: reduce memory footprint
  FT m_f; // value of the implicit function // float precise enough?  
  unsigned char m_type; // INPUT or STEINER
  int m_ind; // Linear index on the triangulation (just affects finite vertices)

// Public methods
public:

	// WARNING: test-> change initialization of the function value m_f to -1.0 instead of 0.0

  Splats_distance_function_vertex_base_3()
    : Vb(), m_f(FT(-1.0)), m_type(0), m_ind( -1 )
  {}

  Splats_distance_function_vertex_base_3(const Point& p)
    : Vb(p), m_f(FT(-1.0)), m_type(0), m_ind( -1 )
  {}

  Splats_distance_function_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p,c), m_f(FT(-1.0)), m_type(0), m_ind( -1 )
  {}

  Splats_distance_function_vertex_base_3(Cell_handle c)
    : Vb(c), m_f(FT(-1.0)), m_type(0), m_ind( -1 )
  {}

  /// Gets/sets the value of the implicit function.
  /// Default value is 0.0.
  FT  f() const { return m_f; }
  FT& f()       { return m_f; }

  /// Gets/sets the type = INPUT or STEINER.
  unsigned char  type() const { return m_type; }
  unsigned char& type()       { return m_type; }

  /// Gets/sets the linear index of the vertex
  int  ind() const { return m_ind; }
  int& ind()       { return m_ind; }
  
	
// Private methods
private:

	/// Copy constructor and operator =() are not implemented.
	Splats_distance_function_vertex_base_3(const Splats_distance_function_vertex_base_3& toCopy) ;
    Splats_distance_function_vertex_base_3& operator =(const Splats_distance_function_vertex_base_3& toCopy);

}; // end of Splats_distance_function_vertex_base_3


///// Helper class:
///// Splats_distance_function_triangulation_default_geom_traits_3
///// changes in a geometric traits class the Point_3 type to
///// Point_with_normal_3<BaseGt>.
/////
///// @heading Parameters:
///// @param BaseGt   Geometric traits class.
//template <class BaseGt>
//struct Splats_distance_function_triangulation_default_geom_traits_3 : public BaseGt
//{
//  typedef Point_with_normal_3<BaseGt> Point_3;
//};


/// The Splats_distance_function_triangulation_3 class
/// provides the interface requested by the Poisson_Splats_distance_function_function class:
/// - Each vertex stores a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
/// The vertex class must derive from Splats_distance_function_vertex_base_3.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Splats_distance_function_vertex_base_3.

template <class Gt,          
          class Tds_ = Triangulation_data_structure_3<Splats_distance_function_vertex_base_3<Gt> > >
class Splats_distance_function_triangulation_3 : public Delaunay_triangulation_3<Gt,Tds_>
{
// Private types
private:

  // Base class
  typedef Delaunay_triangulation_3<Gt,Tds_>  Base;

  // Auxiliary class to build an iterator over input points.
  class Is_steiner_point
  {
  public:
      typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

      bool operator()(const Finite_vertices_iterator& v) const
      {
        return (v->type() == Splats_distance_function_triangulation_3::STEINER);
      }
  };

// Public types
public:

  /// Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
  typedef Gt  Geom_traits;

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  typedef Tds_ Triangulation_data_structure;
  typedef typename Base::Segment      Segment;
  typedef typename Base::Triangle     Triangle;
  typedef typename Base::Tetrahedron  Tetrahedron;
  typedef typename Base::Line         Line;
  typedef typename Base::Ray          Ray;
  typedef typename Base::Object       Object;
  typedef typename Base::Cell_handle   Cell_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell   Cell;
  typedef typename Base::Vertex Vertex;
  typedef typename Base::Facet  Facet;
  typedef typename Base::Edge   Edge;
  typedef typename Base::Cell_circulator  Cell_circulator;
  typedef typename Base::Facet_circulator Facet_circulator;
  typedef typename Base::Cell_iterator    Cell_iterator;
  typedef typename Base::Facet_iterator   Facet_iterator;
  typedef typename Base::Edge_iterator    Edge_iterator;
  typedef typename Base::Vertex_iterator  Vertex_iterator;
  typedef typename Base::Point_iterator Point_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Base::All_cells_iterator       All_cells_iterator;
  typedef typename Base::All_vertices_iterator       All_vertices_iterator;
  typedef typename Base::Locate_type Locate_type;
  /// @endcond

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector; 
  typedef typename Geom_traits::Point_3 Point;  
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// Point type
  enum Point_type {
    INPUT,    ///< Input point.
    STEINER   ///< Steiner point created by Delaunay refinement.
  };

  /// Iterator over input vertices.
  typedef Filter_iterator<Finite_vertices_iterator, Is_steiner_point>
                                                    Input_vertices_iterator;

  /// Iterator over input points.
  typedef Iterator_project<Input_vertices_iterator,
                           Project_point<Vertex> >  Input_point_iterator;

// Public methods
public:

  /// Default constructor.
  Splats_distance_function_triangulation_3()
  {
  }

  // Default copy constructor and operator =() are fine.

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  using Base::points_begin;
  using Base::points_end;
  using Base::number_of_vertices;
  using Base::finite_vertices_begin;
  using Base::finite_vertices_end;
  using Base::all_vertices_begin;
  using Base::all_vertices_end;

  using Base::geom_traits;
  /// @endcond

  /// Gets first iterator over input vertices.
  Input_vertices_iterator input_vertices_begin() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point(),
                                     finite_vertices_begin());
  }
  /// Gets past-the-end iterator over input vertices.
  Input_vertices_iterator input_vertices_end() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point());
  }

  /// Gets iterator over the first input point.
  Input_point_iterator input_points_begin() const
  {
      return Input_point_iterator(input_vertices_begin());
  }
  /// Gets past-the-end iterator over the input points.
  Input_point_iterator input_points_end() const
  {
      return Input_point_iterator(input_vertices_end());
  }

  /// Gets the bounding sphere of all points.
  Sphere bounding_sphere() const
  {
    typedef Min_sphere_of_spheres_d_traits_3<Gt,FT> Traits;
    typedef Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere Traits_sphere;

    // Represents *all* points by a set of spheres with 0 radius
    std::vector<Traits_sphere> spheres;
    spheres.reserve(number_of_vertices());

    for (Point_iterator it=points_begin(), eit=points_end();
         it != eit; ++it)
      spheres.push_back(Traits_sphere(*it,0));

    // Computes min sphere
    Min_sphere ms(spheres.begin(),spheres.end());
    typename Min_sphere::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    FT cx = *coord++;
    FT cy = *coord++;
    FT cz = *coord++;
    return Sphere(Point(cx,cy,cz), ms.radius()*ms.radius());
  }

  /// Gets the bounding sphere of input points.
  Sphere input_points_bounding_sphere() const
  {
    typedef Min_sphere_of_spheres_d_traits_3<Gt,FT> Traits;
    typedef Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere Traits_sphere;

    // Represents *input* points by a set of spheres with 0 radius
    std::vector<Traits_sphere> spheres;
    for (Input_point_iterator it=input_points_begin(), eit=input_points_end();
         it != eit;
         ++it)
      spheres.push_back(Traits_sphere(*it,0));

    // Computes min sphere
    Min_sphere ms(spheres.begin(),spheres.end());
    typename Min_sphere::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    FT cx = *coord++;
    FT cy = *coord++;
    FT cz = *coord++;
    return Sphere(Point(cx,cy,cz), ms.radius()*ms.radius());
  }

  /// Insert point in the triangulation.
  /// Default type is INPUT.
  Vertex_handle insert(const Point& p,
                       Point_type type = INPUT,
                       Cell_handle start = Cell_handle())
  {
    Vertex_handle v = Base::insert(p, start);
    v->type() = type;
    return v;
  }

  /// Insert the [first, beyond) range of points in the triangulation using a spatial sort.
  /// Default type is INPUT.
  ///
  /// @commentheading Template Parameters:
  /// @param InputIterator iterator over input points.
  /// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Point_3.
  ///        It can be omitted if InputIterator value_type is convertible to Point_3.
  /// @param NormalPMap is a model of boost::ReadablePropertyMap with a value_type = Vector_3.
  ///
  /// @return the number of inserted points.

  // This variant requires all parameters.
  template <typename InputIterator>
  int insert(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator over the input points.
    Point_type type = INPUT)
  {
    int n = number_of_vertices();
	
    std::vector<Point> points;
    for( InputIterator it = first; it != beyond; ++it )
    {
        points.push_back( *it ) ;
    }
	
    // Spatial sorting
    std::random_shuffle( points.begin(), points.end() ) ;
    spatial_sort (points.begin(), points.end(), geom_traits());
	
    // Insert in triangulation
    Cell_handle hint;
    for (typename std::vector<Point>::const_iterator p = points.begin();
         p != points.end(); ++p)
    {
      Vertex_handle v = insert(*p, type, hint);
      hint = v->cell();
    }

    return number_of_vertices() - n;
  }

  /// Delaunay refinement callback:
  /// insert STEINER point in the triangulation.
  template <class CellIt>
  Vertex_handle
  insert_in_hole(const Point& p, CellIt cell_begin, CellIt cell_end,
	         Cell_handle begin, int i,
                 Point_type type = STEINER)
  {
      Vertex_handle v = Base::insert_in_hole(p, cell_begin, cell_end, begin, i);
      v->type() = type;
      return v;
  }

  /// Index unconstrained vertices following the order of Finite_vertices_iterator.
  /// @return the number of unconstrained vertices.
  unsigned int index_unconstrained_vertices()
  {
		// Why?
	  return 0 ;
  }

}; // end of Splats_distance_function_triangulation_3


} //namespace CGAL

#endif // CGAL_SPLATS_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
