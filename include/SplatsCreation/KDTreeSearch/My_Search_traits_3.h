#ifndef CGAL_MY_SEARCH_TRAITS_3_H
#define CGAL_MY_SEARCH_TRAITS_3_H

namespace CGAL {


  template <class K>

  class My_Search_traits_3 {

  public:
    
    typedef typename K::Cartesian_const_iterator_3 Cartesian_const_iterator_d;
    typedef typename K::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_d;
    typedef typename K::Point_3 Point_d;
    typedef typename K::Iso_cuboid_3 Iso_box_d;
    typedef typename K::Sphere_3 Sphere_d;
	typedef typename K::Triangle_3 Triangle_3;
	typedef typename K::Point_3 Point_3;
	typedef typename K::Segment_3 Segment_3;
    typedef typename K::Line_3 Line_3;
    typedef typename K::Vector_3 Vector_3;
    typedef typename K::Construct_iso_cuboid_3 Construct_iso_box_d;

    typedef typename K::Construct_min_vertex_3 Construct_min_vertex_d;
    typedef typename K::Construct_max_vertex_3 Construct_max_vertex_d;
    typedef typename K::Construct_center_3 Construct_center_d;
    typedef typename K::Compute_squared_radius_3 Compute_squared_radius_d;
    typedef typename K::FT FT;
  
    Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
       return Construct_cartesian_const_iterator_d();
    }
  
  };

  
} // namespace CGAL
#endif // MY_SEARCH_TRAITS_3_H