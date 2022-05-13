//
// Based on: CGAL/Unsigned_distance_in_band_mesh_cell_criteria_3.h
//
// Author: Ricard Campos

#ifndef CGAL_UNSIGNED_DISTANCE_IN_BAND_MESH_CRITERIA_3_H
#define CGAL_UNSIGNED_DISTANCE_IN_BAND_MESH_CRITERIA_3_H

#include <iostream>
#include <functional>

namespace CGAL {

template <	class Tr,
			class Function >
class Unsigned_distance_in_band_mesh_cell_criteria_3 
{
	// Attributes
	double squared_radius_bound_ ;
	double radius_edge_bound_ ;
	double func_band_bound_ ;
	Function func_ ;

public:

	typedef typename Tr::Geom_traits::FT		FT ;
	typedef typename Tr::Geom_traits::Point_3	Point_3 ;
	
	struct Cell_quality : public std::pair<double, double>
	{
		typedef std::pair<double, double> Base;

		Cell_quality() : Base() {}
		Cell_quality(double _aspect, double _sq_size) : Base(_aspect, _sq_size) {};

		double sq_size() const { return second; }
		double aspect() const { return first; }

		// q1<q2 means q1 is prioritised over q2
		// ( q1 == *this, q2 == q )
		bool operator<(const Cell_quality& q) const
		{
			// std::cout << sq_size() << " - " << q.sq_size() << std::endl ;

			if( sq_size() > 1 )
				if( q.sq_size() > 1 )
					return ( sq_size() > q.sq_size() );
				else
					return true; // *this is big but not q
			else
				if( q.sq_size() >  1 )
					return false; // q is big but not *this
			return( aspect() > q.aspect() );
		}
	};



	inline
	double squared_radius_bound() const 
	{
		return squared_radius_bound_; 
	}

	typedef typename Tr::Cell_handle Cell_handle;

	Unsigned_distance_in_band_mesh_cell_criteria_3( const Function& func, //< Function
											const double radius_edge_bound = 2, //< radius edge ratio bound (ignored if zero)
											const double radius_bound = 0, //< cell radius bound (ignored if zero)
											const double func_approx_bound = 0.001 ) //< function approximation error in the linear approximation provided by the triangulation
	: func_(func),
	  squared_radius_bound_(radius_bound*radius_bound),
	  radius_edge_bound_(radius_edge_bound),
	  func_band_bound_(func_approx_bound)
	{}

	inline 
	void set_squared_radius_bound(const double squared_radius_bound) 
	{ 
		squared_radius_bound_ = squared_radius_bound;
	}

	inline
	double radius_edge_bound() const 
	{
		return radius_edge_bound_; 
	}

	inline 
	void set_radius_edge_bound(const double radius_edge_bound) 
	{
		radius_edge_bound_ = radius_edge_bound;
	}

	inline
	double func_approx_bound() const 
	{
		return func_band_bound_; 
	}

	inline 
	void set_func_approx_bound(const double func_approx_bound) 
	{
		func_band_bound_ = func_approx_bound;
	}

	FT operator()(Point_3 p) const
	{
		return func(p);
	}



	class Is_bad
	{
	protected:
		const double radius_edge_bound_ ;
		const double squared_radius_bound_ ;
		const double func_band_bound_ ;
		Function func_ ;

	public:
		typedef typename Tr::Point Point_3;

		Is_bad( const Function func,
				const double radius_edge_bound, 
				const double squared_radius_bound, 
				const double func_approx_bound )
		: func_(func),
		  radius_edge_bound_(radius_edge_bound),
		  squared_radius_bound_(squared_radius_bound),
		  func_band_bound_(func_approx_bound)
		{}

		bool operator()( const Cell_handle& c,
						 Cell_quality& qual) const
		{
			// Check that the cell is in the band!
			const Point_3& p = c->vertex(0)->point();
			/*if ( func_( p ) > func_band_bound_ ) 
				return false ;*/
			// std::cout << func_( p ) << std::endl ;
			const Point_3& q = c->vertex(1)->point();
			/*if ( func_( q ) > func_band_bound_ ) 
				return false ;*/
			// std::cout << func_( q ) << std::endl ;
			const Point_3& r = c->vertex(2)->point();
			/*if ( func_( r ) > func_band_bound_ ) 
				return false ;*/
			// std::cout << func_( r ) << std::endl ;
			const Point_3& s = c->vertex(3)->point();
			/*if ( func_( s ) > func_band_bound_ ) 
				return false ;*/
			// std::cout << func_( s ) << std::endl ;
			if ( func_( p ) > func_band_bound_ && func_( q ) > func_band_bound_ && func_( r ) > func_band_bound_ && func_( s ) > func_band_bound_ ) 
				return false ;
			
			
			typedef typename Tr::Geom_traits Geom_traits;
			typedef typename Geom_traits::Compute_squared_radius_3 Radius;
			typedef typename Geom_traits::Compute_squared_distance_3 Distance;
			typedef typename Geom_traits::FT FT;

			Radius radius = Geom_traits().compute_squared_radius_3_object();
			Distance distance = Geom_traits().compute_squared_distance_3_object();

			double size = to_double( radius(p, q, r, s) ) ;

			if( squared_radius_bound_ != 0 )
			{
				qual.second = size / squared_radius_bound_;
				// normalized by size bound to deal
				// with size field
				if( qual.sq_size() > 1 )
				{
					// std::cout << "Cell-radius is bad" << std::endl ;
					qual.first = 1; // (do not compute aspect)
					return true;
				}
			}
			if( radius_edge_bound_ == 0 )
			{
				qual = Cell_quality(0,1);
				return false;
			}

			double min_sq_length = CGAL::to_double(distance(p, q));
			min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, r)));
			min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, s)));
			min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, r)));
			min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, s)));
			min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(r, s)));

			qual.first = size / min_sq_length;

			// std::cout << "Radius-edge is bad" << std::endl ;
			return (qual.first > radius_edge_bound_);
		}


		//void barycentric_coordinates(	const Point_3& p,
		//							Cell_handle cell,
		//							FT& a,
		//							FT& b,
		//							FT& c,
		//							FT& d) const
		//{
		//	const Point_3& pa = cell->vertex(0)->point();
		//	const Point_3& pb = cell->vertex(1)->point();
		//	const Point_3& pc = cell->vertex(2)->point();
		//	const Point_3& pd = cell->vertex(3)->point();

		//	FT v = volume(pa,pb,pc,pd);
		//	a = CGAL::abs(volume(pb,pc,pd,p) / v);
		//	b = CGAL::abs(volume(pa,pc,pd,p) / v);
		//	c = CGAL::abs(volume(pb,pa,pd,p) / v);
		//	d = CGAL::abs(volume(pb,pc,pa,p) / v);
		//}

		///// Generates a uniformly distributed random point inside the tetrahedron described by [ p0, p1, p2, p3 ]
		///// Also returns the baricentric coordinates of the point inside the cell
		///// (adapted from: http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html)
		//Point_3 random_point_in_tetrahedron( Point_3 p0, 
		//									 Point_3 p1, 
		//									 Point_3 p2, 
		//									 Point_3 p3,
		//									 FT& a,
		//									 FT& s,
		//									 FT& t,
		//									 FT& u) const
		//{
		//	s = CGAL::default_random.get_double() ;
		//	t = CGAL::default_random.get_double() ;
		//	u = CGAL::default_random.get_double() ;

		//	if ( s+t > 1.0 ) { // cut'n fold the cube into a prism
		//		s = 1.0 - s ;
		//		t = 1.0 - t ;
		//	}
		//	if ( t+u>1.0 ) { // cut'n fold the prism into a tetrahedron
		//		double tmp = u ;
		//		u = 1.0 - s - t ;
		//		t = 1.0 - tmp ;
		//	} else if( s+t+u>1.0 ) {
		//		double tmp = u ;
		//		u = s + t + u - 1.0 ;
		//		s = 1 - t - tmp ;
		//	}
		//	a = 1-s-t-u ; // [ a, s, t, u ] are the barycentric coordinates of the random point.

		//	Point_3 p0a( p0.x()*a, p0.y()*a, p0.z()*a ) ;
		//	Point_3 p1s( p1.x()*s, p1.y()*s, p1.z()*s ) ;
		//	Point_3 p2t( p2.x()*t, p2.y()*t, p2.z()*t ) ;
		//	Point_3 p3u( p3.x()*u, p3.y()*u, p3.z()*u ) ;

		//	return Point_3( p0a.x() + p1s.x() + p2t.x() + p3u.x(), 
		//					p0a.y() + p1s.y() + p2t.y() + p3u.y(), 
		//					p0a.z() + p1s.z() + p2t.z() + p3u.z() ) ;
		//}

		//FT linear_approx( FT fp, FT fq, FT fr, FT fs, FT ba, FT bb, FT bc, FT bd ) const {

		//	// Compute the function value at the vertices of the tetrahedron
		//	// If one of them is a negative value, default to -1 (i.e., mark it as invalid)
		//	if ( fp == 999999 ) {
		//		return 999999 ;
		//	}
		//	if ( fq == 999999 ) {
		//		return 999999 ;
		//	}
		//	if ( fr == 999999 ) {
		//		return 999999 ;
		//	}
		//	if ( fs == 999999 ) {
		//		return 999999 ;
		//	}
		//	
		//	return ba*fp + bb*fq + bc*fr + bd*fs ;

		//}

	
	}; // end Is_bad
	
	Is_bad is_bad_object() const
	{ 
		return Is_bad( func_, radius_edge_bound_, squared_radius_bound_, func_band_bound_ ) ; 
	}

}; // end Unsigned_distance_in_band_mesh_cell_criteria_3



template <typename Tr, typename Func>
std::ostream& operator<< ( std::ostream& os,
							const typename Unsigned_distance_in_band_mesh_cell_criteria_3<Tr, Func>::Cell_quality& q)
{
	return os << q.sq_size() << ", " << q.aspect() << ", " << q.func_approx_bound() ;
}

} // end namespace CGAL

#endif // CGAL_UNSIGNED_DISTANCE_IN_BAND_MESH_CRITERIA_3_H
