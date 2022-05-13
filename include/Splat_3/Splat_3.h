#ifndef SPLAT_3_H
#define SPLAT_3_H

#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
#include <CGAL/Bbox_3.h>

namespace CGAL {

template < class Kernel > 
class Splat_3 {

	typedef typename Kernel::RT										RT ;
	typedef typename Kernel::FT										FT ;
	typedef typename Kernel::Point_3									Point_3 ;
	typedef typename Kernel::Vector_3									Vector_3 ;
	typedef typename Kernel::Sphere_3									Sphere_3 ;
	typedef typename CGAL::Monge_via_jet_fitting< Kernel, Kernel >	Monge_via_jet_fitting ;
	
public:  
	typedef typename Monge_via_jet_fitting::Monge_form				Monge_form ;

	typedef typename CGAL::Splat_3<Kernel> Rep ;

	const Rep& rep() const
	{
		return *this;
	}

	Rep& rep()
	{
		return *this;
	}


	/* Constructors */
	Splat_3( ) {
		m_monge = Monge_form() ;
		m_boundSphere = Sphere_3() ;
	}


	Splat_3( Monge_form monge, const double &sqRadius ) {
		m_monge = monge ;
		m_boundSphere = Sphere_3( monge.origin(), sqRadius ) ;
	}

	Splat_3( Monge_form monge, const Sphere_3& boundSphere ) {
		m_monge( monge() ) ;
		m_boundSphere = boundSphere ; 
	}
	
	Splat_3( const CGAL::Splat_3<Kernel> &s ) {
		m_monge = s.monge() ;
		m_boundSphere = Sphere_3( s.boundSphere() ) ;
	}
	
	/* Methods */
	void set_monge( Monge_form monge ) {
		m_monge = monge ;
	}

	FT squared_radius() const {
		return m_boundSphere.squared_radius() ;
	}

	Point_3 center() const {
		return m_monge.origin() ;
	}

	Monge_form monge() const {
		return m_monge ;
	}

	Sphere_3 boundSphere() const { 
		return m_boundSphere ;
	}

	Point_3 project( const Point_3 &p ) const {
		Point_3 res = m_monge.project( p ) ;

		return res ;
	}
  
	FT algebraicDistance( const Point_3 &p ) const {
		return m_monge.algebraicDistance( p ) ;
	}

	Bbox_3 bbox() const
	{
		return typename Kernel::Construct_bbox_3()( m_boundSphere ) ; 
	}

	bool has_monge() const {
		// If one of the vectors of the coordinate system of the monges is null, then the monge is not constructed properly
		//return m_monge.maximal_principal_direction() != Vector_3( 0.0, 0.0, 0.0  ) ; 
		/*return ( m_monge.maximal_principal_direction() != CGAL::NULL_VECTOR && 
				m_monge.minimal_principal_direction() != CGAL::NULL_VECTOR && 
				m_monge.normal_direction() != CGAL::NULL_VECTOR ) ; */
		return ( m_monge.normal_direction() != CGAL::NULL_VECTOR ) ;
	}

private:
	
	Monge_form m_monge ;  
	Sphere_3 m_boundSphere ;

} ;

}

#endif // SPLAT_3_H