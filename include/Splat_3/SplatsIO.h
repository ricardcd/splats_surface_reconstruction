#ifndef SPLATSIO_H
#define SPLATSIO_H

#include <iostream>
#include <fstream>
#include "Splat_3.h"

namespace SplatsIO {

template< class K >
bool WriteSplats( std::ostream& stream, ///< output stream.
				  const std::vector< typename CGAL::Splat_3<K> >& splats ) {

	typedef typename K::FT		FT ;
	typedef typename K::Point_3	Point_3 ;
	typedef typename K::Vector_3	Vector_3 ;

	typename std::vector< typename CGAL::Splat_3<K> >::const_iterator it ;
	for ( it = splats.begin(); it != splats.end(); ++it ) {
		Point_3 o = it->monge().origin() ;
		Vector_3 d1 = it->monge().maximal_principal_direction() ;
		Vector_3 d2 = it->monge().minimal_principal_direction() ;
		Vector_3 normal = it->monge().normal_direction() ;
		FT rad = CGAL::sqrt( it->squared_radius() ) ;

		stream << o << " " << d1 << " " << d2 << " " << normal << " " << rad ;

		// Coefficients
		std::vector<FT> coefs = it->monge().coefficients() ;
		stream << " " << coefs.size() ;
		for ( int i = 0; i < coefs.size(); i++ ) {
			stream << " " << coefs[i] ;
		}
		stream << std::endl ;
	}

	return !stream.fail() ;

}


template< class K >
bool ReadSplats(  std::istream& stream, ///< input stream.
				  std::vector< typename CGAL::Splat_3<K> >& splats ) {

	typedef typename K::FT			FT ;
	typedef typename K::Point_3		Point_3 ;
	typedef typename K::Vector_3	Vector_3 ;

	while( stream.good() ) {
		double mo_x, mo_y, mo_z, md1_x, md1_y, md1_z, md2_x, md2_y, md2_z, mn_x, mn_y, mn_z, rad ;
		int	numCoefs ;
		stream	>> mo_x 
				>> mo_y 
				>> mo_z 
				>> md1_x 
				>> md1_y 
				>> md1_z 
				>> md2_x 
				>> md2_y 
				>> md2_z 
				>> mn_x 
				>> mn_y 
				>> mn_z 
				>> rad 
				>> numCoefs ;
		
		typename std::vector< FT > coefs ;
		for( int i = 0; i < numCoefs; i++ ) {
			FT coef ;
			stream >> coef ;
			coefs.push_back( coef ) ;
		}

		Point_3 o( mo_x, mo_y, mo_z ) ;
		Vector_3 d1( md1_x, md1_y, md1_z ) ;
		Vector_3 d2( md2_x, md2_y, md2_z ) ;
		Vector_3 n( mn_x, mn_y, mn_z ) ;
		
		typename CGAL::Splat_3<K>::Monge_form monge( o, d1, d2, n, coefs ) ;
		typename CGAL::Splat_3<K> splat( monge, rad*rad ) ;

		splats.push_back( splat ) ;

	}

	return !stream.eof() ;

}

} // End namespace SplatsIO

#endif // SPLATSIO_H