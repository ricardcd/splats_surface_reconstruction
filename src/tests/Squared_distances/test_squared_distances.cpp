// CGAL includes & redefinitions
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include "SquaredDistance/squared_distance_segment_triangle_3.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_3		Point_3 ;
typedef K::Triangle_3	Triangle_3 ;
typedef K::Segment_3	Segment_3 ;
typedef K::Line_3	Line_3 ;
typedef K::Vector_3	Vector_3 ;

// Std
#include <iostream>
using namespace std ;


int main( int argc, char **argv ) {
	
	//if ( argc < 16 ) {
	//	// Default triangle-segment
	//}

	Point_3 pt1( 0.8147, 0.9134, 0.2785 ) ;
    Point_3 pt2( 0.9058, 0.6324, 0.5469 ) ;
    Point_3 pt3( 0.1270, 0.0975, 0.9575 ) ;
	Triangle_3 tri( pt1, pt2, pt3 ) ;

	// Does not intersect
	/*Point_3 ps1( 1.0, 1.0, 0.0 ) ;
	Point_3 ps2( 0.0, 0.0, 1.0 ) ;*/

	Point_3 ps1( 1.0, 1.0, 0.0 ) ;
	Point_3 ps2( 0.6, 0.5, 0.2 ) ;

	// Intersects
	// Point_3 ps1( 1.0, 1.0, 0.0 ) ;
	// Point_3 ps2( 0.5, 0.5, 1.0 ) ;

    Segment_3 seg( ps1, ps2 ) ;
    Line_3 l(ps1, ps2);
    Vector_3 v = l.to_vector();

	cout << "Triangle = " << tri << endl ;
	cout << "Segment = " << seg << endl ;
//	cout << "Squared distance = " << CGAL::squared_distance( seg, tri ) << endl ;

	return EXIT_SUCCESS ;
}