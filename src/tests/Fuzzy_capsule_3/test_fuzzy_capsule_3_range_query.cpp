//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include "SplatsCreation/KDTreeSearch/Fuzzy_capsule_3.h"
#include "SplatsCreation/KDTreeSearch/My_Search_traits_3.h"
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/IO/read_xyz_points.h>
#include <fstream>

typedef CGAL::Cartesian<double>							K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel		K ;

typedef K::Point_3											Point_3;
typedef K::Segment_3										Segment_3;
typedef CGAL::My_Search_traits_3<K>							Traits;
typedef CGAL::Kd_tree<Traits>								Tree;
typedef CGAL::Fuzzy_sphere<Traits>							Fuzzy_sphere;
typedef CGAL::Fuzzy_iso_box<Traits>							Fuzzy_iso_box;
typedef CGAL::Fuzzy_capsule_3<Traits>						Fuzzy_capsule_3;

int main( int argc, char **argv ) {
  
	// Read the dataset
	std::ifstream in( argv[1] ) ; 
	std::vector< Point_3 > points ;
	if ( !in ) {
		std::cerr << "Cannot open file!" << std::endl ;
	}
	read_xyz_points( in, std::back_inserter( points ) ) ;
	
	Tree tree( points.begin(), points.end() ) ;
	
	// define range query objects
	Point_3 p( atof( argv[2] ), atof( argv[3] ), atof( argv[4] ) ) ;
	Point_3 q( atof( argv[5] ), atof( argv[6] ), atof( argv[7] ) ) ;

	Segment_3 seg( p, q ) ;

	std::cout << "Query segment:" << std::endl ;
	std::cout << seg << std::endl ;

	double rad = atof( argv[8] ) ;	
	std::cout << "Query radius:" << rad << std::endl ;

	Fuzzy_capsule_3 cap( seg, atof( argv[8] ) ) ;
  
	std::ofstream ofc( "insideCapsule.xyz" ) ;
	std::vector< Point_3 > fuzzyQueryResult ;
	tree.search( std::back_inserter( fuzzyQueryResult ), cap ) ;

	std::cout << "Max Radius sq = " << rad*rad << std::endl ;
	std::cout << "Distances:" << std::endl ;
	for ( std::vector< Point_3 >::iterator it = fuzzyQueryResult.begin(); it != fuzzyQueryResult.end(); ++it ) {
		std::cout << CGAL::squared_distance( *it, seg )  << std::endl ;
	}


	std::ostream_iterator<Point_3> out_it(ofc, "\n") ;
	std::copy( fuzzyQueryResult.begin(), fuzzyQueryResult.end(), out_it ) ;
	ofc.close() ;


	std::cout << "Brute force search in progress..." << rad << std::endl ;
	std::ofstream ofcbf( "insideCapsuleBruteForce.xyz" ) ;
	std::vector< Point_3 >::iterator it ;
	double sqRadius = rad * rad ;
	for ( it = points.begin(); it != points.end(); ++it ) {
		double sqDist = CGAL::squared_distance( *it, seg ) ;
		if ( sqDist < sqRadius )
			ofcbf << *it << std::endl ;
	}
	ofcbf.close() ;




	return 0;
}