#include <CGAL/Cartesian.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include "SplatsCreation/KDTreeSearch/Fuzzy_capsule_3.h"
#include <CGAL/Fuzzy_iso_box.h>
#include "SplatsCreation/KDTreeSearch/My_Search_traits_3.h"
#include <fstream>

#include <CGAL/Orthogonal_k_neighbor_search.h>

typedef CGAL::Cartesian<double>							K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel		K ;

typedef K::Point_3											Point_3;
typedef K::Segment_3										Segment_3;
typedef K::Point_3											Point_d;
typedef CGAL::My_Search_traits_3<K>							Traits;
typedef CGAL::Orthogonal_k_neighbor_search< Traits >	Neighbor_search ;
typedef CGAL::Random_points_in_cube_3<Point_d>				Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>		N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits>								Tree;
typedef CGAL::Fuzzy_sphere<Traits>							Fuzzy_sphere;
typedef CGAL::Fuzzy_iso_box<Traits>							Fuzzy_iso_box;
typedef CGAL::Fuzzy_capsule_3<Traits>						Fuzzy_capsule_3;

int main() {

    const int N = 100000;
    // generator for random data points in the square ( (-1000,-1000), (1000,1000) )
    Random_points_iterator rpit( 1000.0 ) ;

    // Insert N points in the tree
    Tree tree(N_Random_points_iterator(rpit,0),
              N_Random_points_iterator(rpit,N));

    // define range query objects
    double  pcoord[3] = { 300, 300, 300 };
    double  qcoord[3] = { 900.0, 900.0, 900.0 };
    Point_3 p(pcoord[0], pcoord[1], pcoord[2]);
    Point_3 q(qcoord[0], qcoord[1], qcoord[2]);
    Segment_3 seg( p, q ) ;
    Fuzzy_sphere fs( p, 700.0, 100.0 ) ;
    Fuzzy_iso_box fib( p, q, 100.0 ) ;
    Fuzzy_capsule_3 cap( seg, 100.0 ) ;

    std::cout << "SPHERE QUERY" << std::endl ;
    std::cout << "------------" << std::endl ;
    std::cout << "points approximately in fuzzy range query" << std::endl;
    std::cout << "with center (300.0, 300.0, 300.0)" << std::endl;
    std::cout << "and fuzzy radius <200.0,400.0> are:" << std::endl;
    std::ofstream ofs( "outSphere.xyz" ) ;
    tree.search(std::ostream_iterator<Point_d>(ofs, "\n"), fs);
    ofs.close() ;

    std::cout << "ISOBOX QUERY" << std::endl ;
    std::cout << "------------" << std::endl ;
    std::cout << "Computing points approximately in fuzzy range query (will be stored in \"outIsoBox.xyz\" file)";
    std::cout << "[<200,4000>,<800,1000>]]^3" << std::endl;
    std::ofstream ofib( "outIsoBox.xyz" ) ;
    tree.search(std::ostream_iterator<Point_d>(ofib, "\n"), fib);
    ofib.close() ;

    std::cout << "CAPSULE QUERY" << std::endl ;
    std::cout << "-------------" << std::endl ;
    std::cout << "points approximately in fuzzy range query " << std::endl;
    std::ofstream ofc( "outCapsule.xyz" ) ;
    tree.search(std::ostream_iterator<Point_d>(ofc, "\n"), cap);
    ofc.close() ;

//  std::cout << "CAPSULE QUERY, 5 Nearest neighbor" << std::endl ;
//  std::cout << "-------------" << std::endl ;
//  Neighbor_search search( tree, seg, 5 ) ;
//  // Find k+1 nearest neighbors on the structure (k+1 because first NN will be the point itself)
//  for( Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
//	  std::cout << it->first << std::endl ;
//  }

    return 0;
}