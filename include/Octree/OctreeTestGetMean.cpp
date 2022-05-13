#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include "CGAL_Octree.h"

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;
typedef K::FT FT;
typedef K::Point_3 Point_3;



int main( int argc, char *argv[] )
{
    if ( argc < 4 ) {
        std::cerr << "Not enough input arguments..." << std::endl ;
        std::cerr << "Usage:" << std::endl ;
        std::cerr << "    OctreeTestGetMean <in_points_xyz_file> <octree_depth> <out_mean_points_xyz>" << std::endl ;
        return -1 ;
    }

    // Reads a .xyz point set file in points.
    std::cout << "Loading data..." << std::flush ;
    std::vector< Point_3 > points ;
    std::ifstream stream( argv[1] ) ;
    if ( !stream ||
         !CGAL::read_xyz_points( stream, std::back_inserter( points ) ) )
    {
        std::cerr << "Error: cannot read input xyz file!" << std::endl;
        return -1 ;
    }
    std::cout << "done." << std::endl ;

    // Create the octree
    CGAL_Octree::CGAL_Octree< K > octree ;
    CGAL_Octree::CubicBounds< K > bounds( points ) ;

    std::cout << "Building tree..." << std::flush ;
    octree.build( points, atoi(argv[2]), bounds ) ;
    std::cout << "done." << std::endl ;

    std::cout << "Computing mean points..." << std::flush ;
    std::vector< Point_3 > meanPts ;
    octree.getMeanPoints( meanPts ) ;
    std::cout << "done (" << meanPts.size() << " pts)." << std::endl ;

    std::cout << "Saving results..." << std::flush ;
    std::ofstream ofs( argv[3] ) ;
    CGAL::write_xyz_points( ofs, meanPts.begin(), meanPts.end() ) ;
    std::cout << "done." << std::endl ;

    return 0 ;
}
