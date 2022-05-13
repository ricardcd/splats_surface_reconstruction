// Std includes
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

// CGAL includes & redefinitions
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	// AABB_tree
#include "Splat_3/Splat_3.h"
#include "Splat_3/SplatsIO.h"
#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
	// General
#include <CGAL/Random.h>
	// Timer (debugging)
#include <CGAL/Timer.h>	
	// Point set processing
#include <CGAL/mst_orient_normals.h>

// CImg include (only used here for easily dealing with input parameters)
#include "CImg.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT								FT ;
typedef K::Triangle_3						Triangle_3 ;
typedef CGAL::Splat_3<K>					Splat_3 ;
typedef K::Point_3							Point_3 ;
typedef K::Vector_3							Vector_3 ;

typedef CGAL::Monge_via_jet_fitting< K, K >			Monge_via_jet_fitting ;
typedef Monge_via_jet_fitting::Monge_form			Monge_form ;

// Point with normal vector stored in a std::pair.
typedef std::pair< Point_3, Vector_3 > OrientedPoint ;

// Other definitions
using namespace std ;



/* Main function */
int main ( int argc, char **argv) {

	// Parse/default input parameters
	const char*		inputSplatsFile			= cimg_option( "-i",		(char*)0,		"Input splats file" ) ;	
	const char*		outputSplatsFile		= cimg_option( "-o",		(char*)0,		"Output oriented splats file" ) ;
	const int		k						= cimg_option( "-k",		15,				"K-nearest neighbors to take into account during MST orientation." ) ;
	
	// Get the file type (*.splat)
	string inputFilePathStr( inputSplatsFile ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	string inputFileExtension = inputFilePathStr.substr( indexFileExtension+1, ( inputFilePathStr.size()-(indexFileExtension+1) ) ) ;	
	if ( strcmp( inputFileExtension.c_str(), "splat" ) != 0 ) {
		cerr << "Unrecognized file type..." << endl ;
		return -1 ;
	}

	CGAL::Timer timer ; // Reusable timer initialization

	// Loading features
	
	// Read the input pointset file	
	cout << "Reading splats from file..." << flush ;	
	timer.start() ;
	// Opening file
	ifstream file( inputSplatsFile, ifstream::in ) ;
	if ( !file ) {
		cerr << "Cannot open input file!" << endl ;	
		return -1 ;
	}

	std::vector< Splat_3 > splats ;
	bool readed = SplatsIO::ReadSplats( file, splats ) ;
	/*if ( !readed ) {
		cerr << "Error while reading the input file!" << endl ;	
		return -1 ;
	}*/
	timer.stop() ;
	cout << "done (" << timer.time() << " s), readed " << splats.size() << " splats." << endl ;
	timer.reset() ;
	
	/* Orient the splats */

	// Create a point-normal set
	std::vector< Splat_3 >::iterator it ;
	std::vector< OrientedPoint > orientedPoints ;
	for ( it = splats.begin(); it != splats.end(); ++it ) {
		OrientedPoint op( it->center(), it->monge().normal_direction() ) ;
		orientedPoints.push_back( op ) ;
	}

	// Orient normals coherently
	std::vector< OrientedPoint >::iterator unoriented_points_begin =	CGAL::mst_orient_normals(	orientedPoints.begin(), orientedPoints.end(),
																									CGAL::First_of_pair_property_map< OrientedPoint >(),
																									CGAL::Second_of_pair_property_map< OrientedPoint >(),
																									k ) ;

	std::vector< OrientedPoint >::iterator itOp ;
	for ( it = splats.begin(), itOp = orientedPoints.begin(); it != splats.end(); ++it, ++itOp ) {
		if ( itOp->second != it->monge().normal_direction() ) {
			// std::cout << "Normal changed!" << it->monge().normal_direction() << "-->" ;
			// it->monge().normal_direction() = -it->monge().normal_direction() ;
			Monge_form m = it->monge() ;
			m.normal_direction() = -m.normal_direction() ;
			m.minimal_principal_direction() = -m.minimal_principal_direction() ;
			std::vector<FT> coefs = m.coefficients() ;
			for ( int i = 0 ; i<coefs.size(); i++ ) {
				coefs[i] = -coefs[i] ;
			}
			m.coefficients() = coefs ;
			it->set_monge( m ) ;
			// std::cout << it->monge().normal_direction() << std::endl ;
		}
	}

	// Save results
	cout << "- Writing results to file..." << flush ;
	ofstream outFile ;
	outFile.open( outputSplatsFile ) ;
	if ( !outFile.is_open() ) {
		cerr << "Cannot open the output file for writing!" << endl ;
		return -1 ;
	}
	SplatsIO::WriteSplats<K>( outFile, splats ) ;
	cout << "done" << std::endl ;

	return 0 ;
}
