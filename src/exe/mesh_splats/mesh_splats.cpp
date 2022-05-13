// Std includes
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

// CGAL includes & redefinitions
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
	// Surface Mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h> 
	// Mesher
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
	// AABB_tree
#include <CGAL/AABB_tree.h>
#include "Splat_3/AABBTree/AABB_Splat_3_primitive.h"
#include "Splat_3/AABBTree/AABB_Splat_3_intersections.h"
#include "Splat_3/AABBTree/AABB_Splat_3_traits.h"
#include "Splat_3/Splat_3.h"
#include "Splat_3/SplatsIO.h"
//#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
	// General
#include <CGAL/Random.h>
	// Timer (debugging)
#include <CGAL/Timer.h>	
	// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include "C2t3/Complex_3_in_triangulation_3_file_writer.h"

// CImg include (only used here for easily dealing with input parameters)
#include "CImg.h"

// Project-specific includes
#include "Splat_3_Mesher/AABB_Tree_Splat_3_Mesher_Oracle.h"
#include "Splat_3_Surface_Mesher/AABB_Tree_Splat_3_Surface_Mesher_Oracle.h"
#include "C2t3/Complex_2_in_triangulation_3_post_processing.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT								FT ;
typedef K::Triangle_3						Triangle_3 ;
typedef CGAL::Splat_3<K>					Splat_3 ;
typedef K::Point_3							Point_3 ;
typedef K::Vector_3							Vector_3 ;

//typedef CGAL::Monge_via_jet_fitting< K, K >			Monge_via_jet_fitting ;
//typedef Monge_via_jet_fitting::Monge_form			Monge_form ;

// AABB Tree related definitions
typedef std::vector< Splat_3 >::iterator					SplatsIterator ;
typedef CGAL::AABB_Splat_3_primitive< K, SplatsIterator >		Primitive ;
typedef CGAL::AABB_Splat_3_traits< K, Primitive >				AABB_Splat_3_traits ;
typedef CGAL::AABB_tree< AABB_Splat_3_traits >			AABBTree ;
typedef AABBTree::Object_and_primitive_id				Object_and_primitive_id;

// Custom oracles
typedef CGAL::AABB_Tree_Splat_3_Mesher_Oracle< K, SplatsIterator > Mesher_Oracle ;
typedef CGAL::AABB_Tree_Splat_3_Surface_Mesher_Oracle< K, SplatsIterator > Surface_Mesher_Oracle ;

// Default triangulation for Mesher
typedef CGAL::Mesh_triangulation_3< Mesher_Oracle >::type	MTr ;
typedef CGAL::Mesh_complex_3_in_triangulation_3< MTr >		C3t3 ;
typedef CGAL::Mesh_criteria_3< MTr >						Mesh_criteria ;

// Default triangulation for Surface Mesher
typedef CGAL::Surface_mesh_default_triangulation_3			Tr ;
typedef CGAL::Complex_2_in_triangulation_3<Tr>				C2t3 ;



using namespace std ;

/* Main function */
int main ( int argc, char **argv) {

	cimg_usage("Error, wrong number of parameters!\nUsage: ./Splat_3_Mesher -ips <input_pointset_file_path> -od <output_off_dir_path> -ofn <(opt)output_off_file_name(will be automatically generated if not specifyed)> -ab <(opt)angular_bound_in_degrees> -rb <(opt)radius_bound> -db <(opt)distance_bound> -mf <(opt)ManifoldnessFlag> -it <(opt)intersectionType> -ds <(opt)distanceSigma> -ms <(opt)mogSigma> -nb <(opt)numBins>\n       Where: ManifoldnessFlag = 1 --> Manifold_tag\n              ManifoldnessFlag = 2 --> Manifold_with_boundary_tag\n              ManifoldnessFlag = 3 --> Non_manifold_tag\nNote: refer to the CGAL documentation, on the 3D Surface Mesh Generation section, for more information about the meaning of these parameters." ) ;

	// Parse/default input parameters
	const char*		inputSplatsFile			= cimg_option( "-i",		(char*)0,		"Input Splats file" ) ;
	const char*		outputFileDir			= cimg_option( "-od",		(char*)0,		"Output mesh directory" ) ;
	const char*		outputFileName			= cimg_option( "-ofn",		(char*)0,		"Output file name" ) ;
	const double	angularBound			= cimg_option( "-ab",		10.0,			"Angular Bound" ) ;
	const double	radiusBound				= cimg_option( "-rb",		0.1,			"Radius Bound" ) ;
	const double	distanceBound			= cimg_option( "-db",		0.1,			"Distance Bound" ) ;
	const int		manifoldFlag			= cimg_option( "-mf",		2,				"Manifoldness Flag" ) ;
	const double	ransacDistThres			= cimg_option( "-rdt",		0.1,			"Ransac distance threshold" ) ;
	const double	smallestRansacSegment	= cimg_option( "-srs",		distanceBound,	"Smallest segment to apply ransac (if segment.length is smaller, no ransac will be applied)" ) ;
	const bool		doOutputRedirection		= cimg_option( "-outTxt",	false,			"Output on-screen redirection to txt file." ) ;
	const double	distanceSigma			= cimg_option( "-ds",		0.25,			"Distance Sigma (only used when -it 1|2)" ) ;
	const int		mesherType				= cimg_option( "-mesher",	0,				"CGAL Mesher type: 0 = Surface Mesher / 1 = Mesher" ) ;
	const bool		doCleanResults			= cimg_option( "-clean",	true,			"Try some cleaning steps on the output surface" ) ;
	
	// Some input parameters check
	if ( manifoldFlag < 0 || manifoldFlag > 2 ) {
		cerr << "ManifoldnessFlag parameter must be 0, 1 or 2!" << endl ;
		return -1 ;
	}	
	if ( ransacDistThres < 0 ) {
		cerr << "Ransac distance threshold parameter must be > 0!" << endl ;
		return -1 ;
	}
	
	// Get the file type ( *.ops = oriented pointset / *.opss = oriented pointset with score )
	string inputFilePathStr( inputSplatsFile ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	string inputFileExtension = inputFilePathStr.substr( indexFileExtension+1, ( inputFilePathStr.size()-(indexFileExtension+1) ) ) ;
	//cout << inputFileExtension << endl ;
	if ( strcmp( inputFileExtension.c_str(), "splat" ) != 0 ) {
		cerr << "Unrecognized file type..." << endl ;
		return -1 ;
	}

	// Build the output file name
	std::string str ; 
	std::ostringstream oss( str ) ;
	std::ostringstream otxt ; 
	std::ostringstream orName ; // Debug only!
	std::ostringstream cleanName ; 
	if ( outputFileName == (char*)0 ) {
		// Get the input file name, to name the output one with the same name
		string inputFilePathStr( inputSplatsFile ) ; 
		size_t indexFileName = inputFilePathStr.find_last_of( "/\\" ) ;
		size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
		string inputFileNameStr = inputFilePathStr.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;
		
		// Build a name containing the parameters used
		oss << outputFileDir << "/" << inputFileNameStr 
									<< "__ab" << angularBound
									<< "_rb" << radiusBound
									<< "_db" << distanceBound
									<< "_mf" << manifoldFlag
									<< "_rdt" << ransacDistThres ;

		if ( doCleanResults ) {
			cleanName << oss.str() << "_clean.off" ;
		}

		if ( doOutputRedirection ) {
			otxt << oss.str() << ".txt" ;
		}

		oss << ".off" ;
		
	}
	else {
		if ( outputFileDir == (char*)0 ) {
			oss << outputFileName ;
		}
		else {
			oss << outputFileDir << "/" << outputFileName ;
		}

		string outputFileNamePath( outputFileName ) ; 
		size_t indexFileName = outputFileNamePath.find_last_of( "/\\" ) ;
		size_t indexFileExtension = outputFileNamePath.find_last_of( "." ) ;
		string outputFileNameStr = outputFileNamePath.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;

		cleanName << outputFileNamePath.substr( 0, indexFileName ) << "/" << outputFileNameStr << "_clean.off" ;

	}
	char* fullOutputFilePath ;
	fullOutputFilePath = new char [ oss.str().size()+1 ] ;
	strcpy ( fullOutputFilePath, oss.str().c_str() ) ;
	
	// Redirect the output if needed
	std::ofstream outTxt;
	if ( doOutputRedirection ) {
		cout << "Redirecting output to a txt file..." << endl ;
		outTxt.open( otxt.str().c_str() ) ;
		std::streambuf *coutbuf = std::cout.rdbuf() ; //save old buf
		std::cout.rdbuf( outTxt.rdbuf() ) ;
	}

	// Print parameters on screen
	cout << "Parameters:" << endl ;
	cout << "  Input file: " << inputSplatsFile << endl ;
	cout << "  Output file: " << fullOutputFilePath << endl ;
	cout << "  Angular Bound = " << angularBound << endl ;
	cout << "  Radius Bound = " << radiusBound << endl ;
	cout << "  Distance Bound = " << distanceBound << endl ;
	cout << "  Manifoldness Flag = " << manifoldFlag << endl ;
	cout << "  Ransac Distance Threshold = " << ransacDistThres << endl ;
	cout << "  Gaussian Sigma = " << distanceSigma << endl ;
		
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
	
	// AABBTree 
	cout << "Creating the AABB Tree..." ;
	CGAL::AABB_tree< AABB_Splat_3_traits >* tree = new CGAL::AABB_tree< AABB_Splat_3_traits >( splats.begin(), splats.end() ) ;
	timer.stop() ;
	cout << "done (" << timer.time() << " s)" << endl ;
	timer.reset() ;
	
	/* Meshing */
	if ( mesherType == 0 ) {
		/* Using Surface_Mesh_3 */
		cout << "Using Surface Mesher" << endl ;

		Tr tr;            // 3D-Delaunay triangulation
		C2t3 c2t3 (tr);   // Initialize 2D-complex in 3D-Delaunay triangulation

		// Add the initial points to the triangulation
		std::vector< Point_3 > initialPoints ;
		{	
			// Insert a set of 10 random points from the original polygon into the triangulation		
			typedef std::vector< Point_3 >::size_type size_type ;
			size_type nb_initial_points = 10 ;
			nb_initial_points = (std::min)( nb_initial_points, splats.size() ) ;
			for( size_type n = 0; n < nb_initial_points; n++ ) {				 
				const int pos = CGAL::default_random.get_int( 0, static_cast< int >( splats.size() ) ) ;
				//const int pos = CGAL::default_random.get_int( 0, 1000 ) ; // WARNING!!! Debug only!!!
				initialPoints.push_back( splats[ pos ].center() ) ;
			}
		}
		
		// defining meshing criteria
		CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angularBound,		// angular bound
															radiusBound,		// radius bound
															distanceBound ) ;	// distance bound
	
		// Construct the oracle
		Surface_Mesher_Oracle oracle( tree,
									  initialPoints,
									  ransacDistThres, 
									  smallestRansacSegment ) ;
	
		// Mesh the surface
		cout << "Meshing the Splats..." << flush ;		
		timer.start() ;
		switch ( manifoldFlag ) {
			case 0:
				CGAL::make_surface_mesh( c2t3, oracle, oracle, criteria, CGAL::Manifold_tag() ) ;
				break ;
			case 1:
				CGAL::make_surface_mesh( c2t3, oracle, oracle, criteria, CGAL::Manifold_with_boundary_tag() ) ;
				break ;
			case 2:
				CGAL::make_surface_mesh( c2t3, oracle, oracle, criteria, CGAL::Non_manifold_tag() ) ;
				break ;
			default: 
				// Is this possible?
				cerr << "ManifoldnessFlag parameter must be 1, 2 or 3!" << endl ;
				return -1 ;
				break ;
		}	
		timer.stop() ;
		cout << "done (" << timer.time() << " s)" << endl ;
		timer.reset() ;
	
		cout << "  Number of resulting vertices on surface " << tr.number_of_vertices() << endl ;
		cout << "  Number of resulting facets on surface " << c2t3.number_of_facets() << endl ;
	
		// Save the results
		cout << "Saving the result..." ;		
		timer.start() ;
		std::ofstream out( fullOutputFilePath ) ; 
		CGAL::output_surface_facets_to_off( out, c2t3 ) ; 
		out.close() ;

		timer.stop() ;
		cout << "done (" << timer.time() << " s)" << endl ;
		timer.reset() ;

		
		if ( manifoldFlag == 2 && doCleanResults ) {
			// Apply battery of cleaning steps (only if non-manifold configurations are found...)
			int numNonManifoldEdges = 0, numNonManifoldVertices = 0 ;
			CGAL::C2t3_PostProcessing::manifoldnessCheck( c2t3, numNonManifoldEdges, numNonManifoldVertices ) ;
			if ( numNonManifoldEdges > 0 || numNonManifoldVertices > 0 ) {
				cout << "Cleaning the result..." << flush ;		
				timer.start() ;
		
				int numSlivers = CGAL::C2t3_PostProcessing::removeSlivers( c2t3 ) ;
				int numEars = CGAL::C2t3_PostProcessing::removeEars( c2t3 ) ;
				int totalNMVert = 0 ;
				int numNMVert = 9999999 ;
				int nmvRounds = 0 ;
				while ( numNMVert > 0 ) {
					numNMVert = CGAL::C2t3_PostProcessing::removeNonManifoldVertices( c2t3 ) ;
					totalNMVert += numNMVert ;
					nmvRounds++ ;
				}	

				timer.stop() ;
				cout << "done (" << timer.time() << " s)" << endl ;
				timer.reset() ;

				cout << "  Removed facets from slivers configurations " << numSlivers << endl ;
				cout << "  Removed facets from ears configurations " << numEars << endl ;
				cout << "  Removed non-manifold vertices " << totalNMVert << "(" << nmvRounds << " iterations)" <<  endl ;
				cout << "  Number of facets after cleaning " << c2t3.number_of_facets() << endl ;
		
				// Manifoldness check
				CGAL::C2t3_PostProcessing::manifoldnessCheck( c2t3, numNonManifoldEdges, numNonManifoldVertices ) ;
				cout << "  Resulting non-Manifold edges " << numNonManifoldEdges << endl ;
				cout << "  Resulting non-Manifold vertices " << numNonManifoldVertices << endl ;

				cout << "Saving the cleaned result..." << flush ;		
				timer.start() ;
				std::ofstream outClean( cleanName.str().c_str() ) ; 
				CGAL::output_surface_facets_to_off( outClean, c2t3 ) ; 
				out.close() ;
			}
		}	

	}
	else {
		/* Mesh_3 */
		cout << "Using Volume Mesher" << endl ;

		MTr tr ;            // 3D-Delaunay triangulation
		// C3t3 c3t3( tr ) ;   // Initialize 2D-complex in 3D-Delaunay triangulation	
		
		// Add the initial points to the triangulation
		std::vector< Point_3 > initialPoints ;
		{	
			// Insert a set of 10 random points from the original polygon into the triangulation		
			typedef std::vector< Point_3 >::size_type size_type ;
			size_type nb_initial_points = 10 ;
			nb_initial_points = (std::min)( nb_initial_points, splats.size() ) ;
			for( size_type n = 0; n < nb_initial_points; n++ ) {				 
				const int pos = CGAL::default_random.get_int( 0, static_cast< int >( splats.size() ) ) ;
				//const int pos = CGAL::default_random.get_int( 0, 1000 ) ; // WARNING!!! Debug only!!!
				initialPoints.push_back( splats[ pos ].center() ) ;
			}
		}
		
		// Define the meshing criteria
		Mesh_criteria criteria( CGAL::parameters::facet_angle = angularBound, 
								CGAL::parameters::facet_size = radiusBound, 
								CGAL::parameters::facet_distance = distanceBound, CGAL::parameters::edge_size = 0.025 ) ;
								// CGAL::parameters::cell_radius_edge_ratio = 3.0, 
								// CGAL::parameters::cell_size = 3.0 ) ;

		// Construct the oracle
		cout << "Creating the oracle..." << flush ;		
		Mesher_Oracle oracle( tree, 	
							  initialPoints,
							  distanceSigma, 
							  ransacDistThres, 
							  smallestRansacSegment ) ; 
		cout << "done" << endl ;

		// Mesh generation
		cout << "Meshing the Splats representation..." << flush ;		
		timer.start() ;
		C3t3 c3t3 = CGAL::make_mesh_3< C3t3 >( oracle, criteria, CGAL::parameters::no_perturb(), CGAL::parameters::no_exude() ) ;
		timer.stop() ;
		cout << "done (" << timer.time() << " s)" << endl ;
		timer.reset() ;

		// Save the results
		cout << "Saving results..." ;		
		timer.start() ;
		std::ofstream out( fullOutputFilePath ) ; 	
		CGAL::output_subdomain_surface_facets_to_off( out, c3t3, 1, 1 ) ; 
		out.close() ;

		timer.stop() ;
		cout << "done (" << timer.time() << " s)" << endl ;
		timer.reset() ;
	
		if ( doOutputRedirection ) {
			outTxt.close() ;
		}
	}

	return 0 ;
}
