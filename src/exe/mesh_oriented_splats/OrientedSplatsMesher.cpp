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
	// Mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h> 
#include "C2t3/Complex_2_in_triangulation_3_post_processing.h"
	// AABB_tree
#include "Splat_3/Splat_3.h"
#include "Splat_3/SplatsIO.h"
#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
#include "SplatsIsosurfaceMesher/OrientedSplatsDistanceFunction.h"
	// General
#include <CGAL/Random.h>
	// Timer (debugging)
#include <CGAL/Timer.h>
    // Boost
#include <boost/program_options.hpp>
// Scale estimation
#include "RobustStatistics/ScaleEstimation.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT								FT ;
typedef K::Triangle_3						Triangle_3 ;
typedef CGAL::Splat_3<K>					Splat_3 ;
typedef K::Point_3							Point_3 ;
typedef K::Vector_3							Vector_3 ;

typedef CGAL::Monge_via_jet_fitting< K, K >			Monge_via_jet_fitting ;
typedef Monge_via_jet_fitting::Monge_form			Monge_form ;

// Default triangulation for Mesher
typedef CGAL::Surface_mesh_default_triangulation_3	Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr>		C2t3;

// Implicit Surface type
typedef OrientedSplatsDistanceFunction<K>									OrientedSplatsDistance ;
typedef CGAL::Implicit_surface_3< K, 
								   OrientedSplatsDistance >			Surface_3 ;

// Other definitions
using namespace std ;
namespace po = boost::program_options;


int removeHallucinatedTriangles2( C2t3 &mesh, const OrientedSplatsDistance &func, double t ) ; // Defined below

/* Main function */
int main ( int argc, char **argv) {

	// Parse/default input parameters
	std::string inputSplatsFile, outputFileName, outputFileDir;
	double angularBound, radiusBound, distanceBound, gaussianRadiusFactor, epsilon;
	int manifoldFlag, k;
	bool doCleanResults = false;

    po::options_description options("Meshes a set of splats, assuming that the input is properly oriented. It computes a signed distance function (SDF) on an adaptive data structure, and then meshes it.");
    options.add_options()
            ("help,h", "Produce help message")
            ("inFile,i", po::value<std::string>(&inputSplatsFile), "Input point set file path.")
            ("outFile,o", po::value<std::string>(&outputFileName), "Output mesh file path (OFF format).")
            ("outDir", po::value<std::string>(&outputFileDir), "Output directory. If outFile is not specified, a file with a default name containig a list of the parameters used will be written in this directory.")
            ("ab", po::value<double>(&angularBound)->default_value(10), "Surface Mesher's Angular bound")
            ("rb", po::value<double>(&radiusBound)->default_value(0.1), "Surface Mesher's Radius bound")
            ("db", po::value<double>(&distanceBound)->default_value(0.1), "Surface Mesher's Distance bound")
            ("mf", po::value<int>(&manifoldFlag)->default_value(2), "Surface Mesher's manifoldness flag (0 == Manifold, 1 == Manifold with boundary, 2 = Non-manifold)")
            ("knn,k", po::value<int>(&k)->default_value(5), "Number of nearby splats to take into account when computing the signed distance function.")
            ("gh", po::value<double>(&gaussianRadiusFactor)->default_value(0.2), "Blending Gaussian H factor (w.r.t. the Bonding Sphere Radius)")
            ("clean", po::bool_switch(&doCleanResults), "If set, will try some cleaning steps on the output surface")
            ("epsilon", po::value<double>(&epsilon)->default_value(0.001), "Difference tolerance when comparing values in the SDF." )
        ;

    // Read parameters from command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << options << "\n";
        return 1;
    }

	// Some input parameters check
    if (outputFileName.empty() && outputFileDir.empty()) {
        std::cout << "[ERROR] Please set either --outFile or --outDir parameters (or both)." << std::endl;
        return -1;
    }
	if ( manifoldFlag < 0 || manifoldFlag > 2 ) {
		cerr << "ManifoldnessFlag parameter must be 0, 1 or 2!" << endl ;
		return -1 ;
	}	
	
	// Get the file type (*.splat)
	string inputFilePathStr( inputSplatsFile ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	string inputFileExtension = inputFilePathStr.substr( indexFileExtension+1, ( inputFilePathStr.size()-(indexFileExtension+1) ) ) ;	
	if ( strcmp( inputFileExtension.c_str(), "splat" ) != 0 ) {
		cerr << "Unrecognized file type..." << endl ;
		return -1 ;
	}

	// Build the output file name
	std::string str ; 
	std::ostringstream oss( str ) ;
	std::ostringstream otxt ; 
	std::ostringstream orName ; // Debug only!
	std::ostringstream osClean ;
	std::ostringstream osClean2 ;
	if ( outputFileName.empty() ) {
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
									<< "_k" << k 
									<< "_grf" << gaussianRadiusFactor ;
									
		osClean << oss.str() << "_clean.off" ;
		osClean2 << oss.str() << "_ch.off" ;

		oss << ".off" ;
	}
	else {
        if ( outputFileDir.empty() ) {
            oss << outputFileName ;
        }
        else {
            oss << outputFileDir << "/" << outputFileName;
            osClean << oss.str() << "_clean.off";
            osClean2 << oss.str() << "_ch.off";
            otxt << oss.str() << ".txt";
        }
	}
	std::string fullOutputFilePath =oss.str();
	
	// Print parameters on screen
	cout << "Parameters:" << endl ;
	cout << "  Input file: " << inputSplatsFile << endl ;
	cout << "  Output file: " << fullOutputFilePath << endl ;
	cout << "  Angular Bound = " << angularBound << endl ;
	cout << "  Radius Bound = " << radiusBound << endl ;
	cout << "  Distance Bound = " << distanceBound << endl ;
	cout << "  Manifoldness Flag = " << manifoldFlag << endl ;
	cout << "  Blending K-NN = " << k << endl ;
	cout << "  Blending Gaussian H (w.r.t BSR) = " << gaussianRadiusFactor << endl ;

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
	
	// Define the surface
	OrientedSplatsDistance distanceFun( splats, 
										k, 
										gaussianRadiusFactor,
										epsilon ) ;

	Surface_3 surface( distanceFun,										// Pointer to function
                       distanceFun.enlarged_bounding_sphere(1.5) ) ;	// Bounding sphere, where the implicit function is computed

	// Mesh generation
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angularBound,		// Angular bound
														radiusBound,		// Radius bound
														distanceBound );	// Distance bound
	
	cout << "Meshing the Splats representation..." << flush ;
	timer.start() ;
	switch ( manifoldFlag ) {
		case 0:
			CGAL::make_surface_mesh( c2t3, 
										surface, 
										criteria, 
										CGAL::Manifold_tag() ) ;
			break ;
		case 1:
			CGAL::make_surface_mesh( c2t3, 
										surface, 
										criteria, 
										CGAL::Manifold_with_boundary_tag() ) ;
			break ;
		case 2:
			CGAL::make_surface_mesh( c2t3, 
										surface, 
										criteria, 
										CGAL::Non_manifold_tag() ) ;
			break ;
		default:
			cerr << "ManifoldnessFlag parameter must be 1, 2 or 3!" << endl ;
			return -1 ;
			break ;
	}
	timer.stop() ;
	cout << "done (" << timer.time() << " s)" << endl ;
	timer.reset() ;

	// Save the results
	cout << "Saving results..." ;		
	timer.start() ;
	std::ofstream out( fullOutputFilePath ) ;
	CGAL::output_surface_facets_to_off( out, c2t3 ) ;
	out.close() ;

	timer.stop() ;
	cout << "done (" << timer.time() << " s)" << endl ;
	timer.reset() ;
	
	double t = 0.0 ;

	// Perform Post processing (if required)
	if ( doCleanResults ) {
		
		cout << "- Cleaning hallucinated triangles..." << flush ;
		timer.reset() ;
		// int numHalluTri = removeHallucinatedTriangles( c2t3, *distanceFun, filterRad ) ;
		int numHalluTri = removeHallucinatedTriangles2( c2t3, distanceFun, 2.5 ) ;
		t = timer.time() ;
		cout << "done (" << t << " s), removed " << numHalluTri << " triangles." << endl ;
				
		cout << "- Saving clean mesh..." << flush ;
		timer.reset() ;
		std::ofstream outClean( osClean2.str().c_str() ) ;
		CGAL::output_surface_facets_to_off( outClean, c2t3 ) ;
		outClean.close() ;
		t = timer.time() ;
		cout << "done (" << t << " s)" << endl ;
		
		int numNonManifoldEdges = 0, numNonManifoldVertices = 0 ;
		if ( !CGAL::C2t3_PostProcessing::manifoldnessCheck< C2t3 >( c2t3, numNonManifoldEdges, numNonManifoldVertices ) ) {

			cout << "- Cleaning resulting mesh..." << flush ;		
			timer.reset() ;
			// Remove sliver triangles
			int numSlivers = CGAL::C2t3_PostProcessing::removeSlivers< C2t3 >( c2t3 ) ;
			// Remove ears
			int numEars = CGAL::C2t3_PostProcessing::removeEars< C2t3 >( c2t3 ) ;
			// Remove non-manifold vertices
			int totalNMVert = 0 ;
			int numNMVert = 9999999 ;
			int nmvRounds = 0 ;
			while ( numNMVert > 0 ) {
				numNMVert = CGAL::C2t3_PostProcessing::removeNonManifoldVertices< C2t3 >( c2t3 ) ;
				totalNMVert += numNMVert ;
				nmvRounds++ ;
			}	

			cout << "done (" << timer.time() << " s)" << endl ;

			cout << "- Saving clean mesh..." << flush ;
			timer.reset() ;
			std::ofstream outClean( osClean.str().c_str() ) ; 
			CGAL::output_surface_facets_to_off( outClean, c2t3 ) ; 
			outClean.close() ;
			cout << "done (" << timer.time() << " s)" << endl ;
		}
		else {
			std::cout << "- Resulting surface is manifold, no cleaning required." << std::endl ;
		}
	}

	return 0 ;
}





int removeHallucinatedTriangles2( C2t3 &mesh, const OrientedSplatsDistance &func, double t ) {

	typedef C2t3::Vertex_iterator							Vertex_iterator ;
	typedef C2t3::Facet_iterator							Facet_iterator ;	

	// Run over all vertices of the 2D complex
	std::vector< FT > values ;
	Vertex_iterator vit ;
	for ( vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) {
		values.push_back( func.unsigned_func_for_scale(vit->point()) ) ;
	}
	
	// Compute scale
	FT scale = ScaleEstimation::msSE( values, 1, 0.1, t ) ;
	FT thres = t*scale ;
	std::cout << "Scale = " << scale << " / Thres = " << thres << " / max = " << *( std::max_element( values.begin(), values.end()) ) << std::endl ;

	// Filter triangles based on the above computed measure
	Tr tr = mesh.triangulation() ;

	Facet_iterator fit ;
	int numRemoved = 0 ;
	for ( fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit ) {
		if ( mesh.face_status( *fit ) != C2t3::NOT_IN_COMPLEX ) {

			Triangle_3 t = tr.triangle( *fit ) ;

			// std::cout << t << std::endl ;
						
			if ( func.unsigned_func_for_scale( t.vertex( 0 ) ) > thres ) {
				mesh.remove_from_complex( *fit ) ;
				numRemoved++ ;
			}
			else if ( func.unsigned_func_for_scale( t.vertex( 1 ) ) > thres ) {
				mesh.remove_from_complex( *fit ) ;
				numRemoved++ ;
			}
			else if ( func.unsigned_func_for_scale( t.vertex( 2 ) ) > thres ) {
				mesh.remove_from_complex( *fit ) ;
				numRemoved++ ;
			}

		}
	}

	return numRemoved ;
}
