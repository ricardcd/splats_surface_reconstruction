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
// Scale estimation
#include "RobustStatistics/ScaleEstimation.h"
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

// Default triangulation for Mesher
typedef CGAL::Surface_mesh_default_triangulation_3	Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr>		C2t3;

// Implicit Surface type
typedef OrientedSplatsDistanceFunction<K>									OrientedSplatsDistance ;
typedef CGAL::Implicit_surface_3< K, 
								   OrientedSplatsDistance >			Surface_3 ;

// Other definitions
using namespace std ;

int removeHallucinatedTriangles2( C2t3 &mesh, const OrientedSplatsDistance &func, double t ) ; // Defined below

/* Main function */
int main ( int argc, char **argv) {

	// Parse/default input parameters
	const char*		inputSplatsFile			= cimg_option( "-i",		(char*)0,		"Input oriented splats file" ) ;	
	const char*		outputFileDir			= cimg_option( "-od",		(char*)0,		"Output mesh directory" ) ;
	const char*		outputFileName			= cimg_option( "-ofn",		(char*)0,		"Output file name" ) ;
	const double	angularBound			= cimg_option( "-ab",		10.0,			"Angular bound" ) ;
	const double	radiusBound				= cimg_option( "-rb",		0.1,			"Radius bound" ) ;
	const double	distanceBound			= cimg_option( "-db",		0.1,			"Distance bound" ) ;
	const int		manifoldFlag			= cimg_option( "-mf",		2,				"Manifoldness Flag" ) ;	
	const int		k						= cimg_option( "-k",		5,				"Number of closest splats to take into account. If set to < 0, then all the splats inside the radial neighborhood defined by gh will be taken into account." ) ;
	const double	gaussianRadiusFactor	= cimg_option( "-gh",		0.2,			"Blending Gaussian H factor (w.r.t. the BSR)" ) ;
	const bool		doOutputRedirection		= cimg_option( "-outTxt",	false,			"Output on-screen redirection to txt file" ) ;		
	const bool		doCleanResults			= cimg_option( "-clean",	true,			"Try some cleaning steps on the output surface" ) ;	
	const double	epsilon					= cimg_option( "-ep",		0.001,			"Error epsilon" ) ;
	
	// Some input parameters check
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
									<< "_k" << k 
									<< "_grf" << gaussianRadiusFactor ;
									
		if ( doOutputRedirection ) {
			otxt << oss.str() << ".txt" ;
		}

		osClean << oss.str() << "_clean.off" ;
		osClean2 << oss.str() << "_ch.off" ;

		oss << ".off" ;
	}
	else {
		oss << outputFileDir << "/" << outputFileName ;		
		osClean << oss.str() << "_clean.off" ;
		osClean2 << oss.str() << "_ch.off" ;
		otxt << oss.str() << ".txt" ;

		oss << ".off" ;
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
	cout << "  Input file: " << inputSplatsFile ;
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


	if ( doOutputRedirection ) {
		outTxt.close() ;
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
