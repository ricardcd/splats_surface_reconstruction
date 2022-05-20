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
#include "CGAL/Monge_via_jet_fitting.h"
#include "SplatsIsosurfaceMesher/SplatsDistanceFunctionGCut.h"
	// General
#include <CGAL/Random.h>
	// Timer (debugging)
#include <CGAL/Timer.h>	
// CImg include (only used here for easily dealing with input parameters)
// #include "CImg.h"
// Boost
#include "boost/program_options.hpp"
// Noise scale estimation
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
typedef SplatsDistanceFunctionGCut<K>									SplatsDistance ;
typedef CGAL::Implicit_surface_3< K, 
								   SplatsDistance >			Surface_3 ;

// Other definitions
using namespace std ;
namespace po = boost::program_options ;



/* Post-processing (declaration) */
// int removeForBoundaries( C2t3 &mesh, const SplatsDistance &func, const FT &uncertainBand ) ;
// int removeHallucinatedTriangles( C2t3 &mesh, const SplatsDistance &func, const FT &radiusPart ) ;
int removeHallucinatedTriangles2( C2t3 &mesh, const SplatsDistance &func, double t ) ;

/* Main function */
int main ( int argc, char **argv) {

	// Parse/default input parameters
	std::string inputSplatsFile, outputFileDir, outputFileName, outputRunTimesFile ;
	double angularBound, radiusBound, distanceBound, gaussianRadiusFactor, functionErrorBound, radiusEdgeBound, ransacThres, confidenceThres, uncertainBand, smoothPrior, areaConstant, stWeight, dataWeight, filterRad, bandRadius, epsilon, tScale ;
    int manifoldFlag, k, numRays, bandKnn, odepth ;
	bool doCleanResults, boundedSurface, signPointsOutsideBand ;

	po::options_description options("Splats distance mesher (GraphCuts)");
	options.add_options()
		( "help", "Produce help message" )
        ("inFile,i", po::value<std::string>(&inputSplatsFile), "Input point set file path.")
        ("outFile,o", po::value<std::string>(&outputFileName), "Output mesh file path (OFF format).")
        ("outDir", po::value<std::string>(&outputFileDir), "Output directory. If outFile is not specified, a file with a default name containig a list of the parameters used will be written in this directory.")
		( "gaussianRadiusFactor", po::value<double>(&gaussianRadiusFactor)->default_value(0.2), "Blending Gaussian H factor. If 0 < bandRadius < 1, it is a multiplicative factor of the Bounding Sphere Radius. Otherwise, if bandRadius > 1, it is a multiplicative radius to the average spacing between points" )
        ( "bandRadius", po::value<double>(&bandRadius)->default_value(0.4), "Size of the band where the distance function will be computed.. If 0 < bandRadius < 1, it is a multiplicative factor of the Bounding Sphere Radius. Otherwise, if bandRadius > 1, it is a multiplicative radius to the average spacing between points" )
        ( "bandKnn", po::value<int>(&bandKnn)->default_value(-1), "Fix the K-NN used in the distance function computation for the secondary band (-1 = not used)" )
        ( "functionErrorBound", po::value<double>(&functionErrorBound)->default_value(1e-3), "Maximum allowed error bound between the real implicit function and its linear approximation in the triangulation" )
		( "radiusEdgeBound", po::value<double>(&radiusEdgeBound)->default_value(2.5), "Radius-edge bound for the implicit function" )
		( "fixKnn", po::value<int>(&k)->default_value(-1), "Number of closest splats to take into account. If set to < 0, then all the splats inside the radial neighborhood defined by gh will be taken into account." )
		( "numRays", po::value<int>(&numRays)->default_value(50), "Number of Random rays to try during in/out check" )
		( "ransacThres", po::value<double>(&ransacThres)->default_value(0.1), "RANSAC-intersection threshold" )
		( "confidenceThres", po::value<double>(&confidenceThres)->default_value(0.7), "Confidence threshold used in initial in/out labelling" )
		( "uncertainBand", po::value<double>(&uncertainBand)->default_value(0.01), "Uncertain band (points under this value in the unsigned function will not be taken into account" )
		( "smoothPrior", po::value<double>(&smoothPrior)->default_value(4.0), "Smooth prior on the unsigned function" )
		( "areaConstant", po::value<double>(&areaConstant)->default_value(1e-5), "Area constant on the unsgned function" )
		( "stWeight", po::value<double>(&stWeight)->default_value(1e5), "Weight to apply to S-T links (multiplies a confidence value ranging from 0 to 1)" )
		( "dataWeight", po::value<double>(&dataWeight)->default_value(1.0), "Weight to apply to Data links" )
		( "boundedSurface", po::value<bool>(&boundedSurface)->default_value(false), "Toggles the bounded surfaces trick" )
		( "outputRunTimesFile", po::value<std::string>(&outputRunTimesFile)->default_value(""), "Output run times file name (if empty, no log will be made)" )
		( "filterRad", po::value<double>(&filterRad)->default_value(0.2), "Hallucinated triangles filter radius" )
		( "epsilon", po::value<double>(&epsilon)->default_value(0.001), "Epsilon value used in distance comparisons" )
		( "signAllPointsOutsideBand", po::value<bool>(&signPointsOutsideBand)->default_value(false), "Toggles the signing confidence on ALL points outside the band (more costly!)" )
		( "hallucinatedTrisFactor", po::value<double>(&tScale)->default_value(2.5), "Multiplicative factor of the scale to eliminate hallucinated triangles" )
        ( "octreeDepth", po::value<int>(&odepth)->default_value(-1), "Octree depth used to filter the input point set before implicit distance creation. If smaller than one, no octree will be used (default = -1)" )
        ( "ab", po::value<double>(&angularBound)->default_value(0.0), "[Surface Mesher] Angular bound" )
        ( "rb", po::value<double>(&radiusBound)->default_value(0.1), "[Surface Mesher] Radius bound" )
        ( "db", po::value<double>(&distanceBound)->default_value(0.1), "[Surface Mesher] Distance bound"  )
        ( "mf", po::value<int>(&manifoldFlag)->default_value(2), "[Surface Mesher] Manifoldness Flag" )
        ( "clean", po::value<bool>(&doCleanResults)->default_value(true), "[Surface Mesher] Try some cleaning steps on the output surface" )
	;
		
	// Read parameters from command line
	po::variables_map vm ;
	po::store( po::parse_command_line(argc, argv, options ), vm ) ;
	// Read the configuration file (if exists)
	ifstream ifs( "params.cfg" ) ;
	po::store( po::parse_config_file( ifs, options ), vm ) ;

	po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << options << "\n";
        return 1;
    }
	
	// --- Debug (Start) ---
	//std::cout << "inputSplatsFile = " << inputSplatsFile << std::endl ;
	//std::cout << "outputFileName = " << outputFileName << std::endl ;
	//std::cout << "outputFileDir = " << outputFileDir << std::endl ;
	//std::cout << "angularBound = " << angularBound << std::endl ;
	//std::cout << "radiusBound = " << radiusBound << std::endl ;
	//std::cout << "distanceBound = " << distanceBound << std::endl ;
	//std::cout << "manifoldFlag = " << manifoldFlag << std::endl ;
	//std::cout << "doOutputRedirection = " << doOutputRedirection << std::endl ;
	//std::cout << "doCleanResults = " << doCleanResults << std::endl ;
	//std::cout << "gaussianRadiusFactor = " << gaussianRadiusFactor << std::endl ;
	//std::cout << "areaConstant = " << areaConstant << std::endl ;
	//std::cout << "functionErrorBound = " << functionErrorBound << std::endl ;
	//std::cout << "radiusEdgeBound = " << radiusEdgeBound << std::endl ;
	// --- Debug (Start) ---

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
									<< "_rq" << numRays 
									<< "_rt" << ransacThres 
									<< "_ct" << confidenceThres
									<< "_ub" << uncertainBand
									<< "_sp" << smoothPrior
									<< "_ac" << areaConstant 
									<< "_dw" << dataWeight
									<< "_stw" << stWeight
									<< "_feb" << functionErrorBound
									<< "_reb" << radiusEdgeBound ;

		osClean << oss.str() << "_c.off" ;
		osClean2 << oss.str() << "_ch.off" ;

		oss << ".off" ;
		
	}
	else {
		oss << outputFileName ;
		osClean << oss.str() << "_c.off" ;
		osClean2 << oss.str() << "_ch.off" ;
	}
	char* fullOutputFilePath ;
	fullOutputFilePath = new char [ oss.str().size()+1 ] ;
	strcpy ( fullOutputFilePath, oss.str().c_str() ) ;
	
	// Print parameters on screen
	cout << "- Parameters:" << endl ;
	cout << "    Input file: " << inputSplatsFile << endl ;
	cout << "    Output file: " << fullOutputFilePath << endl ;
	cout << "    Angular Bound = " << angularBound << endl ;
	cout << "    Radius Bound = " << radiusBound << endl ;
	cout << "    Distance Bound = " << distanceBound << endl ;
	cout << "    Manifoldness Flag = " << manifoldFlag << endl ;
	cout << "    Blending K-NN = " << k << endl ;
	cout << "    Blending Gaussian H (w.r.t BSR) = " << gaussianRadiusFactor << endl ;
	cout << "    Ray queries = " << numRays << endl ;
	cout << "    Ransac threshold = " << ransacThres << endl ;
	cout << "    Confidence threshold = " << confidenceThres << endl ;
	cout << "    Uncertain band = " << uncertainBand << endl ;
	cout << "    Smooth prior = " << smoothPrior << endl ;
	cout << "    Area constant = " << areaConstant << endl ;
	cout << "    S-T weight = " << stWeight << endl ;
	cout << "    Data weight = " << dataWeight << endl ;
	cout << "    Function error bound = " << functionErrorBound << endl ;
	cout << "    Radius-Edge bound = " << radiusEdgeBound << endl ;
    if ( odepth > 0 ) {
        cout << "    Using octree of depth = " << odepth << endl ;
    }
	if ( boundedSurface ) 
		cout << "    Using the open surface trick" << endl ;
	
	// Initialize the timing log, if needed
	bool logTimings = !outputRunTimesFile.empty() ;
	std::ofstream otime ;
	std::ofstream oGlobalTime ;
	CGAL::Timer timer ; // Reusable timer initialization
	CGAL::Timer globalTimer ;
	double t = 0 ; // We store timings in this variable...
	if ( logTimings )
		otime.open( outputRunTimesFile.c_str(), std::ios_base::out ) ;

	globalTimer.start() ;		
	
	// Loading features
	
	// Read the input pointset file	
	cout << "- Reading splats from file..." << flush ;
	timer.start() ;
	// Opening file
	ifstream file( inputSplatsFile.c_str(), ifstream::in ) ;
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
	t = timer.time() ;
	cout << "done (" << t << " s), readed " << splats.size() << " splats." << endl ;
	if ( logTimings ) otime << "Reading splats = " << t << "s" << endl ;
	
	// Define the surface
	cout << "- Computing the signed distance function:" << endl ;
	timer.reset() ;
	SplatsDistance distanceFun( splats, 
								k, 
								gaussianRadiusFactor, 
								ransacThres, 
								numRays, 
								confidenceThres,
								uncertainBand,
								smoothPrior,
								areaConstant, 
								stWeight, 
								dataWeight, 
								functionErrorBound,
								radiusEdgeBound, 
								boundedSurface,
								bandRadius,
								bandKnn,
								epsilon,
                                signPointsOutsideBand,
                                odepth ) ;

	t = timer.time() ;
	if ( logTimings ) {
		otime << "Creating SDF = " << t << "s" << endl ;
		otime << "    - Creating distance function = " << distanceFun.getTimeImplicitFunctionCreation() << "s" << endl ;
		otime << "    - Signing distance function = " << distanceFun.getTimeSigningFunction() << "s" << endl ;
	}

	if ( !distanceFun.correct() ) {
		// Some error during implicit function signing...
		std::cout << "[ERROR!] Cannot mesh the implicit function.\n[ERROR!] ABORTING! " << std::endl ;
		return -1 ;
	}

	// Creating the surface
	Surface_3 surface( distanceFun,										// Pointer to function
					   distanceFun.enlarged_bounding_sphere(1.5) ) ;	// Bounding sphere, where the implicit function is computed

	// Mesh generation
	cout << "- Meshing the Splats representation..." << flush ;
	
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angularBound,		// Angular bound
														radiusBound,		// Radius bound
														distanceBound );	// Distance bound

	timer.reset() ;
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
			cerr << "[ERROR!] ManifoldnessFlag parameter must be 1, 2 or 3!" << endl ;
			return -1 ;
			break ;
	}
	t = timer.time() ;
	cout << "done (" << t << " s)" << endl ;
	if ( logTimings ) otime << "Surface Extraction (SM) = " << t << "s" << endl ;	

	// Save the results
	cout << "- Saving results..." ;		
	timer.reset() ;
	std::ofstream out( fullOutputFilePath ) ;
	CGAL::output_surface_facets_to_off( out, c2t3 ) ;
	out.close() ;

	t = timer.time() ;
	cout << "done (" << t << " s)" << endl ;
	if ( logTimings ) otime << "Saving results = " << t << "s" << endl ;	
	
	// Perform Post processing (if required)
	if ( doCleanResults ) {

		cout << "- Cleaning hallucinated triangles..." << flush ;
		timer.reset() ;
		// int numHalluTri = removeForBoundaries( c2t3, distanceFun, uncertainBand ) ;
		// int numHalluTri = removeHallucinatedTriangles( c2t3, distanceFun, filterRad ) ;
		int numHalluTri = removeHallucinatedTriangles2( c2t3, distanceFun, tScale ) ;
		t = timer.time() ;
		cout << "done (" << t << " s), removed " << numHalluTri << " triangles." << endl ;
		if ( logTimings ) otime << "Removing hallucinated = " << t << "s" << endl ;	

		cout << "- Saving clean mesh..." << flush ;
		timer.reset() ;
		std::ofstream outClean( osClean2.str().c_str() ) ;
		CGAL::output_surface_facets_to_off( outClean, c2t3 ) ;
		outClean.close() ;
		t = timer.time() ;
		cout << "done (" << t << " s)" << endl ;
		if ( logTimings ) otime << "Saving cleaned hallucinated = " << t << "s" << endl ;	

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

			t = timer.time() ;
			cout << "done (" << t << " s)" << endl ;
			if ( logTimings ) otime << "Removing non-manifold = " << t << "s" << endl ;

			cout << "- Saving clean mesh..." << flush ;
			timer.reset() ;
			std::ofstream outClean( osClean.str().c_str() ) ; 
			CGAL::output_surface_facets_to_off( outClean, c2t3 ) ; 
			outClean.close() ;
			t = timer.time() ;
			cout << "done (" << t << " s)" << endl ;
			if ( logTimings ) otime << "Saving clean non-manifold = " << t << "s" << endl ;
		}
		else {
			std::cout << "- Resulting surface is manifold, no cleaning required." << std::endl ;
		}
	}

	if ( logTimings ) {
		if ( logTimings ) otime << "GLOBAL Time = " << globalTimer.time() << "s" << endl ;
		otime.close() ;
	}

	return 0 ;
}



/* Post-processing (implementation) */
//int removeForBoundaries( C2t3 &mesh, const SplatsDistance &func, const FT &uncertainBand ) {
//
//	typedef C2t3::Vertex_iterator							Vertex_iterator ;
//	typedef Tr::Facet										Facet ;
//	typedef Tr::Triangulation_data_structure				Tds ;
//	typedef C2t3::Facet_iterator							Facet_iterator ;	
//
//	// Run over all facets of the triangulation
//	Tr tr = mesh.triangulation() ;
//	Tds tds = tr.tds() ;
//
//	Facet_iterator fit ;
//	int numRemoved = 0 ;
//	for ( fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit ) {
//		if ( mesh.face_status( *fit ) != C2t3::NOT_IN_COMPLEX ) {
//
//			Triangle_3 t = tr.triangle( *fit ) ;
//
//			/*std::cout << "func( t.vertex(0) ) = " << func.distanceToClosestSplatCenter( t.vertex(0) ) << std::endl ;
//			std::cout << "func( t.vertex(1) ) = " << func.distanceToClosestSplatCenter( t.vertex(1) ) << std::endl ;
//			std::cout << "func( t.vertex(2) ) = " << func.distanceToClosestSplatCenter( t.vertex(2) ) << std::endl ;*/
//
//			if ( func.distanceToClosestSplatCenter( t.vertex(0) ) > ( uncertainBand*uncertainBand ) ||
//				 func.distanceToClosestSplatCenter( t.vertex(1) ) > ( uncertainBand*uncertainBand ) ||
//				 func.distanceToClosestSplatCenter( t.vertex(2) ) > ( uncertainBand*uncertainBand ) )
//			{
//				mesh.remove_from_complex( *fit ) ;
//				numRemoved++ ;
//			}
//
//		}
//	}
//
//	return numRemoved ;
//}

/* Post-processing (implementation) */
//int removeHallucinatedTriangles( C2t3 &mesh, const SplatsDistance &func, const FT &radiusPart ) {
//
//	FT rad = radiusPart ;
//	//if ( rad > 1 || rad < 0 ) {
//	//	std::cout << "[WARNING] Hallucinated triangles radius must be between 0 and 1! Defaulting to 0.5..." << std::endl ;
//	//	rad = 0.5 ;
//	//}
//
//	typedef C2t3::Vertex_iterator							Vertex_iterator ;
//	typedef Tr::Facet										Facet ;
//	typedef Tr::Triangulation_data_structure				Tds ;
//	typedef C2t3::Facet_iterator							Facet_iterator ;	
//
//	// Run over all facets of the triangulation
//	Tr tr = mesh.triangulation() ;
//	Tds tds = tr.tds() ;
//
//	Facet_iterator fit ;
//	int numRemoved = 0 ;
//	for ( fit = mesh.facets_begin(); fit != mesh.facets_end(); ++fit ) {
//		if ( mesh.face_status( *fit ) != C2t3::NOT_IN_COMPLEX ) {
//
//			Triangle_3 t = tr.triangle( *fit ) ;
//
//			// std::cout << t << std::endl ;
//						
//			if ( !func.allPointsOnRestrictedBand( t, rad ) )
//			{
//				mesh.remove_from_complex( *fit ) ;
//				numRemoved++ ;
//			}
//		}
//	}
//
//	return numRemoved ;
//}

int removeHallucinatedTriangles2( C2t3 &mesh, const SplatsDistance &func, double t ) {

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
