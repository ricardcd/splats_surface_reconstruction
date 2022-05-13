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
	// Polyhedron
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
	// AABB_tree
#include "Splat_3/Splat_3.h"
#include "Splat_3/SplatsIO.h"
#include "CGAL/Monge_via_jet_fitting.h"
#include "SplatsIsosurfaceMesher/SplatsDistanceFunctionNCut.h"
	// General
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>
	// Timer (debugging)
#include <CGAL/Timer.h>	
// Boost
#include <boost/program_options.hpp>
// Marching cubes
#include "SplatsIsosurfaceMesher/marching_cubes.h"
// Noise scale estimation
#include "RobustStatistics/ScaleEstimation.h"
// #include "SplatsIsosurfaceMesher/GenericFunction.h"


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
typedef SplatsDistanceFunctionNCut<K>									SplatsDistance ;
// typedef GenericFunction<K>									SplatsDistance ;
typedef CGAL::Implicit_surface_3< K, 
								   SplatsDistance >			Surface_3 ;

// Polyhedrons
typedef CGAL::Polyhedron_3<K>							Polyhedron ;

// Other definitions
using namespace std ;
namespace po = boost::program_options;



/* Post-processing (declaration) */
// int removeHallucinatedTriangles( C2t3 &mesh, const SplatsDistance &func, const FT &radiusPart ) ;
int removeHallucinatedTriangles2( C2t3 &mesh, const SplatsDistance &func, double t ) ;
int removeHallucinatedTrianglesMC( Polyhedron& poly, const SplatsDistance &func, const FT &radiusPart ) ;
// bool readDiscretizedDistanceFunction( std::istream& stream, std::vector< Point_3 >& vertices, std::vector< double > funcAtVert ) ;

/* Main function */
int main ( int argc, char **argv) {

	// Parse/default input parameters
	std::string inputFile, outputFileDir, outputFileName, outputRunTimesFile ;
	double angularBound, radiusBound, distanceBound, gaussianRadiusFactor, smoothPow, areaConstant, functionErrorBound, filterRad, bandRadius, defMax, radiusEdgeBound, epsilon, tScale ;
    int manifoldFlag, fixKnn, marchingCubesGridSize, inputType, bandKnn, odepth ;
	bool doOutputRedirection, doCleanResults, useMarchingCubes, contourAtMedian, doSaveFunctionAtVertices ;

	po::options_description options("Splats distance mesher (NCuts)");
	options.add_options()
		( "help", "Produce help message" )
		( "inputFile,i", po::value<std::string>(&inputFile), "Input splats/pts file" )
		( "inputFileType,it", po::value<int>(&inputType)->default_value(0), "Input file type: 0 = splats / 1 = raw points" ) // / 2 = already computed distance function" )
		( "outputFileName,ofn", po::value<std::string>(&outputFileName)->default_value("out.off"), "Output file name" )
		( "outputRunTimesFile,ortf", po::value<std::string>(&outputRunTimesFile)->default_value(""), "Output run times file name (if empty, no log will be made)" )
		( "outputFileDir,od", po::value<std::string>(&outputFileDir)->default_value("./"), "Output mesh directory" ) 
		( "angularBound,ab", po::value<double>(&angularBound)->default_value(0.0), "Mesher: angular bound" )
		( "radiusBound,rb", po::value<double>(&radiusBound)->default_value(0.1), "Mesher: radius bound" )
		( "distanceBound,db", po::value<double>(&distanceBound)->default_value(0.1), "Mesher: distance bound"  )
		( "manifoldFlag,mf", po::value<int>(&manifoldFlag)->default_value(2), "Manifoldness Flag" )
		( "doOutputRedirection,outTxt", po::value<bool>(&doOutputRedirection)->default_value(false), "Output on-screen redirection to txt file" )
		( "doCleanResults,clean", po::value<bool>(&doCleanResults)->default_value(true), "Try some cleaning steps on the output surface" )
		( "gaussianRadiusFactor,gh", po::value<double>(&gaussianRadiusFactor)->default_value(0.2), "Blending Gaussian H factor  (w.r.t. the BSR, if 0 < bandRadius < 1, or w.r.t. AS, if bandRadius > 1)" )
		( "bandRadius,br", po::value<double>(&bandRadius)->default_value(0.4), "Extended band Blending Gaussian H factor  (w.r.t. the BSR, if 0 < bandRadius < 1, or w.r.t. AS, if bandRadius > 1)" )
		( "smoothPow,sp", po::value<double>(&smoothPow)->default_value(4.0), "Smooth prior on the unsigned function" )
		( "areaConstant,ac", po::value<double>(&areaConstant)->default_value(1e-5), "Area constant on the unsigned function." )
		( "functionErrorBound,feb", po::value<double>(&functionErrorBound)->default_value(1e-3), "Maximum allowed error bound between the real implicit function and its linear approximation in the triangulation" )
		( "filterRad,fr", po::value<double>(&filterRad)->default_value(0.5), "Hallucinated triangles filter radius" )
		( "fixKnn,k", po::value<int>(&fixKnn)->default_value(-1), "Fix the K-NN used in distance function computation (-1 = unbounded)" )
		( "bandKnn,bk", po::value<int>(&bandKnn)->default_value(-1), "Fix the K-NN used in the distance function computation for the secondary band (-1 = not used)" )
		( "defMax,dm", po::value<double>(&defMax)->default_value(1000), "Default maximum" )
		( "useMarchingCubes,mc", po::value<bool>(&useMarchingCubes)->default_value(false), "Use Marching Cubes instead of the Surface Mesher" )
		( "marchingCubesGridSize,mcg", po::value<int>(&marchingCubesGridSize)->default_value(500), "Marching Cubes grid size" )
		( "contourAtMedian,cam", po::value<bool>(&contourAtMedian)->default_value(false), "Use median value of implicit function at input points as the contouring isovalue" )
		( "radiusEdgeBound,reb", po::value<double>(&radiusEdgeBound)->default_value(2.5), "Radius-edge bound for the initial coarse implicit function" )
		( "doSaveFunctionAtVertices", po::value<bool>(&doSaveFunctionAtVertices)->default_value(false), "[DEBUG Only] Please, do not use..." )
		( "epsilon,ep", po::value<double>(&epsilon)->default_value(0.001), "Epsilon value used in comparisons" )
		( "T,t", po::value<double>(&tScale)->default_value(2.5), "multiplicative factor of the scale to eliminate hallucinated triangles" )
        ( "octreeDepth,od", po::value<int>(&odepth)->default_value(-1), "Octree depth used to filter the input point set before implicit distance creation. If smaller than one, no octree will be used (default = -1)" )
	;

	// Read parameters from command line
	po::variables_map vm ;
	po::store( po::parse_command_line(argc, argv, options ), vm ) ;
	// Read the configuration file (if exists)
	ifstream ifs( "params.cfg" ) ;
	po::store( po::parse_config_file( ifs, options ), vm ) ;

	po::notify(vm);

	if (vm.count("help")) {
		cout << options << "\n";
		return 1;
	}
	
	// Some input parameters check
	if ( manifoldFlag < 0 || manifoldFlag > 2 ) {
		cerr << "ManifoldnessFlag parameter must be 0, 1 or 2!" << endl ;
		return -1 ;
	}	

	// Get the file type (*.splat)
	string inputFilePathStr( inputFile ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	string inputFileExtension = inputFilePathStr.substr( indexFileExtension+1, ( inputFilePathStr.size()-(indexFileExtension+1) ) ) ;
	if ( ( inputType == 0 && inputFileExtension.compare( "splat" ) != 0 ) ||
		 ( inputType == 1 && ( inputFileExtension.compare( "xyz" ) != 0 && inputFileExtension.compare( "pts" ) != 0 ) ) ) 
	{
		cerr << "[ERROR] Unrecognized input file type <*." << inputFileExtension << ">..." << endl ;
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
		string inputFilePathStr( inputFile ) ; 
		size_t indexFileName = inputFilePathStr.find_last_of( "/\\" ) ;
		size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
		string inputFileNameStr = inputFilePathStr.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;
		
		// Build a name containing the parameters used
		oss << outputFileDir << "/" << inputFileNameStr 
									<< "__ab" << angularBound
									<< "_rb" << radiusBound
									<< "_db" << distanceBound
									<< "_mf" << manifoldFlag 
									<< "_sp" << smoothPow
									<< "_ac" << areaConstant 
									<< "_feb" << functionErrorBound
									<< "_reb" << radiusEdgeBound ;
				
		if ( doOutputRedirection ) {
			otxt << oss.str() << ".txt" ;
		}

		osClean << oss.str() << "_c.off" ;
		osClean2 << oss.str() << "_ch.off" ;

		oss << ".off" ;
		
	}
	else {
		oss << outputFileDir << "/" << outputFileName ;
		osClean << oss.str() << "_c.off" ;
		osClean2 << oss.str() << "_ch.off" ;
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
	cout << "- Parameters:" << endl ;
	cout << "    Input file: " << inputFile << endl ;
	cout << "    Output file: " << fullOutputFilePath << endl ;
	cout << "    Angular Bound = " << angularBound << endl ;
	cout << "    Radius Bound = " << radiusBound << endl ;
	cout << "    Distance Bound = " << distanceBound << endl ;
	cout << "    Manifoldness Flag = " << manifoldFlag << endl ;
	if ( gaussianRadiusFactor > 1 )
		cout << "    Blending Gaussian H (w.r.t AS) = " << gaussianRadiusFactor << endl ;
		if ( bandRadius > gaussianRadiusFactor )
			cout << "    Radius of band (w.r.t AS) = " << bandRadius << endl ;
	else {
		cout << "    Blending Gaussian H (w.r.t BSR) = " << gaussianRadiusFactor << endl ;
		if ( bandRadius > gaussianRadiusFactor )
			cout << "    Radius of band (w.r.t BSR) = " << bandRadius << endl ;
	}
	cout << "    Smooth prior = " << smoothPow << endl ;
	cout << "    Area constant = " << areaConstant << endl ;
	cout << "    Function error bound = " << functionErrorBound << endl ;
	cout << "    Radius-Edge bound = " << radiusEdgeBound << endl ;
	if ( fixKnn > 0 )
		cout << "    Using fixed k-NN = " << fixKnn << endl ;
	if ( bandKnn > 0 )
		cout << "    Secondary band computed using k-NN = " << bandKnn << endl ;
	else
		cout << "    Secondary band default weight = " << defMax << endl ;
	if ( contourAtMedian ) 
		cout << "    Contouring at median" << endl ;
	cout << "    Filter radius (w.r.t. gaussianH) = " << filterRad << endl ;
    if ( odepth > 0 ) {
        cout << "    Using octree of depth = " << odepth << endl ;
    }

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
	
	// Read the input file
	SplatsDistance *distanceFun ;
	if ( inputType == 0 ) {
		/* Splats */
		cout << "- Reading splats from file..." << flush ;
		timer.start() ;
		// Opening file
		ifstream file( inputFile.c_str(), ifstream::in ) ;
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
		// distanceFun = new SplatsDistanceFunctionNCut<K> ( splats, 
		distanceFun = new SplatsDistance ( splats, 
										   gaussianRadiusFactor, 
										   smoothPow,
										   areaConstant, 
										   functionErrorBound,
										   radiusEdgeBound,
										   fixKnn,
										   bandRadius, 
										   defMax,
										   contourAtMedian, 
										   bandKnn,
                                           epsilon,
                                           odepth ) ;
	}
	else {
		/* Raw point set */
		cout << "- Reading points from file..." << flush ;
		timer.start() ;
		// Opening file
		ifstream file( inputFile.c_str(), std::ios_base::in ) ;
		if ( !file ) {
			cerr << "Cannot open input file!" << endl ;	
			return -1 ;
		}
		
		// Read points
		std::vector< Point_3 > points ;
		CGAL::read_xyz_points( file, std::back_inserter(points) ) ;
				
		t = timer.time() ;
		cout << "done (" << t << " s), readed " << points.size() << " points." << endl ;
		if ( logTimings ) otime << "Reading points = " << t << "s" << endl ;
		

		// Define the surface
		cout << "- Computing the signed distance function:" << endl ;
		timer.reset() ;
		// distanceFun = new SplatsDistanceFunctionNCut<K> ( points, 
		distanceFun = new SplatsDistance ( points, 
										   gaussianRadiusFactor, 
										   smoothPow,
										   areaConstant, 
										   functionErrorBound,
										   radiusEdgeBound,
										   fixKnn,
										   bandRadius, 
										   defMax,
										   contourAtMedian,
										   bandKnn, 
										   epsilon ) ;
	}
	t = timer.time() ;
	if ( logTimings ) {
		otime << "Creating SDF = " << t << "s" << endl ;
		otime << "    - Creating distance function = " << distanceFun->getTimeImplicitFunctionCreation() << "s" << endl ;
		otime << "    - Signing distance function = " << distanceFun->getTimeSigningFunction() << "s" << endl ;
	}

	if ( !distanceFun->correct() ) {
		// Some error during implicit function signing...
		std::cout << "[ERROR!] Cannot mesh the implicit function.\n[ERROR!] ABORTING! " << std::endl ;
		return -1 ;
	}

	Surface_3 surface( *distanceFun,										// Pointer to function
                       distanceFun->enlarged_bounding_sphere(1.5) ) ;	// Bounding sphere, where the implicit function is computed

	if ( useMarchingCubes ) {
		cout << "- Meshing the Splats representation (using Marching Cubes)..." << endl ;

		timer.reset() ;
		Polyhedron poly ;
		CGAL::marching_cubes< Surface_3, Polyhedron >(	surface, 
														distanceFun->bounding_box(), 
														marchingCubesGridSize, 
														poly ) ;
		t = timer.time() ;
		if ( logTimings ) otime << "Surface Extraction (MC) = " << t << "s" << endl ;	

		// Save the results
		cout << "- Saving results..." ;
		timer.start() ;
		std::ofstream out( fullOutputFilePath ) ;
		out << poly ;
		out.close() ;
		cout << "done (" << timer.time() << " s)" << endl ;
		
		if ( doCleanResults ) {
			cout << "- Cleaning hallucinated triangles..." << flush ;
			timer.reset() ;
			int numHalluTri = removeHallucinatedTrianglesMC( poly, *distanceFun, filterRad ) ;
			timer.stop() ;
			cout << "done (" << timer.time() << " s), removed " << numHalluTri << " triangles." << endl ;
			timer.reset() ;
			cout << "- Saving clean mesh..." << flush ;
			timer.start() ;
			std::ofstream outClean( osClean2.str().c_str() ) ;
			outClean << poly ;
			outClean.close() ;
			timer.stop() ;
			cout << "done (" << timer.time() << " s)" << endl ;
		}

		if ( logTimings ) {
			if ( logTimings ) otime << "GLOBAL Time = " << globalTimer.time() << "s" << endl ;
			otime.close() ;
		}

		return 0 ;
	}

	cout << "- Meshing the Splats representation (using Surface Mesher)..." << endl ;
	timer.reset() ;

	// Mesh generation
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angularBound,		// Angular bound
														radiusBound,		// Radius bound
														distanceBound );	// Distance bound

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

	timer.reset() ;
	// Save the results
	cout << "- Saving results..." ;		
	timer.reset() ;
	std::ofstream out( fullOutputFilePath ) ;
	CGAL::output_surface_facets_to_off( out, c2t3 ) ;
	out.close() ;

	t = timer.time() ;
	cout << "done (" << t << " s)" << endl ;
	if ( logTimings ) otime << "Saving results = " << t << "s" << endl ;	
	
	//if ( doSaveFunctionAtVertices ) {
	//	typedef C2t3::Vertex_iterator							Vertex_iterator ;
	//	std::ofstream ofav( "./_OUT/fav.xyz" ) ;
	//	Vertex_iterator vit ;
	//	for ( vit = c2t3.vertices_begin(); vit != c2t3.vertices_end(); ++vit ) {
	//		ofav << *vit << " " << (*distanceFun)(vit->point()) << endl ;
	//	}
	//	ofav.close() ;
	//}

	// Perform Post processing (if required)
	if ( doCleanResults ) {
		
		cout << "- Cleaning hallucinated triangles..." << flush ;
		timer.reset() ;
		// int numHalluTri = removeHallucinatedTriangles( c2t3, *distanceFun, filterRad ) ;
		int numHalluTri = removeHallucinatedTriangles2( c2t3, *distanceFun, tScale ) ;
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

	if ( doOutputRedirection )
		outTxt.close() ;


	if ( logTimings ) {
		if ( logTimings ) otime << "GLOBAL Time = " << globalTimer.time() << "s" << endl ;
		otime.close() ;
	}

	return 0 ;
}



/* Post-processing (implementation) */
//int removeHallucinatedTriangles( C2t3 &mesh, const SplatsDistance &func, const FT &radiusPart ) {
//
//	FT rad = radiusPart ;
//	if ( rad > 1 || rad < 0 ) {
//		std::cout << "[WARNING] Hallucinated triangles radius must be between 0 and 1! Defaulting to 0.5..." << std::endl ;
//		rad = 0.5 ;
//	}
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
//
//		}
//	}
//
//	return numRemoved ;
//}


/* Post-processing (implementation) */
int removeHallucinatedTrianglesMC( Polyhedron& poly, const SplatsDistance &func, const FT &radiusPart ) {

	typedef Polyhedron::Facet_iterator		Facet_iterator ;

	FT rad = radiusPart ;
	if ( rad > 1 || rad < 0 ) {
		std::cout << "[WARNING] Hallucinated triangles radius must be between 0 and 1! Defaulting to 0.5..." << std::endl ;
		rad = 0.5 ;
	}

	int numRemoved = 0 ;
	for( Facet_iterator f = poly.facets_begin(); 
		 f != poly.facets_end(); ++f ) 
	{
		//get the three vertices
		Triangle_3 t(	f->halfedge()->vertex()->point(),
						f->halfedge()->next()->vertex()->point(),
						f->halfedge()->prev()->vertex()->point() ) ;

		if ( !func.allPointsOnRestrictedBand( t, rad ) )
		{
			poly.erase_facet( f->halfedge() ) ;
			numRemoved++ ;
		}
	} 

	return numRemoved ;

}

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


//bool readDiscretizedDistanceFunction( std::istream& stream, std::vector< Point_3 >& vertices, std::vector< double > funcAtVert ) {
//	while( stream.good() ) {
//		FT x, y, z, f ;
//		stream	>> x 
//				>> y 
//				>> z 
//				>> f ;
//		vertices.push_back( Point_3( x, y, z ) ) ;
//		funcAtVert.push_back( f ) ;
//	}
//}
