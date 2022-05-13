/* Includes */
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/make_surface_mesh.h>
	// Timer (debugging)
#include <CGAL/Timer.h>
	// Surface Mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
	// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h> 
#include <CGAL/IO/read_xyz_points.h>
// CImg include (only used here for easily dealing with input parameters)
#include "CImg.h"
// Std includes
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
// Project related
#include "PointSetReader/PointSetReader.h" // Reader for pointsets (oriented or not)
#include "PointSetSurfaceMesher/PointSetSurfaceMesherOracle.h"
#include "C2t3/Complex_2_in_triangulation_3_post_processing.h"
#include "PointSetSurfaceMesher/SegmentQueryIntersectionOracle.h"

using namespace std ;

// Main Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

// Default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3		Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr>			C2t3;
typedef CGAL::Simple_cartesian<double>					Simple_cartesian_kernel;
typedef Tr::Triangulation_data_structure				Tds ;
typedef C2t3::Vertex_iterator							Vertex_iterator ;
typedef C2t3::Edge_iterator								Edge_iterator ;
typedef Tr::Cell_handle									Cell_handle ;
typedef Tr::Facet										Facet ;
typedef Tr::Cell_circulator								Cell_circulator ;
typedef Tr::Facet_circulator							Facet_circulator ;

// Custom oracle
typedef CGAL::PointSetSurfaceMesherOracle<K>			PointsetSurfaceMesherOracle ;



/* Main function */
int main ( int argc, char **argv) {

	/* Parse input parameters */
	cimg_usage( "Direct meshing of pointsets." ) ;
	// Basic
	const char*			inputFile = cimg_option(			"-i",					(char*)0,	"Input Point Set file." ) ;
	const char*			outputFile = cimg_option(			"-o",					(char*)0,	"Output file." ) ;
	const char*			outputDir = cimg_option(			"-od",					(char*)0,	"Output directory." ) ;
	const int			degree = cimg_option(				"-deg",					2,			"Local surface degree. Available: 1 (plane) and 2 (Local Bivariate Quadric)." ) ;
	const bool			forceLBQOnSegment = cimg_option(	"-forceSegment",		true,		"Force the Local Bivariate Quadric to follow the query segment." ) ;
	// Scale
	const unsigned int	ransacThresholdType=cimg_option(	"-computeScale",		0,			"Flag indicating how the scale of the noise has to be computed: 0 = Not computed (fixed RANSAC threshold) / 1 = local computation (at each iteration) / 2 = precompute scales before meshing / 3 = Use fixed k-neighborhood / 4 = Load precomputed scales from file" ) ;
	const double		quantile = cimg_option(				"-quantile",			0.1,		"Minimum inlier points expected inside the neighborhood (percentage, starting point to compute scale through MSSE)." ) ;
	const double		probNoOutliers = cimg_option(		"-probNoOutliers",		0.99,		"Probability of getting a sample free from outliers in the LKS procedure." ) ;
	const double		expectFracOutliers=cimg_option(		"-expectFracOutliers",	0.7,		"Expected fraction of outliers for the LKS procedure." ) ;	
	const int			scalesEstimator = cimg_option(		"-scalesEst",			1,			"Scales model estimator (0 = Plane, 1 = LBQ)." ) ;
	const int			scaleK = cimg_option(				"-scaleK",				200,		"Number of neighbors to take into account during scale computation using k-nearest neighbors (used if -computeScale = 3)" ) ;
	// Ransac
	const double		ransacThres = cimg_option(			"-rt",					0.7,		"RANSAC threshold." ) ;
	const double		ransacT	=cimg_option(				"-ransacT",				2.5 ,		"Factor to multiply the scale in order to determine the RANSAC threshold" ) ;
	const int			minInliers = cimg_option(			"-minInliers",			15,			"Minimum number of inliers for the ransac plane fitting in order to take into account this splat." ) ;
	// Capsule search
	const double		capRadPer = cimg_option(			"-capRad",				0.1,		"Radius for query capsule's segment, proportional to the minimum bounding sphere radius." ) ;
	const double		weightHFactor = cimg_option(		"-weightHFactor",		1.0,		"H constant of the Gaussian weighting used during LocalBivariateQuadric fitting, proportional to the search capsule radius." ) ;
	const double		scaleRadFactor=cimg_option(			"-scaleRadFactor",		3.0 ,		"This factor will multiply the query capsule radius (scale computation may need larger neighborhoods)." ) ;
	// Density	
	const bool			computeDensity = cimg_option(		"-computeDensity",		false,		"Density-aware query radius adaptation." ) ;
	// Intersection validation
	const double		intConsFact	= cimg_option(			"-consensusFact",		0.5 ,		"Factor defining the proportion of the current radius where to search for points in consensus with the intersection, that is, falling inside a sphere search of (current radius)*consesusFact" ) ;
	const int			intConsNum	= cimg_option(			"-consensusNum",		3 ,			"Number of points falling at a distance smaller than (current radius)*consesusFact in order to consider the intersection as valid." ) ;
	// Meshing
	const double		angularBound = cimg_option(			"-ab",					10.0,		"Angular Bound." ) ;
	const double		radiusBound = cimg_option(			"-rb",					0.1,		"Radius Bound." ) ;
	const double		distanceBound = cimg_option(		"-db",					0.1,		"Distance Bound." ) ;
	const int			manifoldFlag = cimg_option(			"-mf",					2,			"Manifoldness Flag." ) ;
	const bool		    doCleanResults = cimg_option(		"-clean",				false,		"Try some cleaning steps on the output surface." ) ;	
	// Debug
	const bool			doOutputRedirect = cimg_option(		"-outTxt",				false,		"Output on-screen redirection to txt file." ) ;	
	const bool			doOutputTimes = cimg_option(		"-outTimes",			true,		"Output execution times to txt file." ) ;	
	const bool		    mesherDebugOutput = cimg_option(	"-debugOut",			false,		"Debug output." ) ;	
	const bool		    mesherDebugVisualize = cimg_option(	"-debugVis",			false,		"Debug visualization." ) ;	
	

	// Initialize some variables...
	CGAL::Timer timer ; // Reusable timer initialization
	double logTimes[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;
	SegmentQueryIntersectionOracle<K>* intOraclePtr ;	// Main SplatsCreator object
	bool readed = false ;
	
	// Read the input pointset
	vector< Point_3 > points ; // Input points
	cout << "- Reading input pointset file..." << flush ;	
	timer.start() ;
	std::ifstream in( inputFile ) ; 
	if ( !in ) {
		std::cerr << "[Error] Cannot open the input file..." << std::endl ;
	}
	readed = read_xyz_points( in, std::back_inserter( points ) ) ;
	if ( !readed ) {
		cerr << "[Error] Cannot read the input file..." << endl ;
		return -1 ;
	}
	if ( points.empty() ) {
		cerr << "[Error] There are no points in the input file..." << endl ;
		return -1 ;
	}
	timer.stop() ;
	logTimes[0] = timer.time() ;
	cout << "done (" << timer.time() << " s)" << endl ;
	timer.reset() ;
	in.close() ;

	/* Set the SplatsCreator parameters */
	intOraclePtr = new SegmentQueryIntersectionOracle<K>( points, 
														  degree,
														  forceLBQOnSegment,
														  ransacThres,
														  minInliers,
														  quantile,
														  probNoOutliers,
														  expectFracOutliers,
														  ransacThresholdType,
														  scaleK,
														  capRadPer,
														  weightHFactor,
														  scaleRadFactor, 
														  scalesEstimator,
														  computeDensity,
														  ransacT,
														  intConsFact, 
														  intConsNum,
														  mesherDebugOutput, 
														  mesherDebugVisualize ) ;
		
	// Create the name of the output file
	std::string str ;
	std::ostringstream oss( str ) ;
	std::ostringstream otxt ; // Debug only!	
	std::ostringstream otimes ; // Debug only!	
	std::ostringstream osClean ; // Debug only!
	string inputFileNameStr ;
	string inputFilePathStr( inputFile ) ;
	size_t indexFileName = inputFilePathStr.find_last_of( "/\\" ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	inputFileNameStr = inputFilePathStr.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;
	if ( outputFile == (char*)0 ) {
		oss << outputDir << "/" << inputFileNameStr ;
		// Related to local surface computation
		oss << "_" << intOraclePtr->getParametersString() ;
		// Related to meshing
		oss << "__ab" << angularBound ;
		oss << "_db" << distanceBound ;
		oss << "_rb" << radiusBound ;
		// Related to manifoldness
		oss << "__mf" << manifoldFlag ;

		otxt << oss.str() << "__coutRedirect.txt" ;
		otimes << oss.str() << "__RunTimesLog.txt" ;
		osClean << oss.str() << "__clean.off" ;
		oss << ".off" ;		
	}	
	else {
		oss << outputDir << "/" << outputFile ;		
	}
	cout << "- Output file path: " << oss.str() << endl ;

	// Redirect the output if needed
	std::ofstream outTxt;
	if ( doOutputRedirect ) {
		cout << "- Redirecting output to a txt file..." << endl ;
		outTxt.open( otxt.str().c_str() ) ;
		std::streambuf *coutbuf = std::cout.rdbuf() ; //save old buf
		std::cout.rdbuf( outTxt.rdbuf() ) ;
	}

	// Show useful info on screen
	intOraclePtr->debugShowParameters() ;

	// Prepares the Segment Query Intersection Oracle for the queries
	std::string prepStr ;
	//if ( ransacThresholdType > 3 ) {
	//	// File path points to the precomputed scales file
	//	prepStr = inputScalesFile ;
	//}
	//else {
		prepStr = inputFileNameStr ;
	// }
	timer.start() ;
	intOraclePtr->prepare( outputDir, prepStr ) ;
	timer.stop() ;
	logTimes[1] = timer.time() ;
	std::cout << "- Scales computed/loaded in " << logTimes[1] << " s" << std::endl ;
	timer.reset() ;


	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // Initialize 2D-complex in 3D-Delaunay triangulation	

	// Define meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria( angularBound,		// angular bound
														radiusBound,		// radius bound
														distanceBound ) ;	// distance bound

	// Mesh!
	PointsetSurfaceMesherOracle oracle( intOraclePtr ) ;

	cout << "- Meshing the Point Set..." << flush ;		
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
	logTimes[2] = timer.time() ;
	timer.reset() ;

	// Save the results
	cout << "- Saving results (" << oss.str() << ")..." << flush ;		
	timer.start() ;
	std::ofstream outOr( oss.str().c_str() ) ; 
	CGAL::output_surface_facets_to_off( outOr, c2t3 ) ; 
	outOr.close() ;
	timer.stop() ;
	cout << "done (" << timer.time() << " s)" << endl ;	
	logTimes[3] = timer.time() ;
	timer.reset() ;

	// Perform Post processing (if required)
	if ( doCleanResults ) {
		int numNonManifoldEdges = 0, numNonManifoldVertices = 0 ;
		if ( !CGAL::C2t3_PostProcessing::manifoldnessCheck< C2t3 >( c2t3, numNonManifoldEdges, numNonManifoldVertices ) ) {

			cout << "- Cleaning resulting mesh..." << flush ;		
			timer.start() ;
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

			timer.stop() ;
			cout << "done (" << timer.time() << " s)" << endl ;
			logTimes[4] = timer.time() ;
			timer.reset() ;

			cout << "- Saving clean mesh..." << flush ;		
			timer.start() ;
			std::ofstream outClean( osClean.str().c_str() ) ; 
			CGAL::output_surface_facets_to_off( outClean, c2t3 ) ; 
			outClean.close() ;
			timer.stop() ;
			cout << "done (" << timer.time() << " s)" << endl ;	
			logTimes[5] = timer.time() ;
		}
		else {
			std::cout << "- Resulting surface is manifold, no cleaning needed..." << std::endl ;
		}
	}	


	// Log running times
	if ( doOutputTimes ) {
		std::ofstream outTimes( otimes.str().c_str() ) ;
		outTimes << "- Read input pointsets = " << logTimes[0] << "s" << std::endl ;
		if ( ransacThresholdType > 0 )
			outTimes << "- Scale computation = " << logTimes[1] << "s" << std::endl ;
		outTimes << "- Point set meshing = " << logTimes[2] << "s" << std::endl ;
		outTimes << "- Saving raw mesh = " << logTimes[3] << "s" << std::endl ;
		if ( doCleanResults ) {
			outTimes << "- Cleaning mesh = " << logTimes[4] << "s" << std::endl ;
			outTimes << "- Saving cleaned mesh = " << logTimes[5] << "s" << std::endl ;
		}
	}
		
}
