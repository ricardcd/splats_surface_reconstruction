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
//#define cimg_display 0
//#include "CImg.h"
// Std includes
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
// Boost
#include <boost/program_options.hpp>
// Project related
#include "PointSetReader/PointSetReader.h" // Reader for pointsets (oriented or not)
#include "PointSetSurfaceMesher/PointSetSurfaceMesherOracle.h"
#include "C2t3/Complex_2_in_triangulation_3_post_processing.h"
#include "PointSetSurfaceMesher/SegmentQueryIntersectionOracle.h"

using namespace std ;
namespace po = boost::program_options;

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
	std::string inputFile, outputFile, outputDir;
	double quantile, probNoOutliers, expectFracOutliers, ransacThres, ransacT, capRadPer, weightHFactor, scaleRadFactor, intConsFact, angularBound, radiusBound, distanceBound;
	int degree, scalesEstimator, scaleK, minInliers, intConsNum, manifoldFlag;
	unsigned int ransacThresholdType;
	bool forceLBQOnSegment, computeDensity, doCleanResults, doOutputTimes, mesherDebugOutput;

    po::options_description options("Direct meshing of pointsets." ) ;
    options.add_options()
            ("help,h", "Produce help message")
	        // Basic
            ("inFile,i", po::value<std::string>(&inputFile), "Input point set file path.")
            ("outFile,o", po::value<std::string>(&outputFile)->default_value(""), "Output mesh file path (OFF format).")
            ("outDir", po::value<std::string>(&outputDir)->default_value("."), "Output directory. If outFile is not specified, a file with a default name containig a list of the parameters used will be written in this directory.")
            ("deg", po::value<int>(&degree)->default_value(2), "Local surface degree. Available: 1 (plane) and 2 (Local Bivariate Quadric).")
            ("forceSegment", po::value<bool>(&forceLBQOnSegment)->default_value(true), "Force the Local Bivariate Quadric to follow the query segment.")
	        // Scale
            ("computeScale", po::value<unsigned int>(&ransacThresholdType)->default_value(0), "Flag indicating how the scale of the noise has to be computed using the Modified Selective Statistical Estimator (MSSE). Available options: 0 = Not computed (use a fixed RANSAC threshold) / 1 = local computation (at each iteration) / 2 = precompute scales before meshing using a radial neighborhood / 3 = precompute scales before meshing using a fixed k-neighborhood / 4 = Load precomputed scales from file" )
            ("quantile", po::value<double>(&quantile)->default_value(0.1), "[MSSE] Minimum inlier points expected inside the neighborhood (percentage, starting point to compute scale through MSSE).")
            ("probNoOutliers", po::value<double>(&probNoOutliers)->default_value(0.99), "[MSSE] Probability of getting a sample free from outliers in the Least Kth Squares (LKS) procedure.")
            ("expectFracOutliers", po::value<double>(&expectFracOutliers)->default_value(0.7), "[MSSE] Expected fraction of outliers for the LKS procedure.")
            ("scalesEst", po::value<int>(&scalesEstimator)->default_value(1), "[MSSE] Scales model estimator (0 = Plane, 1 = LBQ).")
            ("scaleK", po::value<int>(&scaleK)->default_value(200), "[MSSE, if computeScale=3] Number of neighbors to take into account during scale computation using k-nearest neighbors")
            ("scaleRadFactor", po::value<double>(&scaleRadFactor)->default_value(3.0), "[MSSE, if computeScale=2] Defines the search radius for computing MSSE based on the capsule radius (--capRad). The search radius will be capRad * scaleRadFactor.")
            // Ransac
            ("fixedRansacThres", po::value<double>(&ransacThres)->default_value(0.7), "Fixed RANSAC threshold (used when computeScale=0).")
            ("ransacThresFactor", po::value<double>(&ransacT)->default_value(2.5), "Factor to multiply the scale in order to determine the RANSAC threshold (used when computeScale>0")
            ("minInliers", po::value<int>(&minInliers)->default_value(15), "Minimum number of inliers for the RANSAC local surface fitting to consider it as valid.")
	        // Capsule search
            ("capRad", po::value<double>(&capRadPer)->default_value(0.1), "Radius for query capsule's segment, proportional to the minimum bounding sphere radius.")
            ("weightHFactor", po::value<double>(&weightHFactor)->default_value(1.0), "H constant of the Gaussian weighting used during LocalBivariateQuadric fitting, proportional to the search capsule radius.")
            // Density
            ("densityAdaptiveCapRad", po::value<bool>(&computeDensity)->default_value(false), "Use a density-adaptive capsule query radius.")
	        // Intersection validation
            ("consensusFact", po::value<double>(&intConsFact)->default_value(0.5), "Factor defining the proportion of the current radius where to search for points in consensus with the intersection, that is, falling inside a sphere search of (current radius)*consesusFact")
            ("consensusNum", po::value<int>(&intConsNum)->default_value(3), "Number of points falling at a distance smaller than (current radius)*consesusFact in order to consider the intersection as valid.")
            // Meshing
            ("ab", po::value<double>(&angularBound)->default_value(10.0), "[Surface Mesher] Angular Bound.")
            ("rb", po::value<double>(&radiusBound)->default_value(0.1), "[Surface Mesher] Radius Bound.")
            ("db", po::value<double>(&distanceBound)->default_value(0.1), "[Surface Mesher] Distance Bound.")
            ("mf", po::value<int>(&manifoldFlag)->default_value(2), "[Surface Mesher] Manifoldness Flag.")
            ("clean", po::value<bool>(&doCleanResults)->default_value(false), "[Surface Mesher] Try some cleaning steps on the output surface.")
            // Debug
            ("outTimes", po::value<bool>(&doOutputTimes)->default_value(true), "Output execution times to txt file.")
            ("debugOut", po::value<bool>(&mesherDebugOutput)->default_value(false), "Activate debug output.")
        ;

    // Read parameters from command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << options << "\n";
        return 1;
    }
    if (outputFile.empty() && outputDir.empty()) {
        std::cout << "[ERROR] Please set either --outFile or --outDir parameters (or both)." << std::endl;
        return -1;
    }
    if (degree < 1 || degree > 2 ) {
        std::cerr << "[ERROR] --deg can only be either 1 (a plane) or 2 (a bivariate quadric)." << std::endl;
        return -1;
    }
    if (ransacThresholdType < 0 || ransacThresholdType > 4) {
        std::cerr << "[ERROR] --computeScale allowed values are 0, 1, 2, 3 or 4 (see help)" << std::endl;
        return -1;
    }

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
														  mesherDebugOutput);
														  //mesherDebugVisualize ) ;
		
	// Create the name of the output file
	std::ostringstream oss;
	std::ostringstream otimes ; // Debug only!	
	std::ostringstream osClean ; // Debug only!
	string inputFileNameStr ;
	string inputFilePathStr( inputFile ) ;
	size_t indexFileName = inputFilePathStr.find_last_of( "/\\" ) ;
	size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
	inputFileNameStr = inputFilePathStr.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;
	if ( outputFile.empty() ) {
		oss << outputDir << "/" << inputFileNameStr ;
		// Related to local surface computation
		oss << "_" << intOraclePtr->getParametersString() ;
		// Related to meshing
		oss << "__ab" << angularBound ;
		oss << "_db" << distanceBound ;
		oss << "_rb" << radiusBound ;
		// Related to manifoldness
		oss << "__mf" << manifoldFlag ;

		otimes << oss.str() << "__RunTimesLog.txt" ;
		osClean << oss.str() << "__clean.off" ;
		oss << ".off" ;		
	}	
	else {
		oss << outputFile ;
	}
	cout << "- Output file path: " << oss.str() << endl ;

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
