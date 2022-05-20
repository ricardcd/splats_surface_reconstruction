/* Includes */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// Timer (debugging)
#include <CGAL/Timer.h>
#include <CGAL/eigen.h>
#include <CGAL/PCA_util.h>
// Reader for pointsets (oriented or not)
#include "PointSetReader/PointSetReader.h"
// Splats creator
#include "Splat_3/Splat_3.h"
#include "SplatsCreation/SplatsCreator.h"
#include "Splat_3/SplatsIO.h"
// Other
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <boost/program_options.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT										FT;
typedef K::Point_3									Point_3;
typedef K::Vector_3									Vector_3;
typedef Splats::SplatsCreator<K>					SplatsCreator;
typedef CGAL::Splat_3<K>							Splat_3;

using namespace std;
namespace po = boost::program_options;

/* Main function */
int main ( int argc, char **argv) {

	/* Parse input parameters */
	std::string inputFile, outputFile, outputDir;
	int inputType, neighType, distType, k, minInliers, minFitPts, scaleK;
	double neighRadius, multFact, maxSize, quantile, p, epsilon, ransacThres;
	unsigned int lssDegree;
	bool saveTime = false, preComputeScale = false;

	po::options_description options("Computes splats from a point sets");
	options.add_options()
            ("help,h", "Produce help message")
			("inFile,i", po::value<std::string>(&inputFile), "Input point set file path.")
			("outFile,o", po::value<std::string>(&outputFile), "Output mesh file path (OFF format).")
			("outDir", po::value<std::string>(&outputDir), "Output directory. If outFile is not specified, a file with a default name containig a list of the parameters used will be written in this directory.")
			("inputType", po::value<int>(&inputType)->default_value(0), "Input type: 0 = Pointset / 1 = Oriented Pointset / 2 = Oriented Pointset With Scores.")
			("neighType", po::value<int>(&neighType)->default_value(0), "Neighborhood type: 0 = k-Nearest Neighbors (default) / 1 = Radially-Nearest Neighbors.")
			("radius", po::value<double>(&neighRadius)->default_value(numeric_limits<double>::infinity()), "Radial Nearest Neighbors distance (only used if neighType == 1)")			
			("knn,k", po::value<int>(&k)->default_value(25), "Number of nearest neighbors to take into account (only used if neighType == 0)")
            ("lssDegree,d", po::value<unsigned int>(&lssDegree)->default_value(2), "Degree of the local smooth surface computed for each splat.")
            ("maxDiscSize", po::value<double>(&maxSize)->default_value(numeric_limits<double>::infinity()), "Maximum size (i.e., radial coverage) of a splat.")
            ("minInliers", po::value<int>(&minInliers)->default_value(15), "Minimum number of inliers for the RANSAC plane fitting (the local reference plane of the splat).")
            ("ransacThres,t", po::value<double>(&ransacThres)->default_value(0.7), "RANSAC threshold")
            ("minFitPts", po::value<int>(&minFitPts)->default_value(6), "Minimum number of fitting points required to generate a splat (corresponds to the minimum number of points needed by RANSAC to fit a local surface of the selected degree).")
            ("saveTime", po::bool_switch(&saveTime), "If set, the run time will be saved in a separate file.")
            ("preComputeScale", po::bool_switch(&preComputeScale), "If set, the scale of the noise will be estimated")
            ("scaleK", po::value<int>(&scaleK)->default_value(500), "Number of nearest neighbors used to compute the scale (only used if --preComputeScale is set).")
            ("scaleQuantile", po::value<double>(&quantile)->default_value(0.3), "Quantile of the Least Kth Squares (LKS) procedure used to compute the scale of the noise (only used if --preComputeScale is set).")
            ("scaleProbNoOutliers", po::value<double>(&p)->default_value(0.99), "Probability of getting a sample free of outliers for the LKS for computing the scale of the noise (only used if --preComputeScale is set).")
            ("scaleExpectedOutliersPercent", po::value<double>(&epsilon)->default_value(0.7), "Expected percentage of outliers for the LKS for computing the scale of the noise (only used if --preComputeScale is set).")
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

	CGAL::Timer timer ; // Reusable timer initialization
	CGAL::Timer globalTimer ; // Global timer
	globalTimer.start() ;

	/* Read the pointset */
	vector< Point_3 > points ;
	vector< Vector_3 > normals ;
	vector< double > scores ;
	bool readed = false ;

	cout << "- Reading input pointset file and creating the Splats Creator structures..." << flush ;	
	timer.start() ;

	SplatsCreator sc ;
	switch( inputType ) {
	case 0:
		// Simple pointsets
		readed = PointSetReader< Point_3, Vector_3 >::readPointsetFile( inputFile.c_str(), points ) ;
		if ( readed ) {
			sc = SplatsCreator( points ) ;
		}
		else {
			cerr << "Error reading the file..." << endl ;
		}
		break ;
	case 1:
		// Oriented pointsets
		readed = PointSetReader< Point_3, Vector_3 >::readOrientedPointsetFile( inputFile.c_str(), points, normals ) ;
		if ( readed ) {
			sc = SplatsCreator( points, normals ) ;
		}
		else {
			cerr << "Error reading the file..." << endl ;
		}
		break ;
	default:
		cerr << "Unknown input type!" << endl ;
		return -1 ;
	}
	
	if ( !readed ) {
		cerr << "Unable to parse the input file!" << endl ;
		return -1 ;
	}
	timer.stop() ;
	cout << "done (" << timer.time() << " s)" << endl ;
	timer.reset() ;

	/* Set Parameters */
	sc.setNeighType( neighType ) ;
	sc.setNeighRadius( neighRadius ) ;
	// sc.setDistType( distType ) ;
	sc.setK( k ) ;
	sc.setScaleK( scaleK ) ;
	sc.setLssDegree( lssDegree ) ;
	sc.setMinInliers( minInliers ) ;
	sc.setMinFitPts( minFitPts ) ;
	sc.setDoComputeScale( preComputeScale ) ; 
	sc.setQuantile( quantile ) ;
	sc.setProbNoOutliers( p ) ;
	sc.setExpectFracOutliers( epsilon ) ;
	sc.setFixedDistanceThreshold( ransacThres ) ;

	// Create the name of the output file
	std::string str ;
	std::ostringstream oss( str ) ;
	std::ostringstream otxt ; // Debug only!	
	std::string baseFileName ;
	if ( outputFile.empty() ) {
		string inputFilePathStr( inputFile ) ;
		size_t indexFileName = inputFilePathStr.find_last_of( "/\\" ) ;
		size_t indexFileExtension = inputFilePathStr.find_last_of( "." ) ;
		string inputFileNameStr = inputFilePathStr.substr( indexFileName+1, ( indexFileExtension - indexFileName - 1 ) ) ;

		oss << outputDir << "/" << inputFileNameStr ;
		oss << "__" << sc.getParametersString() ;

		otxt << oss.str() << "__time.txt" ;
		baseFileName = oss.str() ; // Save the base name of the output file
		oss << ".splat" ;		
	}	
	else {
		if ( outputDir.empty() ) {
			oss << outputFile ;
		}
		else {
			oss << outputDir << "/" << outputFile ;		
		}
	}

	cout << "- Output file path: " << oss.str() << endl ;

	if ( preComputeScale ) {
		if ( !sc.loadPrecomputedScales( outputDir, baseFileName ) ) {
			cout << "- Pre-computing scales..." ;
			timer.start() ;
			std::vector< Splat_3 > splats ;
			sc.precomputeScales( true ) ;
			timer.stop() ;
			cout << "done" << std::endl ;
			cout << "- Saving Pre-computing scales..." ;
			if ( !sc.savePrecomputedScales( outputDir, baseFileName ) ) {
				std::cout << "A problem occurred, and scales file has not been saved..." << std::endl ;
			}
			else {
				cout << "done" << std::endl ;
			}
		}
		else {
			cout << "- Pre-computed scales loaded from file." ;
		}
	}

	cout << "- Computing splats..." ;
	timer.start() ;
	std::vector< Splat_3 > splats ;
	sc.run( splats, true ) ;
	timer.stop() ;
	cout << "done" << std::endl ;

	// Save results
	cout << "- Writing results to file..." << flush ;
	ofstream outFile ;
	outFile.open( oss.str().c_str() ) ;
	if ( !outFile.is_open() ) {
		cerr << "Cannot open the output file for writing!" << endl ;
		return -1 ;
	}
	SplatsIO::WriteSplats<K>( outFile, splats ) ;
	cout << "done" << std::endl ;

	// Save time (if required)
	if ( saveTime ) {
		std::cout << "- Saving time measure..." << std::endl ;
		ofstream outFileTime ;
		outFileTime.open( otxt.str().c_str() ) ;
		if ( !outFileTime.is_open() ) {
			cerr << "Cannot open the output file for writing time!" << endl ;
			return -1 ;
		}
		outFileTime << "Total time = " << timer.time() << " s" << std::endl ;
		outFileTime.close() ;
	}


}

