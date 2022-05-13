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
#include "CImg.h"

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
	cimg_usage( "Compute splats for pointsets." ) ;

	const char*			inputFile = cimg_option( "-i",				(char*)0, "Input Point Set file." ) ;
	const char*			outputFile = cimg_option( "-o",				(char*)0, "Output file." ) ;
	const char*			outputDir = cimg_option( "-od",				(char*)0, "Output file." ) ;
	const int			inputType = cimg_option( "-inputType",		0, "Input type: 0 = Pointset / 1 = Oriented Pointset / 2 = Oriented Pointset With Scores." ) ;
	const int			neighType = cimg_option( "-neighType",		0, "Neighborhood type: 0 = k-NN (default) / 1 = Radial k-NN." ) ;
	const double		neighRadius = cimg_option( "-radius",		numeric_limits<double>::infinity(), "Radial k-NN distance (default = Infinity)." ) ;
	const int			distType = cimg_option( "-distanceType",	0, "Distance type: 0 = k-NN (3D) / 1 = k-NN (2D) / 2 = Radius to enclosing 2D circle." ) ;
	const int			k = cimg_option( "-k",						0, "Number of nearest neighbors to take into account." ) ;
	const double		multFact = cimg_option( "-multFactor",		1.0, "Multiplicative factor" ) ;	
	const unsigned int	lssDegree = cimg_option( "-d",				2, "local smooth surface degree." ) ;
	const double		maxSize = cimg_option( "-maxDiscSize",		numeric_limits<double>::infinity(), "Maximum disc size." ) ;
	const int			minInliers = cimg_option( "-minInliers",	15, "Minimum number of inliers for the ransac plane fitting in order to take into account this splat." ) ;
	const int			minFitPts = cimg_option( "-minFitPts",		6, "Minimum number of points to generate a splat approximation." ) ;
	const bool			saveTime = cimg_option( "-saveTime",		true, "Flag indicating wether the time has to be saved in a separated file (same name as output file, but ended with __time.txt)." ) ;	
	const bool			preComputeScale = cimg_option( "-computeScale", true, "Flag indicating wether the scale of the noise has to be computed" ) ;			
	const int           scaleK = cimg_option( "-scaleK",    500, "Neighborhood to compute the scale from." ) ;			
	const double		quantile = cimg_option( "-quantile",		0.3, "Maximum deviation from computed local smooth surface." ) ;
	const double		p = cimg_option( "-p",						0.99, "Probability of getting a sample free from outliers in the LKS procedure." ) ;
	const double		epsilon = cimg_option( "-epsilon",			0.7, "Expected percentage of outliers for the LKS procedure." ) ;
	const double		ransacThres = cimg_option( "-rt",			0.7, "RANSAC threshold." ) ;

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
		readed = PointSetReader< Point_3, Vector_3 >::readPointsetFile( inputFile, points ) ;
		if ( readed ) {
			sc = SplatsCreator( points ) ;
		}
		else {
			cerr << "Error reading the file..." << endl ;
		}
		break ;
	case 1:
		// Oriented pointsets
		readed = PointSetReader< Point_3, Vector_3 >::readOrientedPointsetFile( inputFile, points, normals ) ;
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
	sc.setDistType( distType ) ;
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
	if ( outputFile == (char*)0 ) {
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
		if ( outputDir == (char*)0 ) {
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

