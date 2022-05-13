#include "Models/FitLine2D/LineModelEstimator.h"
#include "Ransac.h"
#include "LeastKthSquares.h"
#include "ScaleEstimation.h"
#include <random>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2		Point_2 ;

int  main(int argc, char* argv[])
{
	std::vector<double> lineParameters ;
	double distanceThreshold = 0.5 ; //for a point to be on the line it has to be closer than 0.5 units from the line
	LineModelEstimator * lpEstimator = new LineModelEstimator( distanceThreshold ) ;
	std::vector< Point_2 > pointData ;	
	int numForEstimate = 2 ;
	int numSamples = 10000 ;
	int numOutliers = 0 ;
	double desiredProbabilityForNoOutliers = 0.99 ;
	double quantile = 0.3 ;
	double expectedFractionOfOutliers = 0.6 ;	
	double noiseSpreadRadius = 0.03 ;
	double noiseStd = 0.3 ;
	double outlierSpreadRadius = 10 ;
	int i ;
	double newX, newY, dx, dy, norm ;

    //1.Create data with outliers

	//randomly select a direction [dx,dy] and create a line passing through the origin
	//for each point sampled on the line add random noise, finally add outlying 
	//points in the direction of the line normal.

	srand( (unsigned)time( NULL ) ) ; //seed random number generator

	//get random direction
	dx = rand() ;
	dy = rand() ;
	norm = sqrt(dx*dx + dy*dy) ;
	dx/= norm ;
	dy/= norm ;
	dx *= (rand() > RAND_MAX/2 ? 1 : -1) ;
	dy *= (rand() > RAND_MAX/2 ? 1 : -1) ;


	// Add 'numSamples' noisy points 
	std::tr1::ranlux64_base_01 eng ; // Random number generator engine
    std::tr1::normal_distribution<double> normalDist( 0.0, noiseStd ) ; // Specific distribution	
	for( i = 0; i < numSamples; i++ ) {
		// Using uniform noise:
		// newX = i*dx + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1) ;
		// newY = i*dy + noiseSpreadRadius*(double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1) ;

		// Using normal(gaussian) noise
		newX = i*dx + normalDist(eng) ;
		newY = i*dy + normalDist(eng) ;
				
		pointData.push_back( Point_2(newX,newY) ) ;
	}

	//'numOutliers' points
	double centerX = -dy*100;
	double centerY = dx*100;
	for(i=0; i<numOutliers; i++) {
		newX = centerX + outlierSpreadRadius * (double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		newY = centerY + outlierSpreadRadius * (double)rand()/(double)RAND_MAX * (rand() > RAND_MAX/2 ? 1 : -1);
		
		pointData.push_back( Point_2(newX,newY) ) ;
	}
	
	double dotProd;
                        
	//2. Compare least squares approach to Ransac

	std::cout<<"Total number of points used: "<<pointData.size()<<std::endl;
	std::cout<<"Number of outliers: "<<numOutliers<<std::endl;
	// The real line parameters
	std::cout<<"Real line parameters [nx,ny,ax,ay]\n\t [ "<<-dy<<", "<<dx<<", 0, 0 ]"<<std::endl;

	// A least squares estimate of the line parameters
	// Use pointers to all points in the data vector
	std::vector< std::vector< Point_2 >::const_pointer > pointDataPtr ;
	std::vector< Point_2 >::iterator itPts ;
	for ( itPts = pointData.begin(); itPts != pointData.end(); ++itPts ) {
		pointDataPtr.push_back( &(*itPts) ) ;
	}

	lpEstimator->leastSquaresFit(pointDataPtr,lineParameters);
	std::cout << "Parameters size: " << lineParameters.size() << std::endl ;
	std::cout<<"Least squares line parameters: [n_x,n_y,a_x,a_y]\n\t [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
	std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<std::endl;
	dotProd = lineParameters[0]*(-dy) + lineParameters[1]*dx;
	std::cout<<"\tDot product of real and computed line normals[+-1=correct]: "<<dotProd<<std::endl;
	dotProd = (-dy)*(lineParameters[2]) + dx*lineParameters[3];
	std::cout<<"\tCheck if computed point is on real line [0=correct]: "<<dotProd<<std::endl;
	
	//A RANSAC estimate of the line parameters
	clock_t startTime = clock();
	std::vector< std::vector<Point_2>::const_pointer > inliers ;	
	double usedData = Ransac< Point_2, double >::run( lpEstimator, 
													  pointData,													  
													  lineParameters,
													  inliers,
													  desiredProbabilityForNoOutliers ) ;

	clock_t endTime = clock();
	std::cout<<"RANSAC running time: "<<(double)(endTime-startTime)/CLOCKS_PER_SEC<<std::endl;
  	
	if(lineParameters.empty())
	std::cout<<"RANSAC failed to obtain an estimate.\n";
	else
	{
		std::cout<<"RANSAC line parameters [n_x,n_y,a_x,a_y]\n\t [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<std::endl;
		dotProd = lineParameters[0]*(-dy) + lineParameters[1]*dx;
		std::cout<<"\tDot product of real and computed line normals[+-1=correct]: "<<dotProd<<std::endl;
	dotProd = (-dy)*(lineParameters[2]) + dx*lineParameters[3];
		std::cout<<"\tCheck if computed point is on real line [0=correct]: "<<dotProd<<std::endl;
		std::cout<<"\tPercentage of points which were used for final estimate: "<<usedData<<std::endl;
	}

	// LeastKthSquares estimate of the line parameters
	startTime = clock();
	std::vector< double > lksSqResiduals ;		
	bool success = LeastKthSquares< Point_2, double >::run(	lpEstimator, 
															pointData,													  
															lineParameters,
															lksSqResiduals,
															quantile,
															desiredProbabilityForNoOutliers,
															expectedFractionOfOutliers ) ;

	endTime = clock();
	std::cout<<"LeastKthSquares running time: "<<(double)(endTime-startTime)/CLOCKS_PER_SEC<<std::endl;
  	
	if(lineParameters.empty())
	std::cout<<"LeastKthSquares failed to obtain an estimate.\n";
	else
	{
		// Scale estimation
		std::cout << "Original Scale = " << noiseStd << std::endl ;
		double scaleMedianSE = ScaleEstimation::medianSE< double >( lksSqResiduals, 4 ) ;
		std::cout << "Scale MedianSE = " << scaleMedianSE << std::endl ;
		double scaleMadSE = ScaleEstimation::madSE< double >( lksSqResiduals, 4 ) ;
		std::cout << "Scale MadSE    = " << scaleMadSE << std::endl ;
		double scaleMsSE = ScaleEstimation::msSE< double >( lksSqResiduals, 4, quantile ) ;
		std::cout << "Scale MsSE     = " << scaleMsSE << std::endl ;


		std::cout<<"LeastKthSquares line parameters [n_x,n_y,a_x,a_y]\n\t [ "<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]"<<std::endl;
		dotProd = lineParameters[0]*(-dy) + lineParameters[1]*dx;
		std::cout<<"\tDot product of real and computed line normals[+-1=correct]: "<<dotProd<<std::endl;
		dotProd = (-dy)*(lineParameters[2]) + dx*lineParameters[3];
		std::cout<<"\tCheck if computed point is on real line [0=correct]: "<<dotProd<<std::endl;
		std::cout<<"\tPercentage of points which were used for final estimate: "<<usedData<<std::endl;
	}

	pointData.clear();
	
	return EXIT_SUCCESS;
}
