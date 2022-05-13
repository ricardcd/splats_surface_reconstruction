#include <algorithm>
#include <limits>


template< typename DataType, typename ParametersType >
double Ransac< DataType, ParametersType >::run( const ModelEstimator< DataType, ParametersType > *modelEstimator, 
												const std::vector< DataType > &data, 																								
												std::vector< ParametersType > &model,
												std::vector< typename std::vector< DataType >::const_pointer > &inliers,
												const double desiredProbabilityForNoOutliers, 
												const int maxTrials,
												const int maxDataTrials ) {	

	if ( data.size() < (int)modelEstimator->minSamples() ) 
		return -1 ;

	int i ;
	bool someSolution = false ; // Sentinel value allowing detection of solution failure
	int trialcount = 0 ;
	int bestModelNumInliers =  -1 ;
	double N = 1 ;					// Dummy initialisation for number of trials
	bool degenerate = true ;	
	srand( ( unsigned ) time( NULL ) ) ; // Seed for the random number generator
	int numData = (int)data.size() ;
	// std::vector< DataType * > curModelData ; // Vector storing the samples to compute the model at each iteration
	std::vector< typename std::vector<DataType>::const_pointer > curModelData ; // Vector storing the samples to compute the model at each iteration
	std::vector< ParametersType > curModelParameters ; // Vector storing the parameters of the model at each iteration
	bool *inliersFlags = new bool[ numData ] ;
	bool *bestModelInliersFlags = new bool[ numData ] ;

	while ( N > trialcount && trialcount < maxTrials ) {

		bool degenerate = true ;
		int count = 0 ;
		while ( degenerate && count < maxDataTrials ) {
			// Clear all previous data
			curModelData.clear() ;

			// Select a random set of samples to build a model
			//std::cout << "[Ransac] Index selection" << std::endl ;
			for( i = 0; i < (int)modelEstimator->minSamples(); i++ ) {
				int selectedIndex = (int)( ( (float)rand()/(float)RAND_MAX) * numData + 0.5 ) ;
				//std::cout << " " << selectedIndex << std::endl ;//<< ": " << data[selectedIndex].p << data[selectedIndex].n << std::endl ;	
				curModelData.push_back( &( data[selectedIndex] ) ) ;							
			}
			//std::cout << std::endl ;

			// Check if the selected samples are in a degenerate configuration
			//std::cout << "[RANSAC] First Point = " << curModelData[0]->x() << curModelData[0]->y() << curModelData[0]->z() << std::endl ;
			degenerate = modelEstimator->degenerate( curModelData ) ;
			
			if ( !degenerate ) {
				// Fit the model to the selected data points
				modelEstimator->fit( curModelData, curModelParameters ) ;
			}			
			
			// Safeguard against not finding any non-degenerate model for this data in a specifyed number of loops
			count++ ;
		}
		
		if ( degenerate ) {
			// Unable to select a non-degenerate data set, failure
			// std::cout << "[Ransac] Unable to select a non-degenerate data set, failure." << std::endl ;
			return -1.0 ;
		}

		// If we got here, we have some model
		// Compute the set of input data samples that are coherent with the computed model
		std::fill( inliersFlags, inliersFlags+numData, false ) ;
		int numInliers = 0 ;
		typename std::vector< DataType >::const_iterator itData ;		
		i = 0 ;
		for ( itData = data.begin(); itData != data.end(); ++itData, i++ ) {
			//inliersFlags[ i ] = modelEstimator->distance( curModelParameters, *itData ) < distanceThreshold ;				
			inliersFlags[ i ] = modelEstimator->compatible( curModelParameters, *itData ) ;				
			if ( inliersFlags[ i ] ) {
				// Is an inlier
				numInliers++ ;		
			}
		}

		// std::cout << "numInliers = " << numInliers << std::endl ;

		// Update the best model found so far, if needed
		if ( numInliers > bestModelNumInliers ) {   
			// Largest set of inliers so far...
			bestModelNumInliers = numInliers ;    // Record data for this model
			std::copy( inliersFlags, inliersFlags+numData, bestModelInliersFlags ) ;
			// Save best model?

			// Update estimate of N, the number of trials to ensure we pick,
			// with probability p, a data set with no outliers.
			double fracinliers =  static_cast<double>( numInliers ) / static_cast<double>( data.size() ) ;
			double pNoOutliers = 1 - pow( fracinliers, (int)modelEstimator->minSamples() ) ;
			pNoOutliers = std::max<double>( std::numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by -Inf
			pNoOutliers = std::min<double>( 1-std::numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by 0.
			N = log(1-desiredProbabilityForNoOutliers)/log(pNoOutliers) ;

			// std::cout << "N = " << N << ", N(no int) = " << log(1-desiredProbabilityForNoOutliers)/log(pNoOutliers) << ", fracinliers = " << fracinliers << ", numInliers = " << numInliers << std::endl ;

			someSolution = true ;
		}

		trialcount++ ;
	}

	if ( someSolution ) { // && trialcount < maxTrials ) {
		// std::cout << "[Ransac] Some solution found, computing the least-squares model with the inliers." << std::endl ;

		// Build the inliers vector
		for ( i = 0; i < numData; i++ ) {
			if ( bestModelInliersFlags[i] ) {
				inliers.push_back( &( data[i] ) ) ;
			}
		}

		// Least-squares estimate with all the inliers
		modelEstimator->leastSquaresFit( inliers, model ) ;

		//std::cout << "[Ransac] Some solution found, computing the least-squares model with the inliers (DONE)." << std::endl ;

		// Return the percentage of inliers for the best model
		return static_cast< double >( bestModelNumInliers ) / static_cast< double >( data.size() ) ;
	}
	else {
		// std::cout << "[Ransac] No solution found, or maximum number of iterations exceeded, failure." << std::endl ;
		return -1.0 ;
	}

}



template< typename DataType, typename ParametersType >
double Ransac< DataType, ParametersType >::run( const ModelEstimator< DataType, ParametersType > *modelEstimator, 
												const std::vector< DataType > &data, 																								
												std::vector< ParametersType > &model,
												std::vector< typename std::vector< DataType >::const_pointer > &inliers,
												std::vector< typename std::vector< DataType >::const_pointer > &outliers,
												const double desiredProbabilityForNoOutliers, 
												const int maxTrials,
												const int maxDataTrials ) {	

	if ( data.size() < (int)modelEstimator->minSamples() ) 
		return -1 ;

	int i ;
	bool someSolution = false ; // Sentinel value allowing detection of solution failure
	int trialcount = 0 ;
	int bestModelNumInliers =  -1 ;
	double N = 1 ;					// Dummy initialisation for number of trials
	bool degenerate = true ;	
	srand( ( unsigned ) time( NULL ) ) ; // Seed for the random number generator
	int numData = (int)data.size() ;
	// std::vector< DataType * > curModelData ; // Vector storing the samples to compute the model at each iteration
	std::vector< typename std::vector<DataType>::const_pointer > curModelData ; // Vector storing the samples to compute the model at each iteration
	std::vector< ParametersType > curModelParameters ; // Vector storing the parameters of the model at each iteration
	bool *inliersFlags = new bool[ numData ] ;
	bool *bestModelInliersFlags = new bool[ numData ] ;

	while ( N > trialcount && trialcount < maxTrials ) {

		bool degenerate = true ;
		int count = 0 ;
		while ( degenerate && count < maxDataTrials ) {
			// Clear all previous data
			curModelData.clear() ;

			// Select a random set of samples to build a model
			//std::cout << "[Ransac] Index selection" << std::endl ;
			for( i = 0; i < (int)modelEstimator->minSamples(); i++ ) {
				int selectedIndex = (int)( ( (float)rand()/(float)RAND_MAX) * numData + 0.5 ) ;
				//std::cout << " " << selectedIndex << std::endl ;//<< ": " << data[selectedIndex].p << data[selectedIndex].n << std::endl ;	
				curModelData.push_back( &( data[selectedIndex] ) ) ;							
			}
			//std::cout << std::endl ;

			// Check if the selected samples are in a degenerate configuration
			//std::cout << "[RANSAC] First Point = " << curModelData[0]->x() << curModelData[0]->y() << curModelData[0]->z() << std::endl ;
			degenerate = modelEstimator->degenerate( curModelData ) ;
			
			if ( !degenerate ) {
				// Fit the model to the selected data points
				modelEstimator->fit( curModelData, curModelParameters ) ;
			}			
			
			// Safeguard against not finding any non-degenerate model for this data in a specifyed number of loops
			count++ ;
		}
		
		if ( degenerate ) {
			// Unable to select a non-degenerate data set, failure
			std::cout << "[Ransac] Unable to select a non-degenerate data set, failure." << std::endl ;
			return -1.0 ;
		}

		// If we got here, we have some model
		// Compute the set of input data samples that are coherent with the computed model
		std::fill( inliersFlags, inliersFlags+numData, false ) ;
		int numInliers = 0 ;
		typename std::vector< DataType >::const_iterator itData ;		
		i = 0 ;
		for ( itData = data.begin(); itData != data.end(); ++itData, i++ ) {
			//inliersFlags[ i ] = modelEstimator->distance( curModelParameters, *itData ) < distanceThreshold ;				
			inliersFlags[ i ] = modelEstimator->compatible( curModelParameters, *itData ) ;				
			if ( inliersFlags[ i ] ) {
				// Is an inlier
				numInliers++ ;		
			}
		}

		// std::cout << "numInliers = " << numInliers << std::endl ;

		// Update the best model found so far, if needed
		if ( numInliers > bestModelNumInliers ) {   
			// Largest set of inliers so far...
			bestModelNumInliers = numInliers ;    // Record data for this model
			std::copy( inliersFlags, inliersFlags+numData, bestModelInliersFlags ) ;
			// Save best model?

			// Update estimate of N, the number of trials to ensure we pick,
			// with probability p, a data set with no outliers.
			double fracinliers =  static_cast<double>( numInliers ) / static_cast<double>( data.size() ) ;
			double pNoOutliers = 1 - pow( fracinliers, (int)modelEstimator->minSamples() ) ;
			pNoOutliers = std::max<double>( std::numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by -Inf
			pNoOutliers = std::min<double>( 1-std::numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by 0.
			N = log(1-desiredProbabilityForNoOutliers)/log(pNoOutliers) ;

			// std::cout << "N = " << N << ", N(no int) = " << log(1-desiredProbabilityForNoOutliers)/log(pNoOutliers) << ", fracinliers = " << fracinliers << ", numInliers = " << numInliers << std::endl ;

			someSolution = true ;
		}

		trialcount++ ;
	}

	if ( someSolution ) { // && trialcount < maxTrials ) {
		// std::cout << "[Ransac] Some solution found, computing the least-squares model with the inliers." << std::endl ;

		// Build the inliers vector
		for ( i = 0; i < numData; i++ ) {
			if ( bestModelInliersFlags[i] ) {
				inliers.push_back( &( data[i] ) ) ;
			}
		}

		// Build the outliers vector
		for ( i = 0; i < numData; i++ ) {
			if ( !bestModelInliersFlags[i] ) {
				outliers.push_back( &( data[i] ) ) ;
			}
		}

		// Least-squares estimate with all the inliers
		modelEstimator->leastSquaresFit( inliers, model ) ;

		//std::cout << "[Ransac] Some solution found, computing the least-squares model with the inliers (DONE)." << std::endl ;

		// Return the percentage of inliers for the best model
		return static_cast< double >( bestModelNumInliers ) / static_cast< double >( data.size() ) ;
	}
	else {
		if ( !someSolution )
			std::cout << "[Ransac] No solution found, failure." << std::endl ;
		else 
			std::cout << "[Ransac] Maximum number of iterations exceeded, failure." << std::endl ;
		return -1.0 ;
	}

}