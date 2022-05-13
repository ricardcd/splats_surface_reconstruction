
template< typename DataType, typename ParametersType >
bool BiasedLeastKthSquares< DataType, ParametersType >::run(	const ModelEstimator< DataType, ParametersType > *modelEstimator, 
																const std::vector< DataType > &data, 
																const DataType& mainPoint,	// Point that always will be part of the model (should not be part of data!)
																std::vector< ParametersType > &model,
																std::vector< ParametersType > &sqResiduals,																
																const double quantile, 														
																const double desiredProbabilityForNoOutliers, 
																const double expectedFractionOfOutliers, 														
																const int maxDataTrials ) {	

	if ( data.size() < (int)modelEstimator->minSamples() ) {
		// Not enough data!
		return false ;
	}

	int i ;
	bool someSolution = false ; // Sentinel value allowing detection of solution failure
	bool degenerate = true ;	
	srand( ( unsigned ) time( NULL ) ) ; // Seed for the random number generator
	int numData = (int)data.size() ;
    std::vector< typename std::vector< DataType >::const_pointer > curModelData ; // Vector storing the samples to compute the model at each iteration
	std::vector< ParametersType > curModelParameters ; // Vector storing the parameters of the model at each iteration
	double minKthSqResidual = std::numeric_limits<double>::infinity() ;	

	// Compute minimum number of points an structure is (supposed to be) made of inside x
	int k = static_cast<int>( floor( numData * quantile ) ) ;

	// Compute the number of iterations to perform given desiredProbabilityForNoOutliers and expectedFractionOfOutliers
	int numIter = static_cast<int>( ceil( log( 1 - desiredProbabilityForNoOutliers ) / log( 1 - pow( (1 - expectedFractionOfOutliers), (int)modelEstimator->minSamples() ) ) ) ) ;
	
	// std::cout << "[LeastKthSquares] numIter = " << numIter << std::endl ;

	//std::vector< DataType > mainPointVector ;
	//mainPointVector.push_back( mainPoint ) ;

	// Find the main point inside the data vector
    typename std::vector< DataType >::const_iterator itMainPoint = std::find( data.begin(), data.end(), mainPoint ) ;
	if ( itMainPoint == data.end() ) {
		std::cout << "Data vector must contain the mainPoint!" << std::endl ;
		return false ;
	}
	int indMainPoint = itMainPoint - data.begin() ;

	for ( int q = 0; q < numIter; q++ ) {

		bool degenerate = true ;
		int count = 0 ;
		while ( degenerate && count < maxDataTrials ) {
			// Clear all previous data
			curModelData.clear() ;

			// Add the mainPoint to biase the model
			curModelData.push_back( &mainPoint ) ;

			// Select a random set of samples to build a model			
			for( i = 0; i < (int)modelEstimator->minSamples()-1; i++ ) {
			// for( i = 0; i < (int)modelEstimator->minSamples(); i++ ) {
				int selectedIndex = (int)( ( (float)rand()/(float)RAND_MAX) * numData + 0.5 ) ;				
				curModelData.push_back( &( data[selectedIndex] ) ) ;							
			}			

			// Check if the selected samples are in a degenerate configuration
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
			std::cout << "[LeastKthSquares] Unable to select a non-degenerate data set, failure." << std::endl ;
			return false ;
		}

		/* If we got here, we have some model */

		// std::cout << "[LeastKthSquares] some model" << std::endl ;

		// Compute residuals
		std::vector< double > curSqResiduals = modelEstimator->sqResiduals( curModelParameters, data ) ;
		
		// std::cout << "[LeastKthSquares] got residuals!" << std::endl ;

		// Sort them in increasing order
		// std::sort( curSqResiduals.begin(), curSqResiduals.end() ) ;

		// std::cout << "[LeastKthSquares] residuals sorted" << std::endl ;

		// Keep the model with the smallest residual for the current point
		double sqKthResidual = curSqResiduals[ k ] ;
		// double sqKthResidual = curSqResiduals[ indMainPoint ] ;
		// std::cout << "[LeastKthSquares] curSqResiduals[ k ]" << std::endl ;
		if ( sqKthResidual < minKthSqResidual ) {
			// std::cout << "[LeastKthSquares] Min residual UPDATED !" << std::endl ;
			minKthSqResidual = sqKthResidual ;
			model.resize( curModelParameters.size() ) ;
			std::copy( curModelParameters.begin(), curModelParameters.end(), model.begin() ) ;
		}

	}

	// std::cout << "[LeastKthSquares] End of main loop --> Computing final residuals!" << std::endl ;
	// Fit the selected model to the data points again to extract the residuals
	sqResiduals.clear() ;	
	std::vector< double > finalResiduals = modelEstimator->sqResiduals( model, data ) ;
	sqResiduals.resize( finalResiduals.size() ) ;
	std::copy( finalResiduals.begin(), finalResiduals.end(), sqResiduals.begin() ) ;
	
	// std::cout << "[LeastKthSquares] Final residuals computed! sqResiduals.size() = " << sqResiduals.size() << " / finalResiduals.size() = " << finalResiduals.size() << std::endl ;

	return true ;
}
