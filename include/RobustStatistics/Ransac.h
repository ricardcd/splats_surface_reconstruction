#ifndef RANSAC_H
#define RANSAC_H

#include <set>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits>
#include <iostream>

#include "ModelEstimator.h"

template< typename DataType, typename ParametersType >
class Ransac {

public:

	static double run(	const ModelEstimator< DataType, ParametersType > *modelEstimator, 
						const std::vector< DataType > &data, 
						std::vector< ParametersType > &model,
						std::vector< typename std::vector< DataType >::const_pointer > &inliers,
						const double desiredProbabilityForNoOutliers = 0.99, 
						const int maxTrials = 1000000,
						const int maxDataTrials = 1000000 ) ;

	// Convenience function, same as before but also returning the outliers
	static double run(	const ModelEstimator< DataType, ParametersType > *modelEstimator, 
						const std::vector< DataType > &data, 
						std::vector< ParametersType > &model,
						std::vector< typename std::vector< DataType >::const_pointer > &inliers,
						std::vector< typename std::vector< DataType >::const_pointer > &outliers,
						const double desiredProbabilityForNoOutliers = 0.99, 
						const int maxTrials = 1000000,
						const int maxDataTrials = 1000000 ) ;


} ;

// Templated class, needs implementation in the same file
#include "Ransac.cpp"

#endif // RANSAC_H