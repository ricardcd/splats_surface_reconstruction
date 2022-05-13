#ifndef LEASTKTHSQUARES_H
#define LEASTKTHSQUARES_H

#include <set>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <algorithm>

#include "ModelEstimator.h"

template< typename DataType, typename ParametersType >
class LeastKthSquares {

public:

	static bool run(	const ModelEstimator< DataType, ParametersType > *modelEstimator, // Input
						const std::vector< DataType > &data, // Input
						std::vector< ParametersType > &model, // Output 
						std::vector< ParametersType > &sqResiduals, // Output 
						const double quantile = 0.2, // Input 
						const double desiredProbabilityForNoOutliers = 0.99, // Input 
						const double expectedFractionOfOutliers = 0.7, // Input 	
						const int maxDataTrials = 1000000 ) ; // Input 

} ;

// Templated class, needs implementation in the same file
#include "LeastKthSquares.cpp"

#endif // LEASTKTHSQUARES_H