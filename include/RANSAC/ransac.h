/* TODO: This class could be nicer if it was a generic implementation using functors... */

#ifndef  RANSAC_H
#define  RANSAC_H

#include <iostream>
#include <limits>
#include <algorithm>
#include <limits>
#include <math.h>

using namespace std ;

namespace ransac {

	template< class Kernel >
	bool ransacFitCentroid( const std::vector< typename Kernel::Point_3 >& data, 
							const double& distThres, 
							typename Kernel::Point_3& bestModel,
							std::vector< int >& bestModelInliers,
							const int& maxDataTrials = 100, // Max number of attempts to select a non-degenerate data set
							const int& maxTrials = 100		// Maximum number of trials before giving up
						  ) {
		
		typedef typename Kernel::Point_3 Point_3 ;

		double p = 0.99;    // Desired probability of choosing at least one sample free from outliers

		bool someSolution = false ;      // Sentinel value allowing detection of solution failure.
		int trialcount = 0 ;
		int bestscore =  0 ;
		int N = 1 ;            // Dummy initialisation for number of trials.
		bool degenerate = true ;

		while ( N > trialcount && trialcount < maxTrials ) {
			
			// Select at random s datapoints to form a trial model, M.
			// In selecting these points we have to check that they are not in a degenerate configuration.
			int count = 0 ;
			Point_3 model ;
			while ( degenerate && count < maxDataTrials ) {
				// Generate s random indicies in the range 1..npts
				int indA = CGAL::default_random.get_int( 0, static_cast< int >( data.size() ) ) ;
				int indB = CGAL::default_random.get_int( 0, static_cast< int >( data.size() ) ) ;

				// Test that these points are not a degenerate configuration.
				degenerate = ( indA == indB ) ;
	    
				if ( !degenerate ) {
					// Fit model (centroid) to this random selection of data points.
					model = Point_3( ( data[indA].x() + data[indB].x() ) / 2,
									 ( data[indA].y() + data[indB].y() ) / 2,
 									 ( data[indA].z() + data[indB].z() ) / 2 ) ;
						
				}
				
				//Safeguard against being stuck in this loop forever
				count++ ;
			}

			if ( degenerate ) {
				std::cerr << "[RANSAC] Unable to select a nondegenerate data set!" << std::endl ;
			    trialcount++ ;
				break ;
			}
      
			// Once we are out here we should have some kind of model...        
			// Evaluate distances between points and model returning the indices
			// of elements in x that are inliers.  Additionally, if M is a cell
			// array of possible models 'distfn' will return the model that has
			// the most inliers.  After this call M will be a non-cell object
			// representing only one model.
			std::vector< int > inliers ;
			typename std::vector< Point_3 >::const_iterator itData ;
			int i = 0 ;
			for ( itData = data.begin(); itData != data.end(); ++itData, i++ ) {
				if ( CGAL::squared_distance( model, *itData ) <= distThres ) {
					// Is an inlier
					inliers.push_back( i ) ; 
				}
			}

			// Find the number of inliers to this model.
			int ninliers = static_cast<int>( inliers.size() ) ;

			//cout << "[RANSAC] Number of Points / Number of inliers = " << data.size() << " / " << ninliers << endl ;

			if ( ninliers > bestscore ) {   
				// Largest set of inliers so far...
				bestscore = ninliers ;    // Record data for this model
				bestModelInliers = inliers ;
				bestModel = model ;

				// Update estimate of N, the number of trials to ensure we pick,
				// with probability p, a data set with no outliers.
				double fracinliers =  static_cast<double>( ninliers ) / static_cast<double>( data.size() ) ;
				double pNoOutliers = 1 - pow( fracinliers, 2 ) ;
				pNoOutliers = std::max( numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by -Inf
				pNoOutliers = std::min( 1-numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by 0.
				N = static_cast<int>( log(1-p)/log(pNoOutliers) ) ;

				someSolution = true ;
			}

			trialcount++;  
		}

		return ( someSolution && trialcount < maxTrials ) ;

	}



	template< class Kernel >
	bool ransacFitWeightedCentroid( const std::vector< typename Kernel::Point_3 >& data, 
							const std::vector< double >& weights, 
							const double& distThres, 
							typename Kernel::Point_3& bestModel,
							std::vector< int >& bestModelInliers,
							const int& maxDataTrials = 100, // Max number of attempts to select a non-degenerate data set
							const int& maxTrials = 100		// Maximum number of trials before giving up
						  ) {
		
		typedef typename Kernel::Point_3 Point_3 ;

		double p = 0.99;    // Desired probability of choosing at least one sample free from outliers

		bool someSolution = false ;      // Sentinel value allowing detection of solution failure.
		int trialcount = 0 ;
		int bestscore =  0 ;
		int N = 1 ;            // Dummy initialisation for number of trials.
		bool degenerate = true ;

		while ( N > trialcount && trialcount < maxTrials ) {
			
			// Select at random s datapoints to form a trial model, M.
			// In selecting these points we have to check that they are not in a degenerate configuration.
			int count = 0 ;
			Point_3 model ;
			while ( degenerate && count < maxDataTrials ) {
				// Generate s random indicies in the range 1..npts
				int indA = CGAL::default_random.get_int( 0, static_cast<int>( data.size() ) ) ;
				int indB = CGAL::default_random.get_int( 0, static_cast<int>( data.size() ) ) ;

				// Test that these points are not a degenerate configuration.
				degenerate = ( indA == indB ) ;
	    
				if ( !degenerate ) {
					// Fit model (centroid) to this random selection of data points.					
					double sumWeights = weights[indA] + weights[indB] ;
					model = Point_3( ( weights[indA]*data[indA].x() + weights[indB]*data[indB].x() ) / sumWeights,
									 ( weights[indA]*data[indA].y() + weights[indB]*data[indB].y() ) / sumWeights,
 									 ( weights[indA]*data[indA].z() + weights[indB]*data[indB].z() ) / sumWeights ) ;
						
				}
				
				//Safeguard against being stuck in this loop forever
				count++ ;
			}

			if ( degenerate ) {
				std::cerr << "[RANSAC] Unable to select a nondegenerate data set!" << std::endl ;
			    trialcount++ ;
				break ;
			}
      
			// Once we are out here we should have some kind of model...        
			// Evaluate distances between points and model returning the indices
			// of elements in x that are inliers.  Additionally, if M is a cell
			// array of possible models 'distfn' will return the model that has
			// the most inliers.  After this call M will be a non-cell object
			// representing only one model.
			std::vector< int > inliers ;
			//std::vector< Point_3 >::const_iterator itData ;
			int i = 0 ;
			for ( i = 0; i < data.size(); i++ ) {
				//std::cout << "[RANSAC] weights[" << i << "] = "<< weights[i] << "   distance = " << CGAL::squared_distance( model, data[i] ) << "   Mult = " << weights[i] * CGAL::squared_distance( model, data[i] ) << "   distThres = " << distThres << std::endl ;
				//if ( weights[i] * CGAL::squared_distance( model, data[i] ) <= distThres ) {
				if ( CGAL::squared_distance( model, data[i] ) <= distThres ) {
					// Is an inlier
					inliers.push_back( i ) ; 
				}
			}

			// Find the number of inliers to this model.
			int ninliers = static_cast<int>( inliers.size() ) ;

			//cout << "[RANSAC] Number of Points / Number of inliers = " << data.size() << " / " << ninliers << endl ;

			if ( ninliers > bestscore ) {   
				// Largest set of inliers so far...
				bestscore = ninliers ;    // Record data for this model
				bestModelInliers = inliers ;
				bestModel = model ;

				// Update estimate of N, the number of trials to ensure we pick,
				// with probability p, a data set with no outliers.
				double fracinliers =  static_cast<double>( ninliers ) / static_cast<double>( data.size() ) ;
				double pNoOutliers = 1 - pow( fracinliers, 2 ) ;
				pNoOutliers = std::max( numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by -Inf
				pNoOutliers = std::min( 1-numeric_limits< double >::epsilon(), pNoOutliers ) ;		// Avoid division by 0.
				N = static_cast<int>( log(1-p)/log(pNoOutliers) ) ;

				someSolution = true ;
			}

			trialcount++;  
		}

		/*if ( someSolution )
			std::cout << "[RANSAC] Has some solution" << std::endl ;
		else
			std::cout << "[RANSAC] Has not found any solution" << std::endl ;
		if ( trialcount >= maxTrials )
			std::cout << "[RANSAC] Maximum iterations passed" << std::endl ;*/
		return ( someSolution && trialcount < maxTrials ) ;

	}

} // End of namespace ransac

#endif // RANSAC_H