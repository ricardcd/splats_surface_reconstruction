#ifndef SPLATSCREATOR_H
#define SPLATSCREATOR_H

/* Includes */
// Related to this project
#include "Splat_3/Splat_3.h"
// CGAL
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_3.h>
	// Min sphere of spheres
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
	// Spatial searching
#include <CGAL/Orthogonal_k_neighbor_search.h>
//#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include "KDTreeSearch/My_Search_traits_3.h"
#include "KDTreeSearch/Fuzzy_capsule_3.h"
// Robust Statistics-related
#include "RobustStatistics/Ransac.h"
#include "RobustStatistics/LeastKthSquares.h"
#include "RobustStatistics/BiasedLeastKthSquares.h"
#include "RobustStatistics/ScaleEstimation.h"
	// Models
#include "RobustStatistics/Models/FitLocalSmoothSurface/LocalSmoothSurfaceModelEstimator.h"
// Boost
#include <boost/filesystem.hpp>
// Other
#include <iostream>
#include <limits>
#include <algorithm>

namespace Splats {

template < class K >
class SplatsCreator {

public:
	/* Typedefs */
	typedef typename K::FT												FT ;
	typedef typename K::Point_3											Point_3 ;
	typedef typename K::Point_2											Point_2 ;
	typedef typename K::Segment_3										Segment_3 ;
	typedef typename K::Vector_3										Vector_3 ;
	typedef typename K::Plane_3											Plane_3 ;
	typedef typename K::Sphere_3										Sphere_3 ;
	typedef typename CGAL::Splat_3< K >									Splat_3 ;
		
	typedef typename CGAL::Min_sphere_annulus_d_traits_3<K>				TraitsMinSphere ;
	typedef typename CGAL::Min_sphere_d<TraitsMinSphere>				Min_Sphere ;
	typedef typename K::Circle_2										Circle_2 ;

	typedef typename CGAL::Monge_via_jet_fitting<K, K>					Monge_via_jet_fitting ;
	typedef typename Monge_via_jet_fitting::Monge_form					Monge_form ;

	typedef typename CGAL::My_Search_traits_3< K >						TreeTraits ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< TreeTraits >	Neighbor_search ;
	// typedef typename Neighbor_search::Tree								NNTree ;
	typedef typename CGAL::Kd_tree<TreeTraits>							NNTree ;
	typedef typename CGAL::Fuzzy_sphere< TreeTraits >					Sphere_range_query ;
	typedef typename CGAL::Fuzzy_capsule_3< TreeTraits >				Fuzzy_capsule_3_range_query ;
	typedef boost::shared_ptr<NNTree>									pNNTree ;



	/* Enums */
	enum NeighborsType{ KNN = 0, RadialNN } ; // Available types of neighborhood
	// enum DistanceType{ MeanNN = 0 } ; // Available types of neighborhood
    enum InputType{ Pointset = 0, OrientedPointset, OrientedPointsetWithScores } ; // Available types of input

	/* Constructors */
	// Default constructor
	SplatsCreator() ;
	// Pointsets
	SplatsCreator( std::vector< Point_3 > &inputPoints ) ;
	// Oriented Pointsets
	SplatsCreator( std::vector< Point_3 > &inputPoints, std::vector< Vector_3 > &inputNormals ) ;
	// Oriented Pointsets with scores
	SplatsCreator( std::vector< Point_3 > &inputPoints, std::vector< Vector_3 > &inputNormals, std::vector< FT > &inputScores ) ;
	// Read the (oriented) pointsets directly from file
	SplatsCreator( std::string filePath, unsigned int inputType = Pointset ) ;
	
	// Setter of default parameters
	void setDefaultParameters() ;

	/* Getters/Setters for the parameters */
	std::vector< Point_3 > getPts() { return m_pts ; }

	int getInputType() { return m_inputType ; }
	void setInputType( int inputType ) { m_inputType = inputType ; } 

	int	getNeighType() { return m_neighType ; }
	void setNeighType( int neighType ) { 
        if ( neighType != this->KNN ) {
			std::cerr << "[SplatsCreator] Unknown neighbor type! Not modifying..." << std::endl ;
			return ;
		}
		m_neighType = neighType ; 
	}

	double getNeighRadius() { return m_neighRadius ; }
	void setNeighRadius( double neighRadius ) { m_neighRadius = neighRadius ; }
	
	// int	getDistType() { return m_distType ; }
	// void setDistType( int distType ) {
	// 	if ( distType != MeanNN ) {
	// 		std::cerr << "[SplatsCreator] Unknown distance type! Not modifying..." << std::endl ;
	// 		return ;
	// 	}
	// 	m_distType = distType ;
	// }
	
	unsigned int getK() { return m_k ; } 
	void setK( unsigned int k ) { 
		if ( k == 0 || k > (unsigned int)m_numPts ) {
			std::cerr << "[SplatsCreator] Neighborhood must be greater than zero and smaller than the number of input points! Not modifying..." << std::endl ;
			return ;
		}
		
		m_k = k ; 
	} 
	
	unsigned int getScaleK() { return m_scaleK ; } 
	void setScaleK( unsigned int k ) { 
		if ( k == 0 || k > (unsigned int)m_numPts ) {
			std::cerr << "[SplatsCreator] Neighborhood must be greater than zero and smaller than the number of input points! Not modifying..." << std::endl ;
			return ;
		}
		
		m_scaleK = k ; 
	} 

	unsigned int getLssDegree() { return m_lssDegree ; }
	void setLssDegree( unsigned int lssDegree ) {
		m_lssDegree = lssDegree ;
	}

	unsigned int getMinInliers() { return m_minInliers ; }
	void setMinInliers( unsigned int minInliers ) { m_minInliers = minInliers ; }

	unsigned int getMinFitPts() { return m_minFitPts ; }
	void setMinFitPts( unsigned int minFitPts ) { m_minFitPts = minFitPts ; }

	double getQuantile() { return m_quantile ; } 
	void setQuantile( double quantile ) {
		if ( quantile < 0.0 || quantile > 1.0 ) {
			std::cerr << "[SplatsCreator] quantile must be between 0 and 1! Not modifying..." << std::endl ;
			return ;
		}

		m_quantile = quantile ;
	}

	double getProbNoOutliers() { return m_probNoOutliers ; }
	void setProbNoOutliers( double probNoOutliers ) {
		if ( probNoOutliers < 0.0 || probNoOutliers > 1.0 ) {
			std::cerr << "[SplatsCreator] The desired probability of taking a sample free of outliers must be between 0 and 1! Not modifying..." << std::endl ;
			return ;
		}

		m_probNoOutliers = probNoOutliers ;	
	}

	double getExpectFracOutliers() { return m_expectFracOutliers ; }
	void setExpectFracOutliers( double expectFracOutliers ) {
		if ( expectFracOutliers < 0.0 || expectFracOutliers > 1.0 ) {
			std::cerr << "[SplatsCreator] Expected Fraction of Outliers must be between 0 and 1! Not modifying..." << std::endl ;
			return ;
		}

		m_expectFracOutliers = expectFracOutliers ;
	}

	bool doComputeScale() { return m_doComputeScale ; }
	void setDoComputeScale( bool doComputeScale ) {
		m_doComputeScale = doComputeScale ;
	}

	bool areScalesPrecomputed() { return m_precomputedScales ; }
	
	double getFixedDistanceThreshold() { return m_distThres ; }
	void setFixedDistanceThreshold( double distThres ) { 
		m_distThres = distThres ;
	}

	double getMaxNoise() { return m_maxNoise ; }
	void setMaxNoise( double maxNoise ) { 
		m_maxNoise = maxNoise ;
	}

	bool setRansacQuantile( FT ransacQuantile ) {
		if ( ransacQuantile < 0 || ransacQuantile > 1 ) {			
			return false ;
		}

		m_ransacQuantile = ransacQuantile ;
		return true ;
	}
	FT getRansacQuantile() {
		return m_ransacQuantile ;
	}

	// Get a string, sumarizing the parameters. Useful to get the 
	std::string getParametersString() ;

	// Compute the neighborhood of a given point
	std::vector< Point_3 > getNeighbors( const Point_3 &query ) ;

	// Compute the capsule neighborhood of a given segment + radius
	std::vector< Point_3 > getNeighbors( const Segment_3 &queryCapSeg, const double& queryCapRad ) ;

	// Compute a single splat from a set of points
	bool computeSplat( const Point_3 &curPoint, Splat_3& splat ) ;

	// Compute the scale of a set of points (using LKS + MSSE)
	double computeScale( const std::vector< Point_3 > &neighs, const Point_3 &curPoint ) ;

	// Compute a biased scale value, where the current point is always part of the model
	// WARNING: neighs will not contain curPoint after running this function...
	double computeBiasedScale( std::vector< Point_3 > &neighs, const Point_3 &curPoint ) ;

	// Precompute all scale values for all the points
	void precomputeScales( const bool verbose = false ) ;

	// Save precomputed scales to file, for latter use
	void saveScales( const std::string& filePath ) ;

	// Load precomputed scales from file
	bool loadScales( const std::string& filePath ) ;

	// Get an estimate for the scale in a neighborhood given precomputed per-point scale values
	FT getPrecomputedScale( const std::vector< Point_3 >& neighs ) ;

	// Compute a given Model from a set of points (using RANSAC)
	bool computeModel( const std::vector< Point_3 > &neighs, const Point_3 &curPoint, const double distThres, std::vector< Point_3 > &inliers, std::vector< double > &modelParam, Monge_form &monge ) ;

	// Compute size of a splat given the inliers of the computed model
	double computeSize( const std::vector< Point_3 > &inliers, const Point_3 &originalPoint, const Point_3 &curPoint ) ;
	
	// Compute the bounding sphere
	Sphere_3 computeBoundingSphere() ;

	// Scales filename from parameters
	std::string scalesString() ;

	// Loads the precomputed scales from file
	bool loadPrecomputedScales( const std::string &dir, const std::string &fileBaseName ) ;

	// Saves the precomputed scales to file, for later reuse
	bool savePrecomputedScales( const std::string &dir, const std::string &fileBaseName ) ;

	// Run the desired variant of the algorithm
	void run( std::vector< Splat_3 > &splats, const bool verbose = true ) ;
		

private:
	std::vector< Point_3 > m_pts ;
	std::vector< Vector_3 > m_normals ;
	std::vector< FT > m_scores ;		
	int m_inputType ;
	int	m_neighType ;
	double m_neighRadius ;
	// int	m_distType ;
	double m_distThres ; // RANSAC distance threshold
	unsigned int m_k ;
	unsigned int m_lssDegree ;	
	int m_minInliers ;
	int m_minFitPts ;
	int m_numPts ;
	double m_quantile ;
	double m_probNoOutliers ;
	double m_expectFracOutliers ;
	bool m_doComputeScale ;
	bool m_precomputedScales ;
	FT m_ransacQuantile ;
	pNNTree m_pTree ;
	double m_maxNoise ;
	std::map< Point_3, FT > m_scalesMap ;
	FT m_probNoOutliersRANSAC ;
	FT m_expectFracOutliersRANSAC ;
	int m_ransacNumIter ;
	unsigned int m_scaleK ;

} ;



/******** Implementation of the methods ********/

template< class K >
SplatsCreator<K>::SplatsCreator(  ) {
	m_pts = std::vector< Point_3 >() ;	
    m_inputType = Pointset ;
	m_numPts = m_pts.size() ;
	m_normals = std::vector< Vector_3 >() ;
	m_scores = std::vector< FT >() ;
	m_scalesMap = std::map< Point_3, FT >() ;
			
	// Default options
	setDefaultParameters() ;	
}



template< class K >
SplatsCreator<K>::SplatsCreator( std::vector< typename K::Point_3 > &inputPoints ) {
	m_pts = inputPoints ;
    m_inputType = Pointset ;
	m_numPts = m_pts.size() ;
	m_normals = std::vector< Vector_3 >() ;
	m_scores = std::vector< FT >() ;
	m_scalesMap = std::map< Point_3, FT >() ;
			
	// Default options
	setDefaultParameters() ;	
}



template< class K >
SplatsCreator<K>::SplatsCreator( std::vector< Point_3 > &inputPoints, 
							  std::vector< Vector_3 > &inputNormals ) {

	m_numPts = inputPoints.size() ;
	if ( inputNormals.size() != m_numPts ) {
		std::cerr << "Input points and input normals must have the same size!" << std::endl ;
		std::cerr << "Not initializing..." << std::endl ;
		return ;
	}
	m_pts = inputPoints ;
    m_inputType = OrientedPointset ;
	m_normals = inputNormals ;
	m_scores = std::vector< FT >() ;	
	m_scalesMap = std::map< Point_3, FT >() ;
			
	// Default options
	setDefaultParameters() ;
}



template< class K >
SplatsCreator<K>::SplatsCreator( std::vector< Point_3 > &inputPoints, 
							  std::vector< Vector_3 > &inputNormals, 
							  std::vector< FT > &inputScores ) {

	if ( inputPoints.size() != inputNormals.size() ) {
		std::cerr << "Input points and input normals must have the same size!" << std::endl ;
		std::cerr << "Not initializing..." << std::endl ;
		return ;
	}
    m_pts = inputPoints ;
    m_inputType = OrientedPointsetWithScores ;
	m_normals = inputNormals ;
	m_scores = inputScores ;		
			
	// Default options
	setDefaultParameters() ;
}



template< class K >
void SplatsCreator<K>::setDefaultParameters() {
    m_neighType = KNN ;
	m_neighRadius = std::numeric_limits<double>::infinity() ;
    // m_distType = MeanNN ;
	m_distThres = 0.1 ;
	m_k = 50 ;
	m_lssDegree = 2 ;
	m_minInliers = 15 ;
	m_minFitPts = 10 ;
	m_quantile = 0.3 ;
	m_probNoOutliers = 0.99 ;
	m_expectFracOutliers = 0.4 ;
	m_doComputeScale = false ;
	m_precomputedScales = false ;
	m_maxNoise = 0.5 ;
	m_ransacQuantile = -1.0 ; // NOTE: To tell if RANSAC quartile is being use, test whether getRansacQuantile > 0
	m_probNoOutliersRANSAC = 0.99 ;
	m_expectFracOutliersRANSAC = 0.5 ;
	
	// Initialize the search engine
	// m_tree.rebuild( points.begin(), points.end() ) ; // rebuild?
	m_pTree = pNNTree( new NNTree( m_pts.begin(), m_pts.end() ) ) ; // rebuild?

	// Number of RANSAC iterations
	// TODO: Allow modification of the m_probNoOutliersRANSAC and m_expectFracOutliersRANSAC
	int coefs = ( ( m_lssDegree + 1 ) * ( m_lssDegree + 2 ) ) / 2 ;
	m_ransacNumIter = static_cast<int>( ceil( log( 1 - m_probNoOutliersRANSAC ) / log( 1 - pow( (1 - m_expectFracOutliersRANSAC ), coefs ) ) ) ) ;

}



template< class K >
std::vector< typename K::Point_3 > SplatsCreator<K>::getNeighbors( const typename K::Point_3 &query ) {

		/* Compute k-NN for the current point */
		std::vector<long unsigned int> answer ;
		// vector<double> distances ;
		std::vector< Point_3 > neighborPoints ;
		
        if ( m_neighType == KNN ) {
			// K-NN query
			Neighbor_search search( *m_pTree, query, m_k+1 ) ;
			// Find k+1 nearest neighbors on the structure (k+1 because first NN will be the point itself)
            for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
				neighborPoints.push_back( it->first ) ;
				// distances.push_back( std::sqrt(it->second) ) ;
			}
		}
        else if ( m_neighType == RadialNN ) {
			// Range query
			Sphere_range_query searchSphere( query, m_neighRadius ) ;
			m_pTree->search( std::back_inserter( neighborPoints ), searchSphere ) ;			
			//for ( std::vector< Point_3 >::iterator itNP = neighborPoints.begin(); 
			//	  itNP != neighborPoints.end(); 
			//	  ++itNP ) {
			//	distances.push_back( CGAL::sqrt( CGAL::squared_distance( query, *itNP ) ) ) ;				
			//}
			
			if ( std::find( neighborPoints.begin(), neighborPoints.end(), query ) == neighborPoints.end() ) {
				std::cout << "Current point not found on its neighborhood..." << std::endl ;
			}
		}

		return neighborPoints ;
}



// Range query inside a capsule
template< class K >
std::vector< typename K::Point_3 > SplatsCreator<K>::getNeighbors( const typename K::Segment_3 &queryCapSeg, const double& queryCapRad ) {
	
	std::vector< Point_3 > neighborPoints ;	
	Fuzzy_capsule_3_range_query searchCapsule( queryCapSeg, queryCapRad ) ;
	m_pTree->search( std::back_inserter( neighborPoints ), searchCapsule ) ;

	// Brute force search
	//double sqRadius = queryCapRad * queryCapRad ;
	//std::vector< Point_3 >::iterator it ;
	//for ( it = m_pts.begin(); it != m_pts.end(); ++it ) {
	//	double sqDist = CGAL::squared_distance( *it, queryCapSeg ) ;		
	//	if ( sqDist < sqRadius ) 
	//		neighborPoints.push_back( *it ) ;
	//}

	return neighborPoints ;

}



// Compute the scale of a set of points (using LKS + MSSE)
template< class K >
double SplatsCreator<K>::computeScale( const std::vector< typename K::Point_3 > &neighs, const typename K::Point_3 &curPoint ) {
	
	if ( neighs.empty() || neighs.size() < m_minFitPts ) {
		// No neighbors, default return value...
		return 999999.9 ;
	}

	// Declarate the model estimator
	LocalSmoothSurfaceModelEstimator * lssEstimator = new LocalSmoothSurfaceModelEstimator( m_minFitPts, 
																							m_lssDegree, 
																							// curPoint, 
																							Point_3(),
																							99999.0 ) ; // This last value is a dummy, will not be used during scale computation
	
	// Least Kth Squares estimation of a model
	std::vector< double > lksSqResiduals ;					
	std::vector< double > lssParameters0 ;	
	// std::cout << "Num Neighbors = " << neighs.size() << std::endl ;
	// NOTE(TODO): It seems there is a bug when neighs.size() is equal to the minimum number of points to generate a model...
	bool success = LeastKthSquares< Point_3, double >::run(	lssEstimator, 
															neighs,
															lssParameters0,
															lksSqResiduals,
															m_quantile,
															m_probNoOutliers,
															m_expectFracOutliers,
															10000 ) ;
	
	// Scale computation (using msse)	
	// std::cout << "Num Residuals= " << lksSqResiduals.size() << std::endl ;
	/*for ( int i = 0; i < lksSqResiduals.size(); i++ ) 
		std::cout << "Residuals[" << i << "]= " << lksSqResiduals[i] << std::endl ;
	std::cout << "Model dim = " << lssParameters0.size() << std::endl ;*/

	// Select model dimension
	int modelDimension = ( ( m_lssDegree + 1 ) * ( m_lssDegree + 2 ) ) / 2 ;
	// std::cout << "modelDimension = " << modelDimension << std::endl ;
		
	double scale = ScaleEstimation::msSE< double >( lksSqResiduals, modelDimension, m_quantile ) ;
	
	// --- Debug (Start) ---
	// std::cout << "Scale = " << scale << std::endl ;

	//// Compute the residuals
	//std::vector< Point_3 > testData ;
	//testData.push_back( curPoint ) ;

	//// Check if the point is an inlier for the model
	//std::vector< double > sqDist = lssEstimator->sqResiduals( lssParameters0, testData ) ;
	//double maxDist = 2.5*scale ;
	//double maxSqDist = maxDist*maxDist ;

	// All residuals
	/*for ( int i = 0; i < lksSqResiduals.size(); i++ ) 
		std::cout << "Residuals[" << i << "]= " << lksSqResiduals[i] << std::endl ;*/

	//if ( sqDist[0] > maxSqDist ) {
	//	// Point is not part of the inliers to the computed model...outlier?
	//	std::cout << "Outlier: sqDist = " << sqDist[0] << std::endl ;
	//	std::cout << "         Dist = " << CGAL::sqrt( sqDist[0] ) << std::endl ;
	//	std::cout << "         maxDist = " << maxDist << std::endl ;
	//	std::cout << "         maxSqDist = " << maxSqDist << std::endl ;
	//	std::cout << "         residuals at q = " << lksSqResiduals[(int)floor( lksSqResiduals.size() * m_quantile ) - 1] << std::endl ;
	//	
	//	return 999999.9 ;
	//}
	/*else {
		std::cout << "INLIER: sqDist = " << sqDist[0] << std::endl ;
		std::cout << "         maxDist = " << maxDist << std::endl ;
		std::cout << "         maxSqDist = " << maxSqDist << std::endl ;
	}*/
	// --- Debug  (End)  ---

	return scale ;

}



template< class K >
double SplatsCreator<K>::computeBiasedScale( std::vector< typename K::Point_3 > &neighs, const typename K::Point_3 &curPoint ) {
	
	if ( neighs.empty() ) {
		// No neighbors, default return value...
		return 999999.9 ;
	}

	// Remove current point from data
    typename std::vector< Point_3 >::iterator it = std::find( neighs.begin(), neighs.end(), curPoint ) ;
	if ( it == neighs.end() ) {
		// Current point is not part of the neighbors...
		return 999999.9 ;
	}

	// Remove current point from neighbors
	// neighs.erase( it ) ;

	// Declarate the model estimator
	LocalSmoothSurfaceModelEstimator * lssEstimator = new LocalSmoothSurfaceModelEstimator( m_minFitPts, 
																							m_lssDegree, 
																							curPoint, 
																							99999.0 ) ; // This last value is a dummy, will not be used during scale computation
	
	// Least Kth Squares estimation of a model
	std::vector< double > lksSqResiduals ;					
	std::vector< double > lssParameters0 ;	
	// std::cout << "Num Neighbors = " << neighs.size() << std::endl ;
	// NOTE(TODO): It seems there is a bug when neighs.size() is equal to the minimum number of points to generate a model...
	bool success = BiasedLeastKthSquares< Point_3, double >::run(	lssEstimator, 
																	neighs,
																	curPoint,
																	lssParameters0,
																	lksSqResiduals,
																	m_quantile,
																	m_probNoOutliers,
																	m_expectFracOutliers, 
																	10000 ) ;
	
	// Scale computation (using msse)	
	// std::cout << "Num Residuals= " << lksSqResiduals.size() << std::endl ;
	/*for ( int i = 0; i < lksSqResiduals.size(); i++ ) 
		std::cout << "Residuals[" << i << "]= " << lksSqResiduals[i] << std::endl ;
	std::cout << "Model dim = " << lssParameters0.size() << std::endl ;*/

	// Select model dimension
	int modelDimension = ( ( m_lssDegree + 1 ) * ( m_lssDegree + 2 ) ) / 2 ;
	// std::cout << "modelDimension = " << modelDimension << std::endl ;
		
	double scale = ScaleEstimation::msSE< double >( lksSqResiduals, modelDimension, m_quantile ) ;

	// std::cout << "Scale = " << scale << std::endl ;

	return scale ;

}



template< class K >
void SplatsCreator<K>::precomputeScales( const bool verbose ) {
	
	// Compute the scale for each point		
	int verbosePercent = 0 ;
	int verboseLastPercent = 0 ;	
	int i = 0 ;
	double meanNeigh = 0 ;
    if ( verbose ) std::cout << "(" << std::setfill(' ') << std::setw(3) << verbosePercent << "%, " << std::setfill(' ') << std::setw(10) << meanNeigh << " mean neighs)";
    typename std::vector< Point_3 >::iterator it ;
	for ( it = m_pts.begin(); it != m_pts.end(); ++it, i++ ) {

		// Show progress on screen (if necessary)
		if ( verbose ) {
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( m_numPts ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b(" << std::setfill(' ') << std::setw(3) << verbosePercent << "%, " << std::setfill(' ') << std::setw(10) << meanNeigh << " mean neighs)";
				verboseLastPercent = verbosePercent ;			
			}
		}

		// Get neighbors
		// Radial search
		/*std::vector< Point_3 > neighs ;
		Sphere_range_query searchSphere( *it, searchRadius ) ;
		m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;*/
		
		// K-NN search
		std::vector< Point_3 > neighs ;
		std::vector< FT > sqDistances ;
		Neighbor_search search( *m_pTree, *it, m_scaleK+1 ) ; // Find k+1 nearest neighbors on the structure (k+1 because first NN will be the point itself)
        for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it ) {
			neighs.push_back( it->first ) ;
			sqDistances.push_back( it->second ) ;
		}

		// Compute scale
		// FT scale = computeScale( neighs, Point_3() ) ;
		FT scale = computeScale( neighs, *it ) ;
		
		// Compute Biased scale
		// FT scale = computeBiasedScale( neighs, *it ) ;

		// Store scale
		m_scalesMap[ *it ] = scale ;
		
		// Update meanNeighbors
		meanNeigh = meanNeigh + ( ( (double)neighs.size() - meanNeigh ) / ( i+1 ) ) ;
	}
	if ( verbose ) 
		std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b(" << std::setfill(' ') << std::setw(3) << verbosePercent << "%, " << std::setfill(' ') << std::setw(10) << meanNeigh << " mean neighs)";

	// Do some statistics
	//if ( verbose ) {
	//	// Sort the scales
	//	std::vector< FT > scales ;
	//	for ( std::map< Point_3, FT >::iterator itS = m_scalesMap.begin(); itS != m_scalesMap.end(); ++itS ) {
	//		scales.push_back( itS->second ) ;
	//	}
	//	std::sort( scales.begin(), scales.end() ) ;

	//	// --- Debug (Start) ---
	//	//ofstream os( "__scales.txt", ios_base::out ) ;
	//	//for ( std::vector< FT >::iterator itFt = scales.begin(); itFt != scales.end(); ++itFt ) {
	//	//	os << *itFt << std::endl ;
	//	//}
	//	//os.close() ;
	//	// --- Debug  (End)  ---

	//	std::cout << "- Some statistics:" << std::endl ;
	//	std::cout << "    - Noise scale at quantile (" << globalQuantile << ") = " << scales[ ceil( static_cast<double>( scales.size() ) * globalQuantile ) ] << std::endl ;
	//	std::cout << "    - Expected noise scale at outliers (global expected fraction of outliers = " << globalExpectedFractionOfOutliers << ") = " << scales[ ceil( static_cast<double>( scales.size() ) * ( 1 - globalExpectedFractionOfOutliers ) ) ] << std::endl ;
	//}

	// Indicate that the object has precomputed scales
	m_precomputedScales = true ;

}


template< class K >
void SplatsCreator<K>::saveScales( const std::string& filePath ) {

    std::ofstream os( filePath.c_str(), std::ios_base::out )  ;
	if ( !os.good() ) {
		std::cerr << "Error, it wasn't possible to open the output file!" << std::endl ;
		std::cerr << "Scales were not saved!" << std::endl ;
		return ;
	}

	int i = 0 ;
    for ( typename std::vector< Point_3 >::iterator it = m_pts.begin(); it != m_pts.end(); ++it, i++ ) {
		if ( i < m_pts.size()-1 )
			os << m_scalesMap[ *it ] << std::endl ;
		else
			os << m_scalesMap[ *it ] ; // Avoid newline at end of file...
	}
	os.close() ;

	// --- Debug (Start) ---
	//ofstream osP( "__PointsWithScales.txt", ios_base::out )  ;
	//std::vector< Point_3 > points ;
	//for ( std::map< Point_3, FT >::iterator itS = m_scalesMap.begin(); itS != m_scalesMap.end(); ++itS ) {
	//	osP << (*itS).first << std::endl ;
	//}
	//osP.close() ;
	// --- Debug  (End)  ---

}



template< class K >
bool SplatsCreator<K>::loadScales( const std::string& filePath ) {

    std::ifstream is( filePath.c_str(), std::ios_base::out ) ;
	if ( !is.good() ) {
		std::cerr << "Error, it wasn't possible to open the input file!" << std::endl ;
		std::cerr << "Unable to load scales!" << std::endl ;
		return false ;
	}
	
	std::vector< FT > scales ;
	while( is.good() ) {
		double scale ;
		if ( is >> scale )
			scales.push_back( scale ) ;
	}

	is.close() ;

	if ( scales.size() != m_pts.size() ) {
		std::cerr << "Error, number of scales and points must be the same: Scales.size() = " << scales.size() << ", Points.size() = " << m_pts.size() << std::endl ;
		return false ;
	}

	m_scalesMap.clear() ;
	for( int i = 0; i < scales.size(); i++ ) {
		m_scalesMap[ m_pts[ i ] ] = scales[ i ] ;
	}	

	
	// Indicate that the object has precomputed scales
	m_precomputedScales = true ;

	return true ;

}



template< class K >
typename K::FT SplatsCreator<K>::getPrecomputedScale( const std::vector< Point_3 >& neighs ) {
	
	std::vector< FT > sqScales ;
    for ( typename std::vector< Point_3 >::const_iterator it = neighs.begin(); it != neighs.end(); ++it ) {
		sqScales.push_back( m_scalesMap[ *it ] * m_scalesMap[ *it ] ) ;	// Square scales before computing scale!
	}

	// --- Debug (Start) ---
	//std::vector< FT > scales ;
	//for ( std::vector< Point_3 >::const_iterator it = neighs.begin(); it != neighs.end(); ++it ) {
	//	scales.push_back( m_scalesMap[ *it ] ) ;
	//}

	//std::cout << "scales:" << std::endl ;
	//for ( int i = 0; i < scales.size(); i++ ) {
	//	std::cout << scales[i] << std::endl ;
	//}
	//std::cout << "Scale (MSSE) = " << ScaleEstimation::msSE< FT >( sqScales, 1, m_quantile ) << std::endl ;	
	//std::cout << "Scale (Median) = " << ScaleEstimation::median< FT >( scales ) << std::endl ;

	//std::sort( sqScales.begin(), sqScales.end() ) ;
	//std::cout << "Scale (at q) = " << sqScales[ static_cast<int>( ceil( neighs.size() * m_quantile ) ) ] << std::endl ;	
	// --- Debug  (End)  ---

	return ScaleEstimation::msSE< FT >( sqScales, 1, m_quantile ) ;
	// return sqScales[ static_cast<int>( ceil( neighs.size() * m_quantile ) ) ] ;

}



// Compute a given Model from a set of points (using RANSAC)
template< class K >
bool SplatsCreator<K>::computeModel( const std::vector< typename K::Point_3 > &neighs, 
									 const Point_3 &curPoint,
									 const double distThres, 
									 std::vector< typename K::Point_3 > &inliers, 
									 std::vector< double > &modelParam, 								  
								     typename SplatsCreator::Monge_form &monge ) {


	inliers.clear() ;
	modelParam.clear() ;
									  
	// Declarate the model estimator
	/*LocalSmoothSurfaceModelEstimator * lssEstimator = new LocalSmoothSurfaceModelEstimator( m_minFitPts, 
																							m_lssDegree, 
																							curPoint, 
																							distThres ) ; */

	LocalSmoothSurfaceModelEstimator * lssEstimator = new LocalSmoothSurfaceModelEstimator( m_minFitPts, 
																							m_lssDegree, 
																							Point_3(), 
																							distThres ) ; 
	
	/* Compute RANSAC Local Smooth Surface */
    std::vector< typename std::vector< Point_3 >::const_pointer > inliersPtr ;
	int usedData = Ransac< Point_3, double >::run(	lssEstimator, 
											 		neighs,									  
													modelParam,
													inliersPtr,
													0.99,
													m_ransacNumIter,
													100 ) ;	

	int minInliers = 0 ;
	// Get the minimum number of inliers (fixed or dependant on neighborhood size)
	if ( m_ransacQuantile > 0.0 ) {
		minInliers = static_cast<int>( ceil( neighs.size() * m_ransacQuantile ) ) ;
		if ( minInliers < m_minInliers ) {
			// Force a minimum number of inliers after all (useful for small neighborhoods)
			minInliers = m_minInliers ;
		}
	}
	else { 		
		minInliers = m_minInliers ;
	}
	
	if ( !inliersPtr.empty() && (int)inliersPtr.size() > minInliers ) { 
				
		for ( int j = 0; j < inliersPtr.size(); j++ ) {
			inliers.push_back( *inliersPtr[ j ] ) ;
		}			
			
		// Recompute the monge from the recovered parameters
		Point_3 origin( modelParam[0], modelParam[1], modelParam[2] ) ;
		Vector_3 d1( modelParam[3], modelParam[4], modelParam[5] ) ;
		Vector_3 d2( modelParam[6], modelParam[7], modelParam[8] ) ;
		Vector_3 n( modelParam[9], modelParam[10], modelParam[11] ) ;
		std::vector< FT > coefs ;
		if ( m_lssDegree > 1 ) {
			coefs.push_back( modelParam[12] ) ;
			coefs.push_back( modelParam[13] ) ;	
		}
		if ( m_lssDegree > 2 ) {
			coefs.push_back( modelParam[14] ) ;
			coefs.push_back( modelParam[15] ) ;	
			coefs.push_back( modelParam[16] ) ;
			coefs.push_back( modelParam[17] ) ;		
		}
		if ( m_lssDegree > 3 ) {
			coefs.push_back( modelParam[18] ) ;		
			coefs.push_back( modelParam[19] ) ;		
			coefs.push_back( modelParam[20] ) ;		
			coefs.push_back( modelParam[21] ) ;		
			coefs.push_back( modelParam[22] ) ;	
		}	
		monge = Monge_form( origin, d1, d2, n, coefs ) ;
		
		return true ;
	}				
	else {
		return false ;
	}

}



template< class K >
void SplatsCreator<K>::run( std::vector< typename CGAL::Splat_3<K> > &splats, const bool verbose ) {
	
	if ( m_pts.size() == 0 ) {
		std::cerr << "No input points on the dataset! Stopping..." << std::endl ;
		return ;
	}
	splats.clear() ;

    typename std::vector< Point_3 >::iterator itP = m_pts.begin() ;
    typename std::vector< Vector_3 >::iterator itN = m_normals.begin() ;
	int verbosePercent = 0 ;
	int verboseLastPercent = 0 ;	
	int i = 0 ;
	if ( verbose ) std::cout << "(" << std::setfill(' ') << std::setw(3) << verbosePercent << "%)" ;
	for ( ; itP != m_pts.end(); ++itP, ++itN, i++ ) {
		
		// Show progress on screen (if necessary)
		if ( verbose ) {
			verbosePercent = (int)floor( ( static_cast<double>( i ) / static_cast<double>( m_numPts ) ) * 100.0 ) ;
			if ( verbosePercent > verboseLastPercent ) {
				std::cout << "\b\b\b\b\b\b(" << std::setfill(' ') << std::setw(3) << verbosePercent << "%)" ;
				verboseLastPercent = verbosePercent ;			
			}
		}

		Splat_3 splat ;
		// std::cout << "Computing Splat" << std::endl ;
		bool success = computeSplat( *itP, splat ) ;

		if ( success ) {
			splats.push_back( splat ) ;
		}

	}
	std::cout << "\b\b\b\b\b\b" << "(100%)" << std::endl ;

}



// Compute a single splat from a set of points
template< class K >
bool SplatsCreator<K>::computeSplat( const typename K::Point_3 &curPoint, typename CGAL::Splat_3<K>& splat ) {

	// Retain current point		
	Point_3 originalPoint = curPoint ;

	// Compute neighbors
	std::vector< Point_3 > neighs = getNeighbors( curPoint ) ;

	if ( neighs.size() > m_minInliers ) {
			
		// Compute scale (if required to)
		double distThres = 0.0 ;
		//if ( m_doComputeScaleIter ) {
		//	double scale = computeScale( neighs, curPoint ) ;
		//	distThres = scale * 2.5 ; // Assuming Gaussian noise
		//} 
		//else {
			if ( m_precomputedScales ) {
				distThres = m_scalesMap[ curPoint ] * 2.5 ;
			}
			else {
			    distThres = m_distThres ;
			}
		//}

		// Compute model
		std::vector< Point_3 > inliers ;
		std::vector< double > modelParam ;
		Monge_form monge ;
		bool success = computeModel( neighs, curPoint, distThres, inliers, modelParam, monge ) ;

		if ( success ) {

			// Update current point to be projected on the surface
			Point_3 newPoint = monge.origin() ; // Update current point

			double size = computeSize( inliers, curPoint, newPoint ) ;
			
			// Create the splat
			splat = Splat_3( monge, size*size ) ;

			return true ;

		}
		else {
			return false ;
		}

	}
	else {
		return false ;
	}
}



// Compute size of a splat given the inliers of the computed model
template< class K >
double SplatsCreator<K>::computeSize(	const std::vector< typename K::Point_3 > &inliers, 
										const typename K::Point_3 &originalPoint,
										const typename K::Point_3 &curPoint ) {
	
	// For now, just KNN distance is available...

	// Compute mean distance to the k-NN
	double total = 0.0 ;
	double minDist = std::numeric_limits<double>::infinity() ;
	int numValid = 0 ;
	for ( int j = 0; j < inliers.size(); j++ ) { 
		if ( inliers[j] != originalPoint ) {
							
			double dist = CGAL::sqrt( CGAL::squared_distance( curPoint, inliers[j] ) ) ;
			if ( dist < minDist ) {
				minDist = dist ;
			}
			total += dist ;
			numValid++ ;
		}
	}	   

	return total / numValid ;

}



template< class K >
typename K::Sphere_3 SplatsCreator<K>::computeBoundingSphere() {

	typedef typename CGAL::Min_sphere_of_spheres_d_traits_3< K, FT > Traits ;
	typedef typename CGAL::Min_sphere_of_spheres_d< Traits > Min_sphere ;
	typedef typename Traits::Sphere Traits_sphere ;
	typedef typename Min_sphere::Cartesian_const_iterator  Cartesian_const_iterator ;

	// Points as 0-sized spheres
	std::vector<Traits_sphere> spheres ;
    for ( typename std::vector<Point_3>::iterator it = m_pts.begin(); it != m_pts.end(); it++)
		spheres.push_back( Traits_sphere( *it,0 ) ) ;
				
	// Compute min sphere
	Min_sphere ms(spheres.begin(),spheres.end()) ;

	// retrieve the center (3D)
	Cartesian_const_iterator cit = ms.center_cartesian_begin() ;
	FT cx, cy, cz ;
	cx = *cit ;
	cy = *(++cit) ;
	cz = *(++cit) ;
	Point_3 boundingSphereCenter = Point_3( cx, cy, cz ) ;

	// Build the sphere
	return Sphere_3( boundingSphereCenter, ms.radius()*ms.radius() ) ;
}



template< class K >
std::string SplatsCreator<K>::scalesString() {
	std::ostringstream ossScales ;
	ossScales << "__Scales" ;					
	ossScales << "_sc-k" << m_scaleK 
			  << "_sc-q" << m_quantile 
	     	  << "_sc-p" << m_probNoOutliers 
			  << "_sc-e" << m_expectFracOutliers 
			  << "_sc-d" << m_lssDegree ; // Splat degree
	ossScales << ".txt" ;

	return ossScales.str() ;
}



template < class K >
bool SplatsCreator<K>::savePrecomputedScales( const std::string &dir, const std::string &fileBaseName ) {
	if ( !m_precomputedScales ) {
        std::cout << "[SplatsCreator] A problem occurred during scales saving!" << std::endl ;
		return false ;
	}

	std::ostringstream ossScalesFile ;
	ossScalesFile << dir << "/" << fileBaseName << scalesString() ;

    std::ofstream os( ossScalesFile.str().c_str(), std::ios_base::out )  ;
	if ( !os.good() ) {		
		return false ;
	}

	int i = 0 ;
    for ( typename std::vector< Point_3 >::iterator it = m_pts.begin(); it != m_pts.end(); ++it, i++ ) {
		if ( i < m_pts.size()-1 )
			os << m_scalesMap[ *it ] << std::endl ;
		else
			os << m_scalesMap[ *it ] ; // Avoid newline at end of file...
	}
	os.close() ;

	return true ;
}


template < class K >
bool SplatsCreator<K>::loadPrecomputedScales( const std::string &dir, const std::string &fileBaseName ) {
	std::ostringstream ossScalesFile ;
	ossScalesFile << dir << "/" << fileBaseName << scalesString() ;

	if ( !boost::filesystem::exists( ossScalesFile.str() ) ) {
		return false ;
	}

    std::ifstream is( ossScalesFile.str().c_str() )  ;
	if ( !is.good() ) {
		std::cout << "[SplatsCreator::loadPrecomputedScales] Unable to open input Scales file!" << std::endl ;
		return false ;
	}

	std::vector< FT > scales ;
	while( is.good() ) {
		double scale ;
		if ( is >> scale )
			scales.push_back( scale ) ;
	}

	is.close() ;

	if ( scales.size() != m_pts.size() ) {
		std::cerr << "Error, number of scales and points must be the same: Scales.size() = " << scales.size() << ", Points.size() = " << m_pts.size() << std::endl ;
		return false ;
	}

	m_scalesMap.clear() ;
	for( int i = 0; i < scales.size(); i++ ) {
		m_scalesMap[ m_pts[ i ] ] = scales[ i ] ;
	}	

	
	// Indicate that the object has precomputed scales
	m_precomputedScales = true ;

	return true ;
}


template < class K >
std::string SplatsCreator<K>::getParametersString() {
	
	std::string str ;
	std::ostringstream oss( str ) ;
	
	oss << "it" << m_inputType << "_" ;
	oss << "nt" << m_neighType << "_" ;
    if ( m_neighType == KNN )
		oss << "nk" << m_k << "_" ;
	else
		oss << "nr" << m_neighRadius << "_" ;
	
	// oss << "dt" << m_distType << "_" ;	
	oss << "sd" << m_lssDegree << "_" ;
	oss << "sp" << m_minFitPts << "_" ;	
	
	if ( m_doComputeScale ) {
		oss << "sc-q" << m_quantile << "_" ;
		oss << "sc-p" << m_probNoOutliers << "_" ;
		oss << "sc-e" << m_expectFracOutliers << "_" ; // 'e' stands for epsilon
		oss << "sc-k" << m_scaleK << "_" ; 		
	}
	else {
		oss << "rt" << m_distThres << "_" ; // RANSAC distance threshold
	}
	oss << "ri" << m_minInliers ;

	return oss.str() ;

}



} // End namespace Splats

#endif // SPLATSCREATOR_H
