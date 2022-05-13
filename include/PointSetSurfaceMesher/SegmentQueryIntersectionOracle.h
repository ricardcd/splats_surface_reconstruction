#ifndef SEGMENTQUERYINTERSECTIONORACLE_H
#define SEGMENTQUERYINTERSECTIONORACLE_H


/* Includes */
// CGAL
	// Spatial searching
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include "SplatsCreation/KDTreeSearch/My_Search_traits_3.h"
#include "SplatsCreation/KDTreeSearch/Fuzzy_capsule_3.h"
#include "GeneralDistanceModels/Euclidean_distance_segment_3_point.h"
	// Min sphere of spheres
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
	// Multiple return value
#include <CGAL/Object.h>
// Robust Statistics-related
#include "RobustStatistics/Ransac.h"
#include "RobustStatistics/LeastKthSquares.h"
#include "RobustStatistics/ScaleEstimation.h"
// Models
#include "LocalSurfaceFit/LBQModelEstimator.h"
#include "LocalSurfaceFit/WLBQModelEstimator.h"
#include "LocalSurfaceFit/WLBQSegmentModelEstimator.h"
#include "RobustStatistics/Models/FitPlane/PlaneModelEstimator.h"
#include "RobustStatistics/ModelEstimator.h"
// Boost
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
// Other
#include <iostream>
#include <limits>
#include <algorithm>

//#include "SingleStepViewer/SingleStepViewer.h"




template < class K >
class SegmentQueryIntersectionOracle {

public :

	typedef typename K::FT												FT ;
	typedef typename K::Point_3											Point_3 ;
	typedef typename K::Point_2											Point_2 ;
	typedef typename K::Ray_3											Ray_3 ;
	typedef typename K::Segment_3										Segment_3 ;
	typedef typename K::Vector_3										Vector_3 ;
	typedef typename K::Plane_3											Plane_3 ;
	typedef typename K::Sphere_3										Sphere_3 ;
	typedef CGAL::Object												Object ;
	
	typedef typename CGAL::Min_sphere_of_spheres_d_traits_3< K, FT >	Traits ;
	typedef typename CGAL::Min_sphere_of_spheres_d< Traits >			Min_sphere ;
	typedef typename Traits::Sphere										Traits_sphere ;
	typedef typename Min_sphere::Cartesian_const_iterator				Cartesian_const_iterator ;

	typedef typename CGAL::My_Search_traits_3< K >						TreeTraits ;
	typedef typename CGAL::Orthogonal_k_neighbor_search< TreeTraits >	Neighbor_search ;
	typedef typename CGAL::Kd_tree<TreeTraits>							NNTree ;
	typedef typename CGAL::Fuzzy_sphere< TreeTraits >					Sphere_range_query ;
	typedef typename CGAL::Fuzzy_capsule_3< TreeTraits >				Fuzzy_capsule_3_range_query ;
	typedef boost::shared_ptr<NNTree>									pNNTree ;
	typedef CGAL::K_neighbor_search<TreeTraits, 
								Euclidean_distance_segment_3_point<K> >	Segment_nearest_neighbor_search ;

	typedef Fitting::LocalBivariateQuadric<K>							LocalBivariateQuadric ;



	// Enums (Available types of scale of noise computation)
	enum RansacThreshold{ Fixed = 0,					// Fixed ransac threshold (user-input) 
						  IterativeScaleComputation,	// Scale is computed at each query segment (VERY computationally expensive)
						  PrecomputeScales,				// Scales are precomputed before intersection queries are computed
						  UsePrecomputedScales } ;		// Use already precomputed scales



	// Constructors
	SegmentQueryIntersectionOracle() ;
	SegmentQueryIntersectionOracle( const std::vector< Point_3 > &inputPoints ) ;
	SegmentQueryIntersectionOracle( const std::vector< Point_3 > &inputPoints, 
									const int &degree,
									const bool &forceLBQOnSegment,
									const FT& distThres,
									const int& minInliers,
									const FT& quantile,
									const FT& probNoOutliers,
									const FT& expectFracOutliers,									
									const int& ransacThresholdType,
									const int& scaleK,
									const FT& capsuleRadiusPercent,
									const FT& weightsGaussianH,
									const FT& scaleRadiusFactor,
									const int &scaleEstimator,
									const bool &computeDensity,
									const double &ransacT,
									const double &intConsFact, 
									const int &intConsNum,
									const bool &debugOutput,
									const bool &debugVisualize ) ;
	
	// Sets the default parameters
	void setDefaultParameters() ;

	// Set debug flags
	void setDebugOutput( const bool &flag ) { m_debugOutput = flag ; }
//	void setDebugVisualize( const bool &flag ) { m_debugVisualize = flag ; }

	// Get all points
	std::vector< Point_3 > getPts() { return m_pts ; }

	// Get the neighbors inside the capsule defined by the query segment and the user-defined radius
	std::vector< Point_3 > getNeighborsOnCapsule( const Segment_3 &queryCapSeg, const double &capsuleRadius ) ;

	// Check the validity of an intersection point by checking that at least "numAgreeing" points from "neighs" are within the radial neighborhood described by "intersection" and "distThres"
	bool intersectionHasConsensus( const std::vector< Point_3 > &neighs, const Point_3 &intersection, const double &distThres, const int &numAgreeing = 1 ) ;

	// Compute the scale of a set of points (using LKS + MSSE)
	FT computeScale( const std::vector< Point_3 > &neighs, const Point_3 &curPoint ) ;

	// Precompute all scale values for all the points
	void precomputeScalesAndDensities( const FT& searchRadius, const bool verbose = false ) ;

	// Save precomputed scales to file, for latter use
	void saveScales( const std::string& filePath ) ;

	// Save precomputed densities to file, for latter use
	void saveDensities( const std::string& filePath ) ;

	// Load precomputed scales from file
	bool loadScales( const std::string& filePath ) ;

	// Load precomputed densities from file
	bool loadDensities( const std::string& filePath ) ;
	
	// Compute the model
	bool computeModel(  const std::vector< Point_3 > &neighs, 
						const Segment_3 &segment, 
						const double &distThres, 
						std::vector< Point_3 > &inliers, 
						std::vector< Point_3 > &outliers, 
						std::vector< double > &modelParam, 
						LocalBivariateQuadric &lbq,
						const double &searchRadius = 0.0 ) ;
	
	bool computePlaneModel( const std::vector< Point_3 > &neighs, 
							const Segment_3 &segment, 
							const double &distThres, 
							std::vector< Point_3 > &inliers, 
							std::vector< Point_3 > &outliers, 
							std::vector< double > &modelParam, 
							Plane_3 &plane,
							const double &searchRadius = 0.0 ) ;

	// Get an estimate for the scale in a neighborhood given precomputed per-point scale values
	FT getPrecomputedScale( const std::vector< Point_3 > &neighs ) ;

	// Compute the bounding sphere
	Sphere_3 computeBoundingSphere() ;

	// Compute the query capsule's radius from a percent of the bounding sphere radius
	FT computeCapsuleRadiusFromBBPercent( const FT &bBPercent ) ;

	// Parameters string (useful for creating filenames for different executions)
	std::string getParametersString() ;
	
	// Useful debug info to show on screen before starting processing...
	void debugShowParameters() ;

	// Debug: save invalid intersection information
	void debugInvalidIntersect( const int &m_iter, const Point_3 &intersectionPoint, const Segment_3 &segment, const LocalBivariateQuadric &lbq, const std::vector< Point_3 > &neighs, const std::vector< Point_3 > &inliers ) const ;

	/* Main function */
	// Function to execute previous to the meshing
	void prepare( const std::string &outputDir = std::string(), const std::string &precomputedScalesFilePath = std::string() ) ;

	// Computes the intersection point for the query segment, if any
	Object intersection( const Segment_3 &querySegment ) ;

	// Get a value for the radius from the precomputed local density (Using a set of precomputed points) 
	FT getAdaptiveRadius( const std::vector< Point_3 > &pts ) ;

	// Get a value for the radius from the precomputed local density (Using nearest neighbor to the query segment)
	FT getAdaptiveRadiusNN( const Segment_3 &querySegment ) ;
	
	// Computes the intersection point for the query segment, if any
	Object intersection( const Ray_3 &queryRay ) ;

	// Convenience function, computes the factorial of a number
	// From Knuth's "The Art of Computer Programming, 3rd Edition, Volume 2: Seminumerical Algorithms"
	// Modifications from StackOverflow.com: http://stackoverflow.com/questions/1838368/calculating-the-amount-of-combinations
	static unsigned long long gcd(unsigned long long x, unsigned long long y) {
		while (y != 0)
		{
			unsigned long long t = x % y;
			x = y;
			y = t;
		}
		return x;
	}

	static unsigned long long choose(unsigned long long n, unsigned long long k) {
		if (k > n)
			throw std::invalid_argument("invalid argument in choose");
		unsigned long long r = 1;
		for (unsigned long long d = 1; d <= k; ++d, --n)
		{
			unsigned long long g = gcd(r, d);
			r /= g;
			unsigned long long t = n / (d / g);
			if (r > std::numeric_limits<unsigned long long>::max() / t)
			   throw std::overflow_error("overflow in choose") ;
			r *= t;
		}
		return r;
	}

private:

	// The intersection function is divided in the following 3:
	std::vector< Point_3 > intersection_getNeighbors( const Segment_3& querySegment,
													  FT &newRadius ) ;
	FT intersection_getRansacThreshold( const Segment_3& querySegment, 
										const std::vector< Point_3 > &pts ) ;
	Object intersection_computeIntersection( const Segment_3 &querySegment,
											 const std::vector< Point_3 > &pts, 
											 const double &distThres,
											 const double &radius, 
											 int &trialCount ) ;
	
	/* Attributes */
	std::vector< Point_3 > m_pts ;
	int m_degree ;
	int m_numPts ;
	FT m_distThres ; // RANSAC distance threshold
	int m_minInliers ;
	FT m_quantile ;
	FT m_probNoOutliers ;
	FT m_expectFracOutliers ;
	bool m_hasPrecomputedScales ;
	pNNTree m_pTree ;
	std::map< Point_3, FT > m_scalesMap ;
	FT m_capsuleRadius ;
	FT m_capsuleRadiusPercent ;
	FT m_scaleRadiusFactor ;
	int m_ransacThresholdType ;
	bool m_debugOutput ;
//	bool m_debugVisualize ;
	unsigned int m_iter ; // Debug iterator
//	SingleStepViewer<K> *m_stepViewer ;
	FT m_weightsGaussianHFactor ;
	int m_scaleEstimator ;
	int m_scaleK ;
	bool m_computeDensity ;
	std::map< Point_3, FT > m_densitiesMap ;	
	double m_ransacT ;
	double m_newRadius ;
	double m_intersectionConsensusFactor ;
	int m_intersectionConsensusNumber ;
	int m_ransacNumIter ;
	bool m_forceLBQOnSegment ;
	
} ; // End of class SegmentQueryIntersectionOracle<K>



/******** Implementation of the methods ********/
template< class K >
SegmentQueryIntersectionOracle<K>::SegmentQueryIntersectionOracle() {
	m_pts = std::vector< Point_3 >() ;	
	// m_inputType = InputType::Pointset ;
	m_numPts = m_pts.size() ;
	//m_normals = std::vector< Vector_3 >() ;
	// m_scores = std::vector< FT >() ;
	m_scalesMap = std::map< Point_3, FT >() ;
			
	// Default options
	setDefaultParameters() ;	
	m_iter = 0 ; // Debug iterator

	boost::filesystem::path dir( "./_OUT" ) ;
	if (boost::filesystem::create_directory( dir ) )
		std::cout << "Created Debug Output Directory" << "\n";

}



template< class K >
SegmentQueryIntersectionOracle<K>::SegmentQueryIntersectionOracle( const std::vector< typename K::Point_3 > &inputPoints ) {
	m_pts = inputPoints ;
	//m_inputType = InputType::Pointset ;
	m_numPts = m_pts.size() ;
	//m_normals = std::vector< Vector_3 >() ;
	// m_scores = std::vector< FT >() ;
	m_scalesMap = std::map< Point_3, FT >() ;
			
	// Default options
	setDefaultParameters() ;	
	m_iter = 0 ; // Debug iterator

	boost::filesystem::path dir( "./_OUT" ) ;
	if (boost::filesystem::create_directory( dir ) )
		std::cout << "Created Debug Output Directory" << "\n";

}



template < class K >
SegmentQueryIntersectionOracle<K>::SegmentQueryIntersectionOracle( const std::vector< Point_3 > &inputPoints, 
																   const int &degree,
																   const bool &forceLBQOnSegment,
																   const FT& distThres,
																   const int& minInliers,
																   const FT& quantile,
																   const FT& probNoOutliers,
																   const FT& expectFracOutliers,
																   const int& ransacThresholdType,
																   const int& scaleK,
																   const FT& capsuleRadiusPercent,
																   const FT& weightsGaussianHFactor,
																   const FT& scaleRadiusFactor,
																   const int &scaleEstimator,
																   const bool &computeDensity,
																   const double &ransacT,
																   const double &intConsFact, 
																   const int &intConsNum,
																   const bool &debugOutput,
																   const bool &debugVisualize ) :
m_pts( inputPoints ),
m_degree( degree ),
m_forceLBQOnSegment( forceLBQOnSegment ),
m_distThres( distThres ),
m_minInliers( minInliers ),
m_quantile( quantile ),
m_probNoOutliers( probNoOutliers ),
m_expectFracOutliers( expectFracOutliers ),
m_scaleK( scaleK ),
m_capsuleRadiusPercent( capsuleRadiusPercent ),
m_weightsGaussianHFactor( weightsGaussianHFactor ),
m_scaleRadiusFactor( scaleRadiusFactor ), 
m_hasPrecomputedScales( false ),
m_scaleEstimator( scaleEstimator ),
m_computeDensity( computeDensity ),
m_ransacT( ransacT ),
m_intersectionConsensusFactor( intConsFact ),
m_intersectionConsensusNumber( intConsNum ),
m_iter(0),  // Debug iterator
m_debugOutput( debugOutput )
//m_debugVisualize( debugVisualize )
{
	m_capsuleRadius = computeCapsuleRadiusFromBBPercent( capsuleRadiusPercent ) ;
	m_ransacThresholdType = ransacThresholdType ;
	if ( m_ransacThresholdType < 0 || m_ransacThresholdType > 4 ) {
		m_ransacThresholdType = 0 ; // Default: Fixed threshold
	}
	m_numPts = m_pts.size() ;

	// Initialize the search engine
	m_pTree = pNNTree( new NNTree( m_pts.begin(), m_pts.end() ) ) ;

	if( m_debugOutput ) {
		boost::filesystem::path dir( "./_OUT" ) ;
		if (boost::filesystem::create_directory( dir ) )
			std::cout << "Created Debug Output Directory" << "\n";
	}

//	if( m_debugVisualize ) {
//		Sphere_3 bs = computeBoundingSphere() ;
//		m_stepViewer = new SingleStepViewer<K>( bs.center(), CGAL::sqrt( bs.squared_radius() ), m_capsuleRadius ) ;
//	}

	if ( m_degree == 2 ) {
		m_ransacNumIter = static_cast<int>( ceil( log( 1 - m_probNoOutliers ) / log( 1 - pow( (1 - m_expectFracOutliers ), 6 ) ) ) ) ;
	}
	else {
		m_ransacNumIter = static_cast<int>( ceil( log( 1 - m_probNoOutliers ) / log( 1 - pow( (1 - m_expectFracOutliers ), 3 ) ) ) ) ;
	}
}



template< class K >
void SegmentQueryIntersectionOracle<K>::setDefaultParameters() {

	m_distThres = 0.1 ;
	m_minInliers = 15 ;
	m_quantile = 0.2 ;
	m_probNoOutliers = 0.99 ;
	m_expectFracOutliers = 0.4 ;
	// m_doComputeScaleIter = true ;
	m_hasPrecomputedScales = false ;
	m_capsuleRadiusPercent = 0.1 ;
	m_capsuleRadius = computeCapsuleRadiusFromBBPercent( m_capsuleRadiusPercent ) ;
	// m_weightsGaussianH = m_capsuleRadius ; // weightsGaussianHFactor = 1 
	m_scaleRadiusFactor = 1.0 ;
	m_numPts = m_pts.size() ;
	m_scaleEstimator = 1 ; // Default scale estimator is the LBQ
	m_computeDensity = false ;
	m_scaleK = 200 ;
	m_intersectionConsensusFactor = 0.5 ;
	m_intersectionConsensusNumber = 3 ;

	// Initialize the search engine
	m_pTree = pNNTree( new NNTree( m_pts.begin(), m_pts.end() ) ) ;

//	if( m_debugVisualize )
//		m_stepViewer = new SingleStepViewer<K>() ;

}



// Range query inside a capsule
template< class K >
std::vector< typename K::Point_3 > SegmentQueryIntersectionOracle<K>::getNeighborsOnCapsule( const typename K::Segment_3 &queryCapSeg, const double& queryCapRad ) {
	
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

	
	// TODO: There is some kind of bug when dealing with large query segments (in the order of squared_length() > 1e18 )
	//       The number of neighbors is not the same when using KD-Tree range search and the brute force approach.
	//		 Need to check this problem, but apparently it does not happen in any other case...


	// Test results against brute force search
	//double sqRadius = queryCapRad * queryCapRad ;
	//std::vector< Point_3 > neighborPoints2 ;	
	//std::vector< Point_3 >::iterator it ;
	//for ( it = m_pts.begin(); it != m_pts.end(); ++it ) {
	//	double sqDist = CGAL::squared_distance( *it, queryCapSeg ) ;		
	//	if ( sqDist < sqRadius ) 
	//		neighborPoints2.push_back( *it ) ;
	//}

	//std::cout << "KdTree = " << neighborPoints.size() << " / BF = " << neighborPoints2.size() << std::endl ;
	//if ( neighborPoints.size() != neighborPoints2.size() ) {
	//	std::cout << "Number of neighbors is not the same!!" <<std::endl ;		
	//	cin.get() ;
	//}

	return neighborPoints ;

}



// Compute the scale of a set of points (using LKS + MSSE)
template< class K >
typename K::FT SegmentQueryIntersectionOracle<K>::computeScale( const std::vector< typename K::Point_3 > &neighs, const typename K::Point_3 &curPoint ) {
	
	if ( neighs.empty() || neighs.size() < m_minInliers ) {
		// No neighbors, default return value...
		return 999999.9 ;
	}

	ModelEstimator< Point_3, double > *mEstimator ;
	int modelDimension = 0 ;
	if ( m_degree == 1 ) {
		// Plane
		modelDimension = 3 ;		
		mEstimator = new PlaneModelEstimator( 999999999.9 ) ;		// This value (distance threshold) is a dummy, will not be used during scale computation
	}
	else {
		// Local Bivariate Quadric				
		if ( m_scaleEstimator == 0 ) {
			// Declarate the model estimator
			mEstimator = new PlaneModelEstimator( 999999999.9 ) ;		// This value (distance threshold) is a dummy, will not be used during scale computation
			// Select model dimension
			modelDimension = 3 ;		
		}
		else if ( m_scaleEstimator == 1 ) {
			// Declarate the model estimator
			mEstimator = new LBQModelEstimator( 999999999.9	) ;			// This value (distance threshold) is a dummy, will not be used during scale computation
	 										 		
			// Select model dimension
			modelDimension = 6 ;		
		}
		else {
			// Declarate the model estimator
			mEstimator = new WLBQModelEstimator( curPoint,
												 999999999.9,			// This value (distance threshold) is a dummy, will not be used during scale computation
	 											 m_weightsGaussianHFactor * m_capsuleRadius ) ; 		
			// Select model dimension
			modelDimension = 6 ;		
		}
	}

	// --- Debug (Start) ---
	// std::cout << "Num Neighbors = " << neighs.size() << std::endl ;
	// NOTE(TODO): It seems there is a bug when neighs.size() is equal to the minimum number of points to generate a model...
	// --- Debug  (End)  ---

	// Least Kth Squares estimation of a model
	std::vector< double > lksSqResiduals ;					
	std::vector< double > lssParameters0 ;	
	bool success = LeastKthSquares< Point_3, double >::run(	mEstimator,
															neighs,
															lssParameters0,
															lksSqResiduals,
															m_quantile,
															m_probNoOutliers,
															m_expectFracOutliers,
															10000 ) ;

	// Scale computation (using msse)
	/*std::cout << "Num Residuals= " << lksSqResiduals.size() << std::endl ;
	for ( int i = 0; i < lksSqResiduals.size(); i++ ) 
		std::cout << "Residuals[" << i << "]= " << lksSqResiduals[i] << std::endl ;
	std::cout << "Model dim = " << lssParameters0.size() << std::endl ;*/

	double scale = ScaleEstimation::msSE< double >( lksSqResiduals, modelDimension, m_quantile ) ;
	/*double scale2 = ScaleEstimation::medianSE< double >( lksSqResiduals, modelDimension ) ;
	double scale3 = ScaleEstimation::median< double >( lksSqResiduals ) ;
		
	std::cout << "Scale MSSE = " << scale << std::endl ;
	std::cout << "Scale MedianSE = " << scale2 << std::endl ;
	std::cout << "Scale Median = " << scale3 << std::endl ;*/
	

	return scale ;

}



template< class K >
void SegmentQueryIntersectionOracle<K>::precomputeScalesAndDensities( const FT& searchRadius, const bool verbose ) {
	
	// Compute the scale for each point	
	
	int verbosePercent = 0 ;
	int verboseLastPercent = 0 ;	
	int i = 0 ;
	double meanNeigh = 0 ;
	FT minDensity = std::numeric_limits<double>::infinity(), maxDensity = 0 ;
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
		std::vector< Point_3 > neighs ;
		std::vector< FT > sqDistances ;
		if ( m_ransacThresholdType == 2 ) {
			// Radial neighbor search
			Sphere_range_query searchSphere( *it, searchRadius ) ;
			m_pTree->search( std::back_inserter( neighs ), searchSphere ) ;
		}
		else if ( m_ransacThresholdType == 3 || m_ransacThresholdType == 0 ) {
			// Nearest neighbor search
			Neighbor_search search( *m_pTree, *it, m_scaleK+1 ) ; // Find k+1 nearest neighbors on the structure (k+1 because first NN will be the point itself)
			for( typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
				neighs.push_back( it->first ) ;
				sqDistances.push_back( it->second ) ;
			}
		}

		// Compute scale
		if ( m_ransacThresholdType > 0 ) {
			FT scale = computeScale( neighs, *it ) ;
		
			// Store scale
			m_scalesMap[ *it ] = scale ;
		}
		
		// Store density, update min/max values
		if ( m_computeDensity ) {
			
			FT density = 0 ;
			if ( m_ransacThresholdType == 2 ) {			
				// Compute density (as neighbors' size)
				density = neighs.size() ;
			}
			else {
				// Compute density (as mean distance to K-NN)
				FT totalDist = 0.0 ;
				for( typename std::vector< FT >::iterator itD = sqDistances.begin()+1; itD != sqDistances.end(); ++itD) {
					totalDist += CGAL::sqrt( *itD ) ;
				}
				density = totalDist / m_scaleK ;
			}

			m_densitiesMap[ *it ] = density ;
			if ( density < minDensity )
				minDensity = density ;
			if ( density > maxDensity )
				maxDensity = density ;
		}

		// Update meanNeighbors
		meanNeigh = meanNeigh + ( ( (double)neighs.size() - meanNeigh ) / ( i+1 ) ) ;
	}
	if ( verbose ) 
		std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b(" << std::setfill(' ') << std::setw(3) << "100%, " << std::setfill(' ') << std::setw(10) << meanNeigh << " mean neighs)" << std::endl ;
	
	// Scale densities between 0-1
	if ( m_computeDensity ) {
		typename std::map< Point_3, FT >::iterator itD ;
		for ( itD = m_densitiesMap.begin(); itD != m_densitiesMap.end(); ++itD ) {
			itD->second = ( itD->second - minDensity ) / ( maxDensity - minDensity ) ;
		}
	}

	// Mark that the object has precomputed scales
	if ( m_ransacThresholdType > 1 ) {
		m_hasPrecomputedScales = true ;
	}

}



template< class K >
void SegmentQueryIntersectionOracle<K>::saveScales( const std::string& filePath ) {

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
	//std::ofstream osP( "__PointsWithScales.txt", ios_base::out )  ;
	//std::vector< Point_3 > points ;
	//for ( std::map< Point_3, FT >::iterator itS = m_scalesMap.begin(); itS != m_scalesMap.end(); ++itS ) {
	//	osP << (*itS).first << std::endl ;
	//}
	//osP.close() ;
	// --- Debug  (End)  ---

}



template< class K >
void SegmentQueryIntersectionOracle<K>::saveDensities( const std::string& filePath ) {

	std::ofstream os( filePath.c_str(), std::ios_base::out )  ;
	if ( !os.good() ) {
		std::cerr << "Error, it wasn't possible to open the output file!" << std::endl ;
		std::cerr << "Densities were not saved!" << std::endl ;
		return ;
	}

	int i = 0 ;
	for ( typename std::vector< Point_3 >::iterator it = m_pts.begin(); it != m_pts.end(); ++it, i++ ) {
		if ( i < m_pts.size()-1 )
			os << m_densitiesMap[ *it ] << std::endl ;
		else
			os << m_densitiesMap[ *it ] ; // Avoid newline at end of file...
	}
	os.close() ;

}



template< class K >
bool SegmentQueryIntersectionOracle<K>::loadScales( const std::string& filePath ) {

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
	m_hasPrecomputedScales = true ;

	return true ;

}


template< class K >
bool SegmentQueryIntersectionOracle<K>::loadDensities( const std::string& filePath ) {

	std::ifstream is( filePath.c_str(), std::ios_base::out ) ;
	if ( !is.good() ) {
		std::cerr << "Error, it wasn't possible to open the input file!" << std::endl ;
		std::cerr << "Unable to load densities!" << std::endl ;
		return false ;
	}
	
	std::vector< FT > densities ;
	while( is.good() ) {
		double density ;
		if ( is >> density )
			densities.push_back( density ) ;
	}

	is.close() ;

	if ( densities.size() != m_pts.size() ) {
		std::cerr << "Error, number of scales and points must be the same: Scales.size() = " << densities.size() << ", Points.size() = " << m_pts.size() << std::endl ;
		return false ;
	}

	m_densitiesMap.clear() ;
	for( int i = 0; i < densities.size(); i++ ) {
		m_densitiesMap[ m_pts[ i ] ] = densities[ i ] ;
	}	
		
	return true ;

}



// Compute a given Model from a set of points (using RANSAC)
template< class K >
bool SegmentQueryIntersectionOracle<K>::computeModel( const std::vector< typename K::Point_3 > &neighs, 
													  const Segment_3 &segment,
													  const double &distThres, 
													  std::vector< typename K::Point_3 > &inliers, 
													  std::vector< typename K::Point_3 > &outliers, 
													  std::vector< double > &modelParam, 								  
													  typename SegmentQueryIntersectionOracle<K>::LocalBivariateQuadric &lbq,
													  const double &searchRadius ) {


	inliers.clear() ;
	outliers.clear() ;
	modelParam.clear() ;
									  
	// Declarate the model estimator
	double curSearchRadius = searchRadius ;
	if ( searchRadius == 0.0 ) {
		curSearchRadius = m_capsuleRadius ;
	}

	ModelEstimator< Point_3, double > *mEstimator ;
	
	if ( m_forceLBQOnSegment ) {
		if ( m_debugOutput )
			std::cout << "    LBQ forced to follow the segment" << std::endl ;
		mEstimator = new WLBQSegmentModelEstimator( segment, distThres, m_weightsGaussianHFactor * curSearchRadius ) ; 
	}
	else {
		if ( m_debugOutput )
			std::cout << "    Free LBQ" << std::endl ;
		mEstimator = new LBQModelEstimator( distThres ) ; 
	}

	// Using combinatorics, the number of iterations is guided by the combinations without repetition of numNeighs in groups of 6
	//int ransacMaxIter = 5000000 ; //  NumIter( P=0.99, epsilon=0.9 ) = log( 1-0.99 ) / log( 1-(1-0.9)^6 ) = 4605168 iter
	//int numNeighs = neighs.size() ;
	//unsigned long long ransacNumIter = 0 ;
	//try {
	//	ransacNumIter = choose( numNeighs, 6 ) ;
	//}
	//catch ( const std::overflow_error& e ) {
	//	ransacNumIter = ransacMaxIter ;
	//}

	//// However, if this number is too big, we clip it
	//std::cout << "    RANSAC num. iterations (unfiltered)= " << ransacNumIter << std::endl ;
	//if ( ransacNumIter > std::numeric_limits<int>::max() || ransacNumIter > ransacMaxIter ) 
	//	ransacNumIter = ransacMaxIter ;

	//std::cout << "    RANSAC num. iterations = " << ransacNumIter << std::endl ;

	// int ransacNumIter = 5000000 ; //  NumIter( P=0.99, epsilon=0.9 ) = log( 1-0.99 ) / log( 1-(1-0.9)^6 ) = 4605168 iter
	// int ransacNumIter = 7000 ; //  NumIter( P=0.99, epsilon=0.7 ) = log( 1-0.99 ) / log( 1-(1-0.7)^6 ) = 7000 iter

	int ransacNumIter = m_ransacNumIter ;

	/* Compute RANSAC Local Smooth Surface */
	std::vector< typename std::vector< Point_3 >::const_pointer > inliersPtr ;	
	std::vector< typename std::vector< Point_3 >::const_pointer > outliersPtr ;	
	int usedData = Ransac< Point_3, double >::run(	mEstimator, 
											 		neighs,									  
													modelParam,
													inliersPtr,
													outliersPtr,
													0.99,
													(int)ransacNumIter,
													100 ) ;	

	int minInliers = 0 ;
	// Get the minimum number of inliers (fixed or dependant on neighborhood size)
	//if ( m_quantile > 0.0 ) {
	//	minInliers = static_cast<int>( ceil( neighs.size() * m_quantile ) ) ;
	//	if ( minInliers < m_minInliers ) {
	//		// Force a minimum number of inliers after all (useful for small neighborhoods)
	//		minInliers = m_minInliers ;
	//	}
	//}
	//else { 		
		minInliers = m_minInliers ;
	// }

	if ( m_debugOutput ) {
		std::cout << "    Inliers/MinInliers = " << inliersPtr.size() << " / " << minInliers << std::endl ;
	}
	
	if ( !inliersPtr.empty() && (int)inliersPtr.size() > minInliers ) { 
				
		for ( int j = 0; j < inliersPtr.size(); j++ ) {
			inliers.push_back( *inliersPtr[ j ] ) ;
		}

		for ( int j = 0; j < outliersPtr.size(); j++ ) {
			outliers.push_back( *outliersPtr[ j ] ) ;
		}

		// Recompute the lbq from the recovered parameters			
		lbq = WLBQSegmentModelEstimator::parameters2LBQ( modelParam ) ;
				
		return true ;
	}				
	else {
		return false ;
	}

}



template< class K >
bool SegmentQueryIntersectionOracle<K>::computePlaneModel(	const std::vector< typename K::Point_3 > &neighs, 
															const Segment_3 &segment,
															const double &distThres, 
															std::vector< typename K::Point_3 > &inliers, 
															std::vector< typename K::Point_3 > &outliers, 
															std::vector< double > &modelParam, 								  
															typename K::Plane_3 &plane,
															const double &searchRadius ) {


	inliers.clear() ;
	outliers.clear() ;
	modelParam.clear() ;
									  
	// Declarate the model estimator
	double curSearchRadius = searchRadius ;
	if ( searchRadius == 0.0 ) {
		curSearchRadius = m_capsuleRadius ;
	}
		
	PlaneModelEstimator *mEstimator = new PlaneModelEstimator( distThres ) ; 
	
	int ransacNumIter = m_ransacNumIter ;

	/* Compute RANSAC Local Smooth Surface */
	std::vector< typename std::vector< Point_3 >::const_pointer > inliersPtr ;	
	std::vector< typename std::vector< Point_3 >::const_pointer > outliersPtr ;	
	int usedData = Ransac< Point_3, double >::run(	mEstimator, 
											 		neighs,									  
													modelParam,
													inliersPtr,
													outliersPtr,
													0.99,
													(int)ransacNumIter,
													100 ) ;	

	int minInliers = 0 ;
	// Get the minimum number of inliers (fixed or dependant on neighborhood size)
	//if ( m_quantile > 0.0 ) {
	//	minInliers = static_cast<int>( ceil( neighs.size() * m_quantile ) ) ;
	//	if ( minInliers < m_minInliers ) {
	//		// Force a minimum number of inliers after all (useful for small neighborhoods)
	//		minInliers = m_minInliers ;
	//	}
	//}
	//else { 		
		minInliers = m_minInliers ;
	// }

	if ( m_debugOutput ) {
		std::cout << "    Inliers/MinInliers = " << inliersPtr.size() << " / " << minInliers << std::endl ;
	}
	
	if ( !inliersPtr.empty() && (int)inliersPtr.size() > minInliers ) { 
				
		for ( int j = 0; j < inliersPtr.size(); j++ ) {
			inliers.push_back( *inliersPtr[ j ] ) ;
		}

		for ( int j = 0; j < outliersPtr.size(); j++ ) {
			outliers.push_back( *outliersPtr[ j ] ) ;
		}

		// Recompute the lbq from the recovered parameters			
		plane = Plane_3( modelParam[0], modelParam[1], modelParam[2], modelParam[3] ) ;
				
		return true ;
	}				
	else {
		return false ;
	}

}



template< class K >
typename K::FT SegmentQueryIntersectionOracle<K>::getPrecomputedScale( const std::vector< Point_3 >& neighs ) {
	
	std::vector< FT > sqScales ;
	for ( typename std::vector< Point_3 >::const_iterator it = neighs.begin(); it != neighs.end(); ++it ) {
		sqScales.push_back( m_scalesMap[ *it ] * m_scalesMap[ *it ] ) ;	// Square scales before computing scale!		
	}

	// --- Debug (Start) ---
	//std::sort( sqScales.begin(), sqScales.end() ) ;
	//std::cout << "Scales " << std::endl ;
	//for ( std::vector< FT >::iterator itS = sqScales.begin(); itS != sqScales.end(); ++itS ) {
	//	std::cout << *itS << std::endl ;
	//}
	//std::cout << std::endl ;
	// --- Debug  (End)  ---

	return ScaleEstimation::msSE< FT >( sqScales, 1, m_quantile ) ;

}



template< class K >
typename K::Sphere_3 SegmentQueryIntersectionOracle<K>::computeBoundingSphere() {

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
typename K::FT SegmentQueryIntersectionOracle<K>::computeCapsuleRadiusFromBBPercent( const FT& bBPercent ) {

	return CGAL::sqrt( computeBoundingSphere().squared_radius() ) * bBPercent ;
}



template < class K >
std::string SegmentQueryIntersectionOracle<K>::getParametersString() {
	
	std::string str ;
	std::ostringstream oss( str ) ;
	
	oss << "cr" << m_capsuleRadiusPercent << "_" ;	
	
	if ( m_ransacThresholdType > 0 ) {
		if ( m_ransacThresholdType == 2 ) 
			oss << "sc-rad" << "_" ;
		else if ( m_ransacThresholdType == 3 ) 
			oss << "sc-k" << m_scaleK << "_" ;

		oss << "sc-q" << m_quantile << "_" ;
		oss << "sc-p" << m_probNoOutliers << "_" ;
		oss << "sc-e" << m_expectFracOutliers << "_" ; // 'e' stands for epsilon
		oss << "sc-f" << m_scaleRadiusFactor << "_" ; // 'f' stants for factor
	}
	else {
		oss << "rt" << m_distThres << "_" ; // RANSAC distance threshold
	}
	if ( m_computeDensity )
		oss << "d_" ;
	oss << "ri" << m_minInliers << "_" ;
	oss << "if" << m_intersectionConsensusFactor << "_" ;
	oss << "in" << m_intersectionConsensusNumber ;

	return oss.str() ;

}



// Useful debug info to show on screen before starting processing...
template < class K >
void SegmentQueryIntersectionOracle<K>::debugShowParameters() {

	std::cout << "- Segment Query Intersection Oracle Parameters:" << std::endl ;
	std::cout << "    - Num. Points = " << m_pts.size() << std::endl ;
	if ( m_ransacThresholdType == 0 ) {
		std::cout << "    - Fixed RANSAC threshold = " << m_distThres << std::endl ;
	}
	else if ( m_ransacThresholdType == 1 ) {
		std::cout << "    - Iterative Scale Computation" << std::endl ;
	}
	else if ( m_ransacThresholdType == 2 || m_ransacThresholdType == 3 ) {
		std::cout << "    - Precompute Scales before meshing" << std::endl ;
	}
	else if ( m_ransacThresholdType == 4 ) {
		std::cout << "    - Scales will be loaded from a previous computation" << std::endl ;
	}

	if ( m_ransacThresholdType > 1 ) {
		if ( m_scaleEstimator == 0 )
			std::cout << "    - Scales estimator is a plane" << std::endl ;		
		else if ( m_scaleEstimator == 1 )
			std::cout << "    - Scales estimator is a Local Bivariate Quadric" << std::endl ;
		else
			std::cout << "    - Scales estimator is a Weighted Local Bivariate Quadric" << std::endl ;							
	}

	std::cout << "    - Minimum inliers = " << m_minInliers << std::endl ;
	if ( m_ransacThresholdType > 0 )
		std::cout << "    - RansacT = " << m_ransacT << std::endl ;
	std::cout << "    - Quantile = " << m_quantile << std::endl ;
	std::cout << "    - LKS probability of no outliers = " << m_probNoOutliers << std::endl ;
	std::cout << "    - LKS expected fraction of outliers = " << m_expectFracOutliers << std::endl ;
	std::cout << "    - Number of lks/ransac iterations = " << m_ransacNumIter << std::endl ;
	std::cout << "    - Query capsule radius = " << m_capsuleRadius << std::endl ;
	std::cout << "    - Gaussian weighting H Factor = " << m_weightsGaussianHFactor << std::endl ;	
	if ( m_ransacThresholdType > 0 && m_ransacThresholdType < 3 || m_computeDensity ) 		
		std::cout << "    - Enlargement factor during scale/density computation = " << m_scaleRadiusFactor << std::endl ; 

}



template < class K >
void SegmentQueryIntersectionOracle<K>::prepare( const std::string &dir, const std::string& fileName ) {

	/* Build Scales filename from parameters */
	bool computedDensities = false ;
	std::ostringstream ossScales ;
	if ( m_ransacThresholdType > 1 ) {		
		ossScales << dir << "/"
					<< fileName 
					<< "__Scales_" ;					

		if ( m_ransacThresholdType == 2 ) 
			ossScales << "sc-rad" << m_capsuleRadiusPercent << "_f" << m_scaleRadiusFactor ;			
		else if ( m_ransacThresholdType == 3 ) 
			ossScales << "sc-k" << m_scaleK << "_" ;

		ossScales << "sc-q" << m_quantile << "_" 
				  << "sc-p" << m_probNoOutliers << "_" 
				  << "sc-e" << m_expectFracOutliers ;

		if( m_scaleEstimator == 0 || m_degree == 1 ) {
			ossScales << "_plane" ;
		}
		else if ( m_scaleEstimator == 1 ) {
			ossScales << "_lbq" ;
		}
		else {
			ossScales << "_wlbq" ;
		}
		ossScales << ".txt" ;

		// Try to load them, if already computed
		std::cout << ossScales.str() << std::endl ;
		if ( boost::filesystem::exists( ossScales.str() ) ) {
			std::cout << "- Loading precomputed scales..." << std::endl ;
			// Load precomputed scales
			if ( !this->loadScales( ossScales.str() ) ) {
				std::cerr << "Error, cannot open input scales file! Precomputing them..." << std::endl ;				
			}
			m_hasPrecomputedScales = true ;
		}

	}

	/* Build Densities filename from parameters */
	std::ostringstream ossDensities ;
	if ( m_computeDensity ) {		
		ossDensities  << dir << "/"
					  << fileName 
					  << "__Densities_" ;					  
		if ( m_ransacThresholdType == 2 ) 
			ossDensities << "rad" << m_capsuleRadiusPercent << "_f" << m_scaleRadiusFactor ;
		else if ( m_ransacThresholdType == 3 || m_ransacThresholdType == 0 ) 
			ossDensities << "k" << m_scaleK ;	
		
		ossDensities << ".txt" ;

		// Try to load them, if already computed
		if ( boost::filesystem::exists( ossDensities.str() ) ) {
			std::cout << "- Loading precomputed densities..." << std::endl ;
			// Load precomputed scales
			if ( !this->loadDensities( ossDensities.str() ) ) {
				std::cerr << "Error, cannot open input scales file! Precomputing them..." << std::endl ;
				
			}
			computedDensities = true ;
		}
	}

	// Compute global scale?
	if ( !m_hasPrecomputedScales && ( m_ransacThresholdType == 2 || m_ransacThresholdType == 3 ) ) {
		if ( m_computeDensity ) 
			std::cout << "- Computing scale (std of noise) and density for each point " << std::flush ;
		else
			std::cout << "- Computing scale (std of noise) for each point " << std::flush ;

		// Search radius is scaleCapRadFactor times the capsules' radius used during meshing
		FT searchRadius = m_capsuleRadius * m_scaleRadiusFactor ;
		this->precomputeScalesAndDensities( searchRadius, true ) ;

		// Save the scales computed for later reuse
		std::cout << "- Saving scales at: " << ossScales.str() << std::endl ;

		this->saveScales( ossScales.str() ) ;

		if ( m_computeDensity ) {
			// Densities have been computed also, save them
			this->saveDensities( ossDensities.str() ) ;
			computedDensities = true ;
		}

	}

	if ( !computedDensities && m_computeDensity ) {
		// Scales haven't been computed, so densities haven't been computed either
		// Compute densities
		FT searchRadius = m_capsuleRadius * m_scaleRadiusFactor ;
		std::cout << "- Computing a density measure for each point " << std::flush ;
		this->precomputeScalesAndDensities( searchRadius, true ) ;
	}

	// Compute global density? (should have been computed during the scale computation)
	if ( !computedDensities && m_computeDensity ) {
		// Save densities to file		
		this->saveDensities( ossDensities.str() ) ;
	}

}



template < class K >
void SegmentQueryIntersectionOracle<K>::debugInvalidIntersect(  const int& m_iter,
																const Point_3& intersectionPoint,
																const Segment_3& segment,
																const typename SegmentQueryIntersectionOracle<K>::LocalBivariateQuadric &lbq,																
																const std::vector< Point_3 >& neighs, 
																const std::vector< Point_3 >& inliers ) const {
			
	std::cout << "        Invalid..." << std::endl ;

	std::ostringstream invalidIntersectFileName ;						
	invalidIntersectFileName << "./_OUT/__invalidIntersect_" << m_iter << ".xyz" ;
	std::ofstream ofsi( invalidIntersectFileName.str(), std::ios_base::out ) ;
	ofsi << intersectionPoint << std::endl ;
	ofsi.close() ;

	std::ostringstream invalidIntersectInliersFileName ;						
	invalidIntersectInliersFileName << "./_OUT/__invalidIntersectInliers_" << m_iter << ".xyz" ;
	std::ofstream ofsii( invalidIntersectInliersFileName.str(), std::ios_base::out ) ;
	for ( int i = 0; i < inliers.size(); i++ ) {
		ofsii << inliers[i] << std::endl ;
	}
	ofsii.close() ;

	std::ostringstream invalidIntersectNeighFileName ;						
	invalidIntersectNeighFileName << "./_OUT/__invalidIntersectNeigh_" << m_iter << ".xyz" ;
	std::ofstream ofsin( invalidIntersectNeighFileName.str(), std::ios_base::out ) ;
	for ( int i = 0; i < neighs.size(); i++ ) {
		ofsin << neighs[i] << std::endl ;
	}
	ofsin.close() ;

	std::ostringstream invalidIntersectSegmentFileName ;						
	invalidIntersectSegmentFileName << "./_OUT/__invalidIntersectSegment_" << m_iter << ".xyz" ;
	std::ofstream ofsis( invalidIntersectSegmentFileName.str(), std::ios_base::out ) ;
	ofsis << segment << std::endl ;
	ofsis.close() ;

	std::ostringstream invalidIntersectSplatFileName ;						
	invalidIntersectSplatFileName << "./_OUT/__invalidIntersectSplat_" << m_iter << ".splat" ;
	std::ofstream ofsisp( invalidIntersectSplatFileName.str(), std::ios_base::out ) ;
	lbq.print( ofsisp ) ;
	ofsisp.close() ;
						
}
	
	

template < class K >
bool SegmentQueryIntersectionOracle<K>::intersectionHasConsensus( const std::vector< Point_3 > &neighs, 
																  const Point_3 &intersection, 
																  const double &distThres, 
																  const int &numAgreeing ) 
{

	int numPtsInRadialNeigh = 0 ;
	FT sqDist = distThres*distThres ;
	typename std::vector< Point_3 >::const_iterator it = neighs.begin() ;
	bool finish = false ;
	while ( it != neighs.end() && numPtsInRadialNeigh < numAgreeing ) {

		if ( CGAL::squared_distance( intersection, *it ) < sqDist ) {
			numPtsInRadialNeigh++ ;
		}
		
		++it ;
	}
	
	if ( m_debugOutput ) 
		std::cout << "    Number of int. consensus pts =  " << numPtsInRadialNeigh << std::endl ;

	return numPtsInRadialNeigh >= numAgreeing ;
}



template < class K >
typename SegmentQueryIntersectionOracle<K>::FT SegmentQueryIntersectionOracle<K>::getAdaptiveRadius( const std::vector< Point_3 > &pts ) {

	// Get minimum scale density
	typename std::vector< Point_3 >::const_iterator it ;
	FT minDensity = std::numeric_limits<double>::infinity() ;
	for ( it = pts.begin(); it != pts.end(); ++it ) {
		FT curDensity = m_densitiesMap[ *it ] ;
		
		if ( curDensity < minDensity ) 
			minDensity = curDensity ;
	}

	if ( m_debugOutput )
		std::cout << "    Minimum density = " << minDensity << std::endl ;

	// Enlarge the radius by a factor depending on the minimum Density on the neighborhood
	return m_capsuleRadius + ( ( ( m_capsuleRadius*m_scaleRadiusFactor ) - m_capsuleRadius ) * ( 1 - minDensity ) ) ;

}



template < class K >
typename SegmentQueryIntersectionOracle<K>::FT SegmentQueryIntersectionOracle<K>::getAdaptiveRadiusNN( const typename K::Segment_3 &querySegment ) {

	// Get nearest neighbor
	Segment_nearest_neighbor_search search( *m_pTree, querySegment, 1 ) ;
	typename Segment_nearest_neighbor_search::iterator it = search.begin() ;

	// Get the density of the point
	FT density = m_densitiesMap[ it->first ] ;
	if ( m_debugOutput )
		std::cout << "    Segment's NN density = " << density << std::endl ;

	// Enlarge the radius by a factor depending on the minimum Density on the neighborhood
	return m_capsuleRadius + ( ( ( m_capsuleRadius*m_scaleRadiusFactor ) - m_capsuleRadius ) * ( 1 - density ) ) ;

}



/* Main function */
template < class K >
typename SegmentQueryIntersectionOracle<K>::Object SegmentQueryIntersectionOracle<K>::intersection( const Segment_3& querySegment ) {

	if( m_debugOutput ) {
		m_iter++ ;
		if ( m_iter == 1 ) std::cout << std::endl ; // Just a small correction from the main() function...
		std::cout << "******* " << m_iter << " *******" << std::endl ;
	}
	
	/* Get the neighbors */
	FT radius = 0.0 ;
	std::vector< Point_3 > pts = intersection_getNeighbors( querySegment, radius ) ;
	if ( pts.size() < m_minInliers ) 
		return Object() ;
		
	/* Compute the scale of a set of points (using LKS + MSSE), if required to */
	FT distThres = intersection_getRansacThreshold( querySegment, pts ) ;
	if ( distThres < 0 ) 
		return Object() ;
	
	/* Compute the local LBQ model and test intersection */
	int trials = 0 ;
	return intersection_computeIntersection( querySegment, pts, distThres, radius, trials ) ;

}



template < class K >
std::vector< typename K::Point_3 > SegmentQueryIntersectionOracle<K>::intersection_getNeighbors( const typename K::Segment_3& querySegment, typename K::FT &newRadius ) {

	/* Compute points inside capsule */

	// Adapt radius to density?
	if ( m_computeDensity ) {
	
		/* Having into account a radial neighborhood */

		//std::vector< Point_3 > pts0 = this->getNeighborsOnCapsule( querySegment, m_capsuleRadius ) ;
	
		//if ( pts0.size() > 0 ) {

		//	// Update the capsule radius
		//	newRadius = this->getAdaptiveRadius( pts0 ) ;

		//	if( m_debugOutput ) {
		//		std::cout << "    Previous/New radius = " << m_capsuleRadius << " / " << newRadius << std::endl ;
		//	}
		//}
		//else {
		//	newRadius = m_capsuleRadius ;
		//}

		/* Having into account the density of the query segment nearest neighbor */
		newRadius = this->getAdaptiveRadiusNN( querySegment ) ;		
		m_newRadius = newRadius ;
		if( m_debugOutput ) {
			std::cout << "    Radius (min/max) = " << newRadius << "(" << m_capsuleRadius << "/" << m_capsuleRadius * m_scaleRadiusFactor << ")" << std::endl ;
		}

	}
	else {
		newRadius = m_capsuleRadius ;
		m_newRadius = newRadius ;
	}

	// Search with the final radius
	std::vector< Point_3 > pts = this->getNeighborsOnCapsule( querySegment, newRadius ) ;

	// if ( pts.size() < m_minInliers ) {
		//if( m_debugOutput ) {
			/*std::ofstream ofsi( "__segmentsWithoutOrWithFewNeigh.xyz", std::ios_base::app ) ;
			ofsi << queryCapSeg << std::endl ;
			ofsi.close() ;*/
		//}	

		// return Object() ; // It is done in the intersection function...
	//}

	if( m_debugOutput ) {
		/*std::ofstream ofs( "./_OUT/__neighPoints.xyz", std::ios_base::out ) ;
		for ( int i = 0; i < pts.size(); i++ ) {
			ofs << pts[i] << std::endl ;
		}*/
		// std::cout << "Press any key to continue" << std::endl ;
		// std::cin.get() ;
		// ofs.close() ;
		std::cout << "    Number of neighbors = " << pts.size() << std::endl ;

		// Show the squared size of the segment
		std::cout << "    Segment's squared length = " << querySegment.squared_length() << std::endl ;
	}

	return pts ;

}



template < class K >
typename K::FT SegmentQueryIntersectionOracle<K>::intersection_getRansacThreshold( const typename K::Segment_3 &querySegment, const std::vector< typename K::Point_3 > &pts ) {
	double distThres = 0 ;
	if ( m_ransacThresholdType == 1 ) {
		// Iterative scale estimation
	
		// Compute scale using more points than with the model estimation
		std::vector< Point_3 > ptsScale = this->getNeighborsOnCapsule( querySegment, m_capsuleRadius*m_scaleRadiusFactor ) ;

		FT scale = this->computeScale( ptsScale, Point_3() ) ;
		if( m_debugOutput ) {
			std::cout << "    Computed Scale (Now) = " << scale << std::endl ;
		}
		distThres = scale * m_ransacT ;
		if( m_debugOutput ) {
			std::cout << "    Distance threshold = " << distThres << std::endl ;
		}
		else {
			return -1.0 ; // Correct?
		}
	}
	else {
		// Check if there are pre-computed scales in the SplatsCreator object
		if ( m_hasPrecomputedScales ) {
			FT scale = this->getPrecomputedScale( pts ) ;	
			distThres = scale * m_ransacT ;

			if( m_debugOutput ) {
				std::cout << "    Computed Scale = " << scale << std::endl ;
			}

			distThres = scale * m_ransacT ;
			if( m_debugOutput ) {
				std::cout << "    Distance threshold = " << distThres << std::endl ;				
			}

		}
		else {

			if( m_debugOutput ) {
				std::cout << "    Fixed distance threshold = " << m_distThres << std::endl ;
			}

			// Otherwise, fixed distance threshold
			distThres = m_distThres ;
		}
	}

	return distThres ;
}



template < class K >
typename SegmentQueryIntersectionOracle<K>::Object SegmentQueryIntersectionOracle<K>::intersection_computeIntersection( const typename K::Segment_3 &querySegment,
																														const std::vector< typename K::Point_3 > &pts, 
																														const double &distThres,
																														const double &radius,
																														int &trialCount ) {

	// Check the trial number
	if ( trialCount > 2 ) {
		if ( m_debugOutput )
			std::cout << "    (Maximum number of trials exceeded)" << std::endl ;
		return Object() ;
	}

	trialCount++ ; // Update the number of trials

	std::vector< Point_3 > inliers ;
	std::vector< Point_3 > outliers ;
	std::vector< double > modelParam ;
	
			
	if ( m_degree == 1 ) {
		// Plane
		Plane_3 plane ;
		bool ok = this->computePlaneModel( pts, querySegment, distThres, inliers, outliers, modelParam, plane, radius ) ;

		if ( !ok ) {
			if( m_debugOutput ) {					
				std::cout << "    No model..." << std::endl ;

//				if ( m_debugVisualize ) {
//					// m_stepViewer->drawStep( querySegment, pts, inliers, Point_3( 99999999.9, 99999999.9, 99999999.9 ), lbq ) ;
//				}

				/*std::ostringstream noModelFileName ;
				noModelFileName << "./_OUT/__NoModel_" << m_iter << ".xyz" ;
				std::ofstream ofs( noModelFileName.str(), std::ios_base::app ) ;
				for ( int i = 0; i < pts.size(); i++ ) {
					ofs << pts[i] << std::endl ;
				}*/
				//// std::cout << "Press any key to continue" << std::endl ;
				//// std::cin.get() ;
				//ofs.close() ;
			}

			return Object() ;
		}
		else {
			Object int1 = CGAL::intersection( plane, querySegment ) ;
			if ( const Point_3 *ipoint1 = CGAL::object_cast< Point_3 >( &int1 ) ) {
				if ( intersectionHasConsensus( inliers, *ipoint1, m_newRadius*m_intersectionConsensusFactor, m_intersectionConsensusNumber ) ) {

					if( m_debugOutput ) {
						std::cout << "--> Intersection: " << ipoint1->x() << ", " << ipoint1->y() << ", " << ipoint1->z() << std::endl ;
						std::ofstream ofsi( "./_OUT/__intersect.xyz", std::ios_base::app ) ;
						ofsi << *ipoint1 << std::endl ;
						ofsi.close() ;
					}
//					if ( m_debugVisualize ) {
//						// m_stepViewer->drawStep( querySegment, pts, inliers, *ipoint1, lbq ) ;
//					}

					return int1 ;
				}
				else {

					if( m_debugOutput ) {
						std::cout << "--> No Intersection, hasn't got enough consensus" << std::endl ;
						// Second chance of intersection for this neighborhood, using the outliers of the model...
					
					}
//					if ( m_debugVisualize ) {
//						// m_stepViewer->drawStep( querySegment, pts, inliers, Point_3( 99999999.9, 99999999.9, 99999999.9 ), lbq ) ;
//					}


					if ( outliers.size() > m_minInliers ) {
						if( m_debugOutput ) {		
							std::cout << "    Check if there are more structures that might intersect the segment..." << std::endl ;							
						}
						// Second chance of intersection for this neighborhood, using the outliers of the model...
						return intersection_computeIntersection( querySegment, outliers, distThres, radius, trialCount ) ;
					}
					else {
						return Object() ;
					}

				}
			}

			return Object() ;
		}
	}
	else {
		LocalBivariateQuadric lbq ;
		bool ok = this->computeModel( pts, querySegment, distThres, inliers, outliers, modelParam, lbq, radius ) ;
		if ( !ok ) {

			if( m_debugOutput ) {					
				std::cout << "    No model..." << std::endl ;

//				if ( m_debugVisualize ) {
//					m_stepViewer->drawStep( querySegment, pts, inliers, Point_3( 99999999.9, 99999999.9, 99999999.9 ), lbq ) ;
//				}

				//std::ostringstream noModelFileName ;
				//noModelFileName << "./_OUT/__NoModel_" << m_iter << ".xyz" ;
				//std::ofstream ofs( noModelFileName.str(), std::ios_base::app ) ;
				//for ( int i = 0; i < pts.size(); i++ ) {
				//	ofs << pts[i] << std::endl ;
				//}
				//// std::cout << "Press any key to continue" << std::endl ;
				//// std::cin.get() ;
				//ofs.close() ;
			}

			return Object() ;
		}

		if( m_debugOutput ) {		
			std::cout << "    Inliers = " << inliers.size() << std::endl ;
			std::ofstream ofsIn( "./_OUT/__inlierPts.xyz", std::ios_base::out ) ;
			for ( int i = 0; i < inliers.size(); i++ ) {
				ofsIn << inliers[i] << std::endl ;
			}
			std::cout << "    Distance threshold = " << distThres << std::endl ;
			//// std::cout << "Press any key to continue" << std::endl ;
			//// std::cin.get() ;
			//ofsIn.close() ;
		}

		/* Compute the intersection between the segment and the lbq	*/
		Object int1, int2 ;
		if( lbq.intersection( querySegment, int1, int2 ) ) {

			if ( m_forceLBQOnSegment ) {
				if ( const Point_3 *ipoint2 = CGAL::object_cast< Point_3 >( &int2 ) ) {
					std::cerr << "[SegmentQueryIntersectionOracle] Error! Just one intersection should be detected along the Z axis of the fitting basis (the query segment)!" << std::endl ;
				}
			}

			if ( const Point_3 *ipoint1 = CGAL::object_cast< Point_3 >( &int1 ) ) {
				// Single intersection
											
				// Check the validity of the intersection point (must be supported by some neighbor points inside a radial neighborhood of distThres size
				// if ( intersectionHasConsensus( inliers, *ipoint1, distThres, 3 ) ) {
				// if ( intersectionHasConsensus( inliers, *ipoint1, m_capsuleRadius/2.0, 3 ) ) {
				if ( intersectionHasConsensus( inliers, *ipoint1, m_newRadius*m_intersectionConsensusFactor, m_intersectionConsensusNumber ) ) {

					if( m_debugOutput ) {
						std::cout << "--> Intersection: " << ipoint1->x() << ", " << ipoint1->y() << ", " << ipoint1->z() << std::endl ;
						std::ofstream ofsi( "./_OUT/__intersect.xyz", std::ios_base::app ) ;
						ofsi << *ipoint1 << std::endl ;
						ofsi.close() ;
					}
//					if ( m_debugVisualize ) {
//						m_stepViewer->drawStep( querySegment, pts, inliers, *ipoint1, lbq ) ;
//					}

					return int1 ;
				}
				else {

					if ( !m_forceLBQOnSegment ) {
						if ( const Point_3 *ipoint2 = CGAL::object_cast< Point_3 >( &int2 ) ) {
							// Still there is another chance
							if ( intersectionHasConsensus( inliers, *ipoint2, m_newRadius*m_intersectionConsensusFactor, m_intersectionConsensusNumber ) ) {

								if( m_debugOutput ) {
									std::cout << "--> Intersection (2): " << ipoint2->x() << ", " << ipoint2->y() << ", " << ipoint2->z() << std::endl ;
									std::ofstream ofsi( "./_OUT/__intersect.xyz", std::ios_base::app ) ;
									ofsi << *ipoint2 << std::endl ;
									ofsi.close() ;
								}
//								if ( m_debugVisualize ) {
//									m_stepViewer->drawStep( querySegment, pts, inliers, *ipoint2, lbq ) ;
//								}

								return int2 ;

							}
						}
					}

					if( m_debugOutput ) {
						std::cout << "--> No Intersection, hasn't got enough consensus" << std::endl ;
						// Second chance of intersection for this neighborhood, using the outliers of the model...
					
					}
//					if ( m_debugVisualize ) {
//						m_stepViewer->drawStep( querySegment, pts, inliers, Point_3( 99999999.9, 99999999.9, 99999999.9 ), lbq ) ;
//					}


					if ( outliers.size() > m_minInliers ) {
						if( m_debugOutput ) {		
							std::cout << "    Check if there are more structures that might intersect the segment..." << std::endl ;							
						}
						// Second chance of intersection for this neighborhood, using the outliers of the model...
						return intersection_computeIntersection( querySegment, outliers, distThres, radius, trialCount ) ;
					}
					else {
						return Object() ;
					}

				}

			}
			else {
				// Should never happen!
				std::cout << "This should never happen!" << std::endl ;
				return Object() ;
			}
		}
		else {
			if( m_debugOutput ) {		
				std::cout << "--> No Intersection" << std::endl ;							
			}
//			if ( m_debugVisualize ) {
//				m_stepViewer->drawStep( querySegment, pts, inliers, Point_3( 99999999.9, 99999999.9, 99999999.9 ), lbq ) ;
//			}

			if ( outliers.size() > m_minInliers ) {
				if( m_debugOutput ) {		
					std::cout << "    Check if there are more structures that might intersect the segment..." << std::endl ;							
				}
				// Second chance of intersection for this neighborhood, using the outliers of the model...
				return intersection_computeIntersection( querySegment, outliers, distThres, radius, trialCount ) ;
			}
			else {
				return Object() ;
			}
		}
	}

	return Object() ;
}



#endif // SEGMENTQUERYINTERSECTIONORACLE_H