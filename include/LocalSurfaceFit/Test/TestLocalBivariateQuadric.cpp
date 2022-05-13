/* Includes */
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/centroid.h>
// Project-related
#include "LocalBivariateQuadric.h"
// CImg (just for parsing input parameters)
#include "CImg.h"
// Std
#include <algorithm>
#include <fstream>
using namespace std ;

// Main Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::FT				FT ;
typedef K::Point_3			Point_3 ;
typedef K::Point_2			Point_2 ;
typedef K::Vector_3			Vector_3 ;
typedef Fitting::LocalBivariateQuadric<K>		LocalBivariateQuadric ;



/* Main function */
int main ( int argc, char **argv) {

	/* Parse input parameters */
	cimg_usage( "Simple program to test the fitting of local quadrics." ) ;

	const char*			inputFile = cimg_option(		"-i",				(char*)0,	"Input Point Set file." ) ;
	const char*			outputFile = cimg_option(		"-o",				(char*)0,	"Output file." ) ;
	const double		h			= cimg_option(		"-h",				1.0,		"h (spacing) weight parameter." ) ;

	

	// Read input points
	std::ifstream in( inputFile ) ; 
	if ( !in ) {
		std::cerr << "Cannot open input!" << std::endl ;
	}

	std::vector< Point_3 > points ;
	bool readed = read_xyz_points( in, std::back_inserter( points ) ) ;
	if ( !readed ) {
		cerr << "Error reading the file..." << endl ;
		return -1 ;
	}

	// Compute the centroid of the XY plane (where the fit will be centered)
	std::vector< Point_2 > pts2D ;
	for( std::vector< Point_3 >::iterator itP = points.begin(); itP != points.end(); ++itP ) {
		pts2D.push_back( Point_2( itP->x(), itP->y() ) ) ;
	}
	Point_2 centroid = CGAL::centroid( pts2D.begin(), pts2D.end() ) ;

	/* Compute the LS fit */
	Point_3 pt( centroid.x(), centroid.y(), 0 ) ;
	Vector_3 vec( 0, 0, 1 ) ;
	LocalBivariateQuadric lbq( pt, vec ) ;
	lbq.fit( points ) ;

	// Show results
	std::vector< FT > coeffs = lbq.coefficients() ;
	cout << "Computed Coefficients of the equation f(x,y) = Ax^2 + By^2 + Cx + Dy + Exy + F:" << endl ;
	cout << "Least-Squares:" << endl ;
	cout << "  A = " << coeffs[0] << endl ;
	cout << "  B = " << coeffs[1] << endl ;
	cout << "  C = " << coeffs[2] << endl ;
	cout << "  D = " << coeffs[3] << endl ;
	cout << "  E = " << coeffs[4] << endl ;
	cout << "  F = " << coeffs[5] << endl ;

	/* Compute the weighted LS fit */
	std::vector<FT> weights ;
	int i = 0 ;
	for ( std::vector< Point_2 >::iterator itP2 = pts2D.begin(); itP2 != pts2D.end(); ++itP2, i++ ) {
		FT w = LocalBivariateQuadric::weightGaussian( CGAL::sqrt( ( itP2->x() - centroid.x() ) * ( itP2->x() - centroid.x() ) + 
																	 ( itP2->y() - centroid.y() ) * ( itP2->y() - centroid.y() ) ), h ) ;

		weights.push_back( w ) ;
	}
	lbq.weightedFit( points, weights ) ;
			
	// Show results
	std::vector< FT > weightedCoeffs = lbq.coefficients() ;
	cout << "Weighted Least-Squares at centroid:" << endl ;
	cout <<	"  Centroid = " << centroid << endl ;
	cout <<	"  h = " << h << endl ;
	cout << "  A = " << weightedCoeffs[0] << endl ;
	cout << "  B = " << weightedCoeffs[1] << endl ;
	cout << "  C = " << weightedCoeffs[2] << endl ;
	cout << "  D = " << weightedCoeffs[3] << endl ;
	cout << "  E = " << weightedCoeffs[4] << endl ;
	cout << "  F = " << weightedCoeffs[5] << endl ;

	// First line --> Least Squares
	// Second line --> Weighted Least Squares
	ofstream oCoeffs( outputFile ) ;
	oCoeffs << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3] << " " << coeffs[4] << " " << coeffs[5] << endl ;
	oCoeffs << weightedCoeffs[0] << " " << weightedCoeffs[1] << " " << weightedCoeffs[2] << " " << weightedCoeffs[3] << " " << weightedCoeffs[4] << " " << weightedCoeffs[5] ;
	oCoeffs.close() ;

	return 1 ;

} 
