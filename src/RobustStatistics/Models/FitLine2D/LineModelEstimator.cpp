#include "RobustStatistics/Models/FitLine2D/LineModelEstimator.h"
#include <math.h>

LineModelEstimator::LineModelEstimator( double distThres ) : ModelEstimator< Point_2, double >( 2 ), m_distanceThreshold( distThres ) {} 



void LineModelEstimator::fit( const std::vector< std::vector<Point_2>::const_pointer > &data,
							  std::vector< double > &parameters ) const {

	parameters.clear() ;

	if ( data.size() < this->minSamples() ) {
		// Not enough samples to build a model...
		return ;
	}

	/*Line_2 l( data[0], data[1] ) ;

	parameters.push_back( l.a() ) ;
	parameters.push_back( l.b() ) ;
	parameters.push_back( l.c() ) ;*/


	double nx = data[1]->y() - data[0]->y();
	double ny = data[0]->x() - data[1]->x();
	double norm = sqrt(nx*nx + ny*ny);
	
	parameters.push_back(nx/norm);
	parameters.push_back(ny/norm);
	parameters.push_back(data[0]->x());
	parameters.push_back(data[0]->y());		
}



void LineModelEstimator::leastSquaresFit( const std::vector< std::vector<Point_2>::const_pointer > &data,
										  std::vector< double > &parameters ) const {

	/*Line_2 line_fit ;
	CGAL::linear_least_squares_fitting_2( data.begin(), data.end(), line_fit, CGAL::Dimension_tag<0> ) ;*/
	
	double meanX, meanY, nx, ny, norm;
	double covMat11, covMat12, covMat21, covMat22; // The entries of the symmetric covarinace matrix
	int i, dataSize = (int)data.size();

	parameters.clear();
	if(data.size()<this->minSamples())
		return;

	meanX = meanY = 0.0;
	covMat11 = covMat12 = covMat21 = covMat22 = 0;
	for(i=0; i<dataSize; i++) {
		meanX +=data[i]->x();
		meanY +=data[i]->y();

		covMat11	+=data[i]->x() * data[i]->x();
		covMat12	+=data[i]->x() * data[i]->y();
		covMat22	+=data[i]->y() * data[i]->y();
	}

	meanX/=dataSize;
	meanY/=dataSize;

	covMat11 -= dataSize*meanX*meanX;
	covMat12 -= dataSize*meanX*meanY;
	covMat22 -= dataSize*meanY*meanY;
	covMat21 = covMat12;

	if(covMat11<1e-12) {
		nx = 1.0;
	  ny = 0.0;
	}
	else {	    //lamda1 is the largest eigen-value of the covariance matrix 
	           //and is used to compute the eigne-vector corresponding to the smallest
	           //eigenvalue, which isn't computed explicitly.
		double lamda1 = (covMat11 + covMat22 + sqrt((covMat11-covMat22)*(covMat11-covMat22) + 4*covMat12*covMat12)) / 2.0;
		nx = -covMat12;
		ny = lamda1 - covMat22;
		norm = sqrt(nx*nx + ny*ny);
		nx/=norm;
		ny/=norm;
	}
	parameters.push_back(nx);
	parameters.push_back(ny);
	parameters.push_back(meanX);
	parameters.push_back(meanY);

}
	


bool LineModelEstimator::compatible( const std::vector< double > &parameters, 
									 const Point_2 &data ) const {

	double distance = parameters[0]*(data.x()-parameters[2]) + parameters[1]*(data.y()-parameters[3]) ; 
	return distance < m_distanceThreshold ;

}


std::vector< double > LineModelEstimator::sqResiduals(  const std::vector< double > &parameters, 
													  const std::vector< Point_2 > &data ) const {

	std::vector< Point_2 >::const_iterator it ;
	std::vector< double > residuals ;
	for( it = data.begin(); it != data.end(); ++it ) {
		double distance = parameters[0]*(it->x()-parameters[2]) + parameters[1]*(it->y()-parameters[3]) ; 
		residuals.push_back( abs( distance*distance ) ) ;
	}

	return residuals ;
}


bool LineModelEstimator::degenerate( const std::vector< std::vector<Point_2>::const_pointer > &data ) const {
	
	Line_2 l( *data[0], *data[1] ) ;
	return l.is_degenerate() ;

}