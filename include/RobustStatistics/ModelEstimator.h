#ifndef MODEL_ESTIMATOR_H
#define MODEL_ESTIMATOR_H

#include <vector>

template< class DataType, class ParametersType >
class ModelEstimator {
public:	

	ModelEstimator( unsigned int minSamples ) : m_minSamples( minSamples ) {} 

	virtual void fit( const std::vector< typename std::vector< DataType >::const_pointer > &data,
					  std::vector< ParametersType > &parameters ) const = 0 ;

	virtual void leastSquaresFit( const std::vector< typename std::vector< DataType >::const_pointer > &data,
								  std::vector< ParametersType > &parameters ) const = 0 ;
	
	// Used in RANSAC method
	virtual bool compatible( const std::vector< ParametersType > &parameters, 
							 const DataType &data ) const = 0 ;

	// Used in LeastKthSquares method
	virtual std::vector< double > sqResiduals(  const std::vector< ParametersType > &parameters, 
											  const std::vector< DataType > &data ) const = 0 ;

	virtual bool degenerate( const std::vector< typename std::vector< DataType >::const_pointer > &data ) const = 0 ;

	unsigned int minSamples() const { return m_minSamples ; }

private:
	unsigned int m_minSamples ;	

} ;

#endif // MODEL_ESTIMATOR_H