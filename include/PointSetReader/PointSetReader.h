#ifndef POINTSET_READER_H
#define POINTSET_READER_H

#include <iostream>
#include <fstream>

template< typename PointType, typename NormalType >
class PointSetReader {

public:
	/** 
	 * Reads a pointset file
	 */
	// TODO: No error check is performed!
	static bool readPointsetFile(	const char *name, 
									std::vector< PointType >& points ) ;

	/**
	 * Reads an oriented pointset file (points with associated normals) 
	 */
	//TODO: No error check is performed!
	static bool readOrientedPointsetFile(	const char *name, 
									std::vector< PointType >& points, 
									std::vector< NormalType >& normals ) ;

	/**
	* Reads an oriented pointset file with sizes (points with associated normals and sizes)
	*/ 
	//TODO: No error check is performed!
	static bool readOrientedPointsetFileWithScores( const char *name, 
													std::vector< PointType >& points, 
													std::vector< NormalType >& normals,
													std::vector< double >& sizes ) ;					

} ;


/*** Implementation ***/

template< typename PointType, typename NormalType >
bool PointSetReader< PointType, NormalType >::readPointsetFile(	const char *name, 
																std::vector< PointType >& points ) {
									
	std::ifstream file( name, std::ifstream::in ) ;
	if ( !file ) {
		return false ;
	}
		
	while( !file.eof() ) {
		// Read a line
		std::string str ; 
		getline( file, str ) ;
		std::istringstream iss( str ) ;
		double point[3] ;

		iss >> point[0] >> point[1] >> point[2] ;
		
		// Debug
		//std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl ;
		
		points.push_back( PointType( point[0], point[1], point[2] ) ) ;		
	}

	file.close() ;
	return true ;
}



// Reads an oriented pointset file 
// Note that no sizes are readed, this function is different from the one in OrientedPointsetRemeshing.cpp!
// TODO: No error check is performed!
template< typename PointType, typename NormalType >
bool PointSetReader< PointType, NormalType >::readOrientedPointsetFile(	const char *name, 
																		std::vector< PointType >& points, 
																		std::vector< NormalType >& normals ) {
									
	std::ifstream file( name, std::ifstream::in ) ;
	if ( !file ) {
		return false ;
	}
		
	while( !file.eof() ) {
		// Read a line
		std::string str ; 
		getline( file, str ) ;
		std::istringstream iss( str ) ;
		double point[3], normal[3] ;
		
		iss >> point[0] >> point[1] >> point[2] >> normal[0] >> normal[1] >> normal[2] ;
		
		// Debug
		//cout << point[0] << " " << point[1] << " " << point[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << endl ;
		
		points.push_back( PointType( point[0], point[1], point[2] ) ) ;
		normals.push_back( NormalType( normal[0], normal[1], normal[2] ) ) ;		
	}

	file.close() ;
	return true ;

}



template< typename PointType, typename NormalType >
bool PointSetReader< PointType, NormalType >::readOrientedPointsetFileWithScores( const char *name, 
																				 std::vector< PointType >& points, 
																				 std::vector< NormalType >& normals,
																				 std::vector< double >& scores ) {
									
	std::ifstream file( name, std::ifstream::in ) ;
	if ( !file ) {
		return false ;
	}
		
	while( !file.eof() ) {
		// Read a line
		std::string str ; 
		getline( file, str ) ;
		std::istringstream iss( str ) ;
		double point[3], normal[3] ;
		double score ;

		iss >> point[0] >> point[1] >> point[2] >> normal[0] >> normal[1] >> normal[2] >> score ;
		
		// Debug
		//cout << point[0] << " " << point[1] << " " << point[2] << " " << normal[0] << " " << normal[1] << " " << normal[2] << endl ;
		
		points.push_back( PointType( point[0], point[1], point[2] ) ) ;
		normals.push_back( NormalType( normal[0], normal[1], normal[2] ) ) ;		
		scores.push_back( score ) ;
	}

	file.close() ;
	return true ;

}


#endif // POINTSET_READER_H
