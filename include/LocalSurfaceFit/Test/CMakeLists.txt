# CMake compilation instructions for TestLocalBivariateQuadric project
#
# Required libraries:
#	- CGAL
# 	- Boost
#	- CImg
#
# Note: Related projects are assumed to be relative to the current folder, 
#		change the include paths if this is not the case.
# 		Since some OS do not have default paths for libraries/includes, the 
#		use of some (environmental) variables is required in this case to point 
#		to this paths. This script will look for the following (environmental) 
#		variables for each of the required libraries: 
#			- CIMG_INC: Path where the cimg.h include file is located.
# 
# Author: Ricard Campos (rcampos@eia.udg.edu)
# Date: 19/03/2013
#

project( TestLocalBivariateQuadric ) 

cmake_minimum_required( VERSION 2.6 FATAL_ERROR )

# --> Check include directories of external libraries
# CImg (environmental)
if( DEFINED ENV{CIMG_INC} )
	include_directories( BEFORE $ENV{CIMG_INC} )
endif()
# CImg
if( DEFINED CIMG_INC )
	include_directories( BEFORE ${CIMG_INC} )
endif()
# --> Check libraries directories of external libraries

# --> Check required packages
find_package(CGAL COMPONENTS Core)
if( CGAL_FOUND )
	include( ${CGAL_USE_FILE} )
endif()

find_package(Boost COMPONENTS thread filesystem )
if ( Boost_FOUND )
	link_directories ( ${Boost_LIBRARY_DIRS} )
	include_directories ( ${Boost_INCLUDE_DIRS} )
endif()

# Find Eigen
find_package(Eigen3 3.0.91) #(requires 3.1.0-alpha1 or greater)
if ( EIGEN3_FOUND )
	include( ${EIGEN3_USE_FILE} )
endif()

# --> Include directories required
# Related projects
include_directories( BEFORE .. )
include_directories( BEFORE ../../../RobustStatistics )

# Workaround for llebre.udg.edu (2 versions of cgal are installed in default paths)
if( CMAKE_COMPILER_IS_GNUCXX )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /usr/include" )
endif() 

if( CGAL_FOUND AND Boost_FOUND AND EIGEN3_FOUND ) 

	link_directories ( ${CGAL_LIBRARY_DIRS} )
	include_directories ( ${CGAL_INCLUDE_DIRS} )
	
	link_directories ( ${Boost_LIBRARY_DIRS} )
	include_directories ( ${Boost_INCLUDE_DIRS} )
	
	# --> Compile
	ADD_EXECUTABLE( TestLocalBivariateQuadric  			TestLocalBivariateQuadric.cpp )
	ADD_EXECUTABLE( TestLocalBivariateQuadricRANSAC  	TestLocalBivariateQuadricRANSAC.cpp
														../LBQModelEstimator.cpp )
	
	# --> Link the executable
	target_link_libraries( TestLocalBivariateQuadric ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES})				
	target_link_libraries( TestLocalBivariateQuadricRANSAC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES})				
	
else( CGAL_FOUND AND Boost_FOUND AND EIGEN3_FOUND )

	message( STATUS "Missing dependencies!" )
	
endif()



