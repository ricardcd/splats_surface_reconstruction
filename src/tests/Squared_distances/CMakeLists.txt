# CMake compilation instructions for SquaredDistances project

project( SquaredDistances ) 

cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# --> Select build type
SET( CMAKE_BUILD_TYPE Release )

# Workaround for llebre.udg.edu (2 versions of cgal are installed in default paths)
#if( CMAKE_COMPILER_IS_GNUCXX )
#  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /usr/include" )
#endif()

# --> Check required packages
find_package(CGAL COMPONENTS Core)
if( CGAL_FOUND )
	include( ${CGAL_USE_FILE} )
endif()

include_directories( BEFORE .. )

if( CGAL_FOUND ) 
	link_directories ( ${CGAL_LIBRARY_DIRS} )
	include_directories ( ${CGAL_INCLUDE_DIRS} )
	
	# --> Compile
	ADD_EXECUTABLE(test_squared_distances test_squared_distances.cpp)
	
	# --> Link the executable
	target_link_libraries( test_squared_distances ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

