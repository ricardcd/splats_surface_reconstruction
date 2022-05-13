project( CGAL_Octree )

cmake_minimum_required( VERSION 2.6 FATAL_ERROR )

# --> Check required packages
find_package(CGAL COMPONENTS Core)
if( CGAL_FOUND )
	include( ${CGAL_USE_FILE} )
	link_directories ( ${CGAL_LIBRARY_DIRS} )
	include_directories ( ${CGAL_INCLUDE_DIRS} )
endif()

if( CGAL_FOUND )
	add_executable( OctreeTestGetMean OctreeTestGetMean.cpp )
	target_link_libraries( OctreeTestGetMean ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
endif()