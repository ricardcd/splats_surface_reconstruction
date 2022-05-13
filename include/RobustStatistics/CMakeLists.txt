cmake_minimum_required(VERSION 2.8)

PROJECT( RANSAC_EXAMPLE )

INCLUDE_DIRECTORIES (
  ${RANSAC_EXAMPLE_SOURCE_DIR}
)

SET(RANSAC_EXAMPLE_HDRS
    Ransac.h
	ModelEstimator.h    
    Models/FitLine2D/LineModelEstimator.h    
    )

SET(RANSAC_EXAMPLE_SRCS
    Models/FitLine2D/LineModelEstimator.cpp
    RansacFitLineTest.cpp
    )

# --> Check required packages
find_package(CGAL COMPONENTS Core)
if( CGAL_FOUND )
	include( ${CGAL_USE_FILE} )
endif()	
	
ADD_EXECUTABLE( RansacFitLineTest ${RANSAC_EXAMPLE_HDRS} ${RANSAC_EXAMPLE_SRCS} )
target_link_libraries( RansacFitLineTest ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
