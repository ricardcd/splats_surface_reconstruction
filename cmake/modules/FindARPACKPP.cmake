# Simple cmake finder for the ARPACK++ library
# Author: Ricard Campos (rcampos@eia.udg.edu)

SET(ARPACKPP_FOUND FALSE)

# Find the header files
FIND_PATH(ARPACKPP_INCLUDE_DIR arrssym.h
    ${ARPACKPP_INC_DIR}
    $ENV{ARPACKPP_INC_DIR}
    /usr/local/include
    /usr/include
    /usr/include/arpack++
    /usr/include/arpackpp
    /opt/local/include
    /sw/local/include
    /sw/include
    NO_DEFAULT_PATH)

# Find libraries
FIND_LIBRARY(ARPACKPP_LIBRARIES NAMES arpack++ libarpack++ arpackpp libarpackpp
        PATHS
        ${ARPACKPP_LIB_DIR}
        $ENV{ARPACKPP_LIB_DIR}
        /usr/local/lib
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        NO_DEFAULT_PATH)

if(ARPACKPP_LIBRARIES)
    SET(ARPACKPP_FOUND TRUE)
endif(ARPACKPP_LIBRARIES)

IF(ARPACKPP_FOUND)
    IF(NOT ARPACKPP_FIND_QUIETLY)
        MESSAGE(STATUS "Found ARPACK++: libraries at ${ARPACKPP_LIBRARIES}")
        MESSAGE(STATUS "Found ARPACK++: includes at ${ARPACKPP_INCLUDE_DIR}")
    ENDIF(NOT ARPACKPP_FIND_QUIETLY)
ELSE(ARPACKPP_FOUND)
    IF(ARPACKPP_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "ARPACK++ not found")
    ENDIF(ARPACKPP_FIND_REQUIRED)
ENDIF(ARPACKPP_FOUND)