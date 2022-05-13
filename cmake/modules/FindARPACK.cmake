# Simple cmake finder for the ARPACK library (only sets the ARPACK_LIBRARIES variable, no includes as we use ARPACK++ in fact)
# Author: Ricard Campos (rcampos@eia.udg.edu)

SET(ARPACK_FOUND FALSE)

# Find libraries
FIND_LIBRARY(ARPACK_LIBRARIES NAMES arpack libarpack
        PATHS
        ${ARPACK_LIB_DIR}
        $ENV{ARPACK_LIB_DIR}
        /usr/local/lib
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        NO_DEFAULT_PATH)

if(ARPACK_LIBRARIES)
    SET(ARPACK_FOUND TRUE)
endif(ARPACK_LIBRARIES)

IF(ARPACK_FOUND)
    IF(NOT ARPACK_FIND_QUIETLY)
        MESSAGE(STATUS "Found ARPACK: libraries at ${ARPACK_LIBRARIES}")
    ENDIF(NOT ARPACK_FIND_QUIETLY)
ELSE(ARPACK_FOUND)
    IF(ARPACK_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "ARPACK not found")
    ENDIF(ARPACK_FIND_REQUIRED)
ENDIF(ARPACK_FOUND)