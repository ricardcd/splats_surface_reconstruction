# Simple cmake finder for the levmar library
# Author: Ricard Campos (rcampos@eia.udg.edu)

# Find the header files
FIND_PATH(LEVMAR_INCLUDE_DIR levmar.h
        ${LEVMAR_INC_DIR}
        $ENV{LEVMAR_INC_DIR}
        /usr/local/include
        /usr/include
        /opt/local/include
        /sw/local/include
        /sw/include
        NO_DEFAULT_PATH)

FIND_FILE(LEVMAR_LIBRARY liblevmar.a
        ${LEVMAR_LIB_DIR}
        $ENV{LEVMAR_LIB_DIR}
        /usr/local/lib
        /usr/lib
        NO_DEFAULT_PATH)

