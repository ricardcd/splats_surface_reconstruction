project(splats_surface_reconstruction)
cmake_minimum_required( VERSION 2.6 FATAL_ERROR )

# Extent the module path (required for our FindLevmar.cmake)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# -----------------------------------------------------------------------------
# Additional compiler options
# -----------------------------------------------------------------------------
#if ( CMAKE_COMPILER_IS_GNUCXX )
#    ADD_DEFINITIONS(-O3 -frounding-math)
#endif()

add_definitions(-std=c++11 -frounding-math)

# -----------------------------------------------------------------------------
# Required packages
# -----------------------------------------------------------------------------
# CGAL
find_package(CGAL)
include( ${CGAL_USE_FILE} )

# Boost
find_package(Boost 1.52 EXACT COMPONENTS program_options filesystem system thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Eigen
find_package(Eigen3) #(requires 3.1.0-alpha1 or greater)
include(${EIGEN3_USE_FILE})

# --> Check include directories of external libraries
find_package(Levmar REQUIRED)
include_directories(${LEVMAR_INCLUDE_DIR})

# Lapack/Blas (used by levmar)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

if (EIGEN_FOUND)
    # If found, prefer Eigen for CGAL
    include( ${EIGEN_USE_FILE} )
    add_definitions(-DCGAL_EIGEN_ENABLED)
else() # Lapack must be found here, as it is set as required above
    include( ${LAPACK_USE_FILE} )
    add_definitions(-DCGAL_LAPACK_ENABLED)
endif()

# ARPACK
find_package(ARPACK REQUIRED)

# ARPACK++
find_package(ARPACKPP REQUIRED)
include_directories(${ARPACKPP_INCLUDE_DIR})

# -----------------------------------------------------------------------------
# Project-specific options
# -----------------------------------------------------------------------------

# Add the header-only maxflow package (3rd party)
include_directories(3rd_party/maxflow-v3.02)

# Other includes
include_directories(${CMAKE_SOURCE_DIR}/include)

# Follow down the source path
add_subdirectory(src)


