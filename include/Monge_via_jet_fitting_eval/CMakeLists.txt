# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Monge_via_jet_fitting_eval )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories( BEFORE . ) # Needed to redefine Monge_form
  
  # use either Eigen or BLAS/LAPACK
  find_package(Eigen3 3.0.91) #(requires 3.1.0-alpha1 or greater)
  if (NOT EIGEN3_FOUND)
    find_package(LAPACK)
    if(LAPACK_FOUND)
      include( ${LAPACK_USE_FILE} )
    endif(LAPACK_FOUND)
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()
  
  # --> Check include directories of external libraries
  # Levmar (environmental)
  if( DEFINED ENV{LEVMAR_INC} )
	include_directories( BEFORE $ENV{LEVMAR_INC} )
  endif()
  # Levmar
  if( DEFINED LEVMAR_INC )
	include_directories( BEFORE ${LEVMAR_INC} )
  endif()
  # --> Check libraries directories of external libraries
  # Levmar (environmental)
  if( DEFINED ENV{LEVMAR_LIB} )
	link_directories ( $ENV{LEVMAR_LIB} )
  endif()
  # Levmar
  if( DEFINED LEVMAR_LIB )
	link_directories ( ${LEVMAR_LIB} )
  endif()
  # Lapack (environmental)
  if( DEFINED ENV{LAPACK_LIB} )
	link_directories ( $ENV{LAPACK_LIB} )
  endif()
  # Lapack
  if( DEFINED LAPACK_LIB )
	link_directories ( ${LAPACK_LIB} )
  endif()
    
  if(EIGEN3_FOUND OR LAPACK_FOUND)
    # Link with Boost.ProgramOptions (optional)
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
      else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      endif()
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      if ( NOT CGAL_AUTO_LINK_ENABLED )
        link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
      endif()
    endif()

	# --> Compile projection test
    #create_single_source_cgal_program( "Monge_via_jet_fitting_eval.cpp" )    
	
	# --> Compile intersection test
	ADD_EXECUTABLE( Monge_via_jet_fitting_eval Monge_via_jet_fitting_eval.cpp )
	ADD_EXECUTABLE( Monge_via_jet_fitting_intersection Monge_via_jet_fitting_intersection.cpp )
	
	# --> Link the executable
	if ( UNIX )
		target_link_libraries( Monge_via_jet_fitting_eval ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES} levmar lapack blas X11 pthread )
		target_link_libraries( Monge_via_jet_fitting_intersection ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES} levmar lapack blas X11 pthread )
	endif()
	
	if( WIN32 )
		target_link_libraries( Monge_via_jet_fitting_eval ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES} levmar lapack blas libf2c )
		target_link_libraries( Monge_via_jet_fitting_intersection ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_LIBRARIES} levmar lapack blas libf2c )
	endif()

  else(EIGEN3_FOUND OR LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

