add_executable(mesh_splats_ncut SplatsDistanceNCut.cpp)
add_definitions(-fpermissive)
target_link_libraries(mesh_splats_ncut ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${ARPACKPP_LIBRARIES} ${ARPACK_LIBRARIES})
install(TARGETS mesh_splats_ncut DESTINATION bin)