add_executable(mesh_splats mesh_splats.cpp)
target_link_libraries(mesh_splats ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})