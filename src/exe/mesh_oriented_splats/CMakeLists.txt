add_executable(mesh_oriented_splats OrientedSplatsMesher.cpp)
target_link_libraries(mesh_oriented_splats ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS mesh_oriented_splats DESTINATION bin)