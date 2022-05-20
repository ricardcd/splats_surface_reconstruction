add_executable(orient_splats OrientSplats.cpp)
target_link_libraries(orient_splats ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS orient_splats DESTINATION bin)