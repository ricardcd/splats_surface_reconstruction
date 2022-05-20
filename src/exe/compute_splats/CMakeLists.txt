add_executable(compute_splats ComputeSplats.cpp
               ${PROJECT_SOURCE_DIR}/src/RobustStatistics/Models/FitLocalSmoothSurface/LocalSmoothSurfaceModelEstimator.cpp)
target_link_libraries(compute_splats ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS compute_splats DESTINATION bin)