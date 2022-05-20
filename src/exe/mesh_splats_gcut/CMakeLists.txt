add_executable(mesh_splats_gcut SplatsDistanceGCut.cpp
        ${PROJECT_SOURCE_DIR}/src/RobustStatistics/Models/FitCentroid/CentroidModelEstimator.cpp
        ${PROJECT_SOURCE_DIR}/3rd_party/maxflow-v3.02/graph.cpp
        ${PROJECT_SOURCE_DIR}/3rd_party/maxflow-v3.02/maxflow.cpp)
target_link_libraries(mesh_splats_gcut ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS mesh_splats_gcut DESTINATION bin)