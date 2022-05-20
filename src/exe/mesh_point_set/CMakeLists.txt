add_executable(mesh_point_set PointSetMesher.cpp
               ${PROJECT_SOURCE_DIR}/src/RobustStatistics/Models/FitPlane/PlaneModelEstimator.cpp
               ${PROJECT_SOURCE_DIR}/src/LocalSurfaceFit/LBQModelEstimator.cpp
               ${PROJECT_SOURCE_DIR}/src/LocalSurfaceFit/WLBQModelEstimator.cpp
               ${PROJECT_SOURCE_DIR}/src/LocalSurfaceFit/WLBQSegmentModelEstimator.cpp)
target_link_libraries(mesh_point_set ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS mesh_point_set DESTINATION bin)