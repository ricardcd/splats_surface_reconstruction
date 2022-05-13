add_executable(Monge_via_jet_fitting_eval Monge_via_jet_fitting_eval.cpp)
add_executable(Monge_via_jet_fitting_intersection Monge_via_jet_fitting_intersection.cpp)

target_link_libraries(Monge_via_jet_fitting_eval ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}  ${BOOST_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(Monge_via_jet_fitting_intersection ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}  ${BOOST_LIBRARIES} ${LEVMAR_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
