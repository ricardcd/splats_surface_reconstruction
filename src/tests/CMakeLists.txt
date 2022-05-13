if(NOT SPLATS_BUILD_TESTS)
    set(SPLATS_BUILD_TESTS True)
endif()

if (SPLATS_BUILD_TESTS)
    add_subdirectory(Monge_via_jet_fitting_eval)
    add_subdirectory(Fuzzy_capsule_3)
    add_subdirectory(Squared_distances)
endif()