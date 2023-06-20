
add_pfunit_ctest (neptuneClass_tests
    TEST_SOURCES test_neptuneClass.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

if (NOT SKIP_MSIS_2)
    add_pfunit_ctest (atmosphere_tests_msis2
        TEST_SOURCES test_atmosphere_msis2.pf
        LINK_LIBRARIES neptune ${LIBSLAM}
        )
endif()

add_pfunit_ctest (atmosphere_tests_nrlmsis00
    TEST_SOURCES test_atmosphere_nrlmsis00.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

    add_pfunit_ctest (atmosphere_tests_JB2008
    TEST_SOURCES test_atmosphere_JB2008.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )
    add_pfunit_ctest (atmosphere_tests_Exponential
    TEST_SOURCES test_atmosphere_Exponential.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )
    add_pfunit_ctest (solarsystem_tests
    TEST_SOURCES test_solarsystem.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

    add_pfunit_ctest (Thirdbody_tests
    TEST_SOURCES test_thirdbody.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

    add_pfunit_ctest (Radiation_tests
    TEST_SOURCES test_radiation.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

    # add_pfunit_ctest (gravity_tests
    # TEST_SOURCES test_gravity.pf
    # LINK_LIBRARIES neptune ${LIBSLAM}
    # )

add_pfunit_ctest (tides_tests
    TEST_SOURCES test_tides.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )