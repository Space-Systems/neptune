
add_pfunit_ctest (neptuneClass_tests
    TEST_SOURCES test_neptuneClass.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

add_pfunit_ctest (atmosphere_tests
    TEST_SOURCES test_atmosphere.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )

add_pfunit_ctest (tides_tests
    TEST_SOURCES test_tides.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )