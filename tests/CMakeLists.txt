
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

add_pfunit_ctest (tides_tests
    TEST_SOURCES test_tides.pf
    LINK_LIBRARIES neptune ${LIBSLAM}
    )
    