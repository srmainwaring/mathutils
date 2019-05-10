
add_custom_target(mathutils_tests)

set(TESTS
        test_unit
        test_angles
        test_interp1d
        test_lookup_table
        test_interp2d
        test_vector2d
        test_vector3d
        test_transform3d
        test_matrixMN
        test_integration1d
        test_fft
        test_convolution
        test_spline
        test_vectorN
        test_mkl
        )


foreach(test ${TESTS})

    message(STATUS "...add ${test}")

    add_executable(${test} ${test}.cpp)



    target_link_libraries(${test} MathUtils::MathUtils Eigen3::Eigen)

    if (${MATHUTILS_MKL})
        target_link_libraries(${test} ${MKL_LIBRARIES} ${IOMP5_LIBRARY})
        target_include_directories(${test} PRIVATE ${MKL_INCLUDE_DIRS})
        target_compile_definitions(${test} PRIVATE WITH_MKL)
    endif()

    add_dependencies(mathutils_tests ${test})

endforeach()
