

include_directories("..")  # TODO: voir a recuperer la variable qui va bien...
#message(${MATHUTILS_INCLUDE_DIR})

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
        )

message(STATUS "Test programs for MathUtils module...")

foreach(test ${TESTS})
    message(STATUS "...add ${test}")

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} MathUtils)

endforeach()
