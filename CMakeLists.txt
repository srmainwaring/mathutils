
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(MathUtils LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Cmake scripts
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# Languages requirements
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

option(MATHUTILS_BUILD_TESTS "Build the MathUtils tests" ON)
option(MATHUTILS_MKL "Build MathUtils with MKL support" ON)

# Version informations
set(MathUtils_VERSION_MAJOR   "0")
set(MathUtils_VERSION_MINOR   "1")
set(MathUtils_VERSION_PATCH   "0")
set(MathUtils_VERSION
        "${MathUtils_VERSION_MAJOR}.${MathUtils_VERSION_MINOR}.${MathUtils_VERSION_PATCH}")

# MathUtils depends on Eigen3 library
include(Add_Eigen3)

# Adding MKL support
if (${MATHUTILS_MKL})
    message(STATUS ---- Looking for MKL library ----)
    find_package(MKL)

    message(STATUS "   MKL include dirs:   ${MKL_INCLUDE_DIRS}")
    message(STATUS "   MKL libraries:      ${MKL_LIBRARIES}")
    message(STATUS "   IOMP5 library:      ${IOMP5_LIBRARY}")
#    message(STATUS "   MATH library:       ${MATH_LIBRARY}")
    message(STATUS "   MKL library dirs:   ${MKL_LIBRARY_DIRS}")

endif()

# Creating MathUtils target
add_subdirectory(src/MathUtils)


# Adding tests
message(STATUS "Building MathUtils tests : ${MATHUTILS_BUILD_TESTS}")
if (${MATHUTILS_BUILD_TESTS})
    add_subdirectory(tests)
endif()
