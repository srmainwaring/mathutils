
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


# Version informations
set(MathUtils_VERSION_MAJOR   "0")
set(MathUtils_VERSION_MINOR   "1")
set(MathUtils_VERSION_PATCH   "0")
set(MathUtils_VERSION
        "${MathUtils_VERSION_MAJOR}.${MathUtils_VERSION_MINOR}.${MathUtils_VERSION_PATCH}")

# MathUtils depends on Eigen3 library
include(Add_Eigen3)

# Creating MathUtils target
add_subdirectory(src/MathUtils)


# Adding tests
message(STATUS "Building MathUtils tests : ${MATHUTILS_BUILD_TESTS}")
if (${MATHUTILS_BUILD_TESTS})
    add_subdirectory(tests)
endif()
