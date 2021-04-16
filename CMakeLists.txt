#=============================================================================
#
# CMake configuration file for the MathUtils project
#
#=============================================================================
cmake_minimum_required(VERSION 3.11)

project(MathUtils LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

#=============================================================================
# Options
#=============================================================================
option(MATHUTILS_BUILD_TESTS "Build the MathUtils tests" ON)
option(MATHUTILS_BUILD_BOOST_TESTS "Build the MathUtils tests requiring Boost" ON)

#=============================================================================
# Retrieving the current Git revision
#=============================================================================
set(PRE_CONFIGURE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/MathUtils/version.h.in CACHE FILEPATH "" FORCE)
set(POST_CONFIGURE_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/MathUtils/version.h CACHE FILEPATH "" FORCE)
set(GIT_WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE FILEPATH "" FORCE)
set(GIT_STATE_FILE ${CMAKE_CURRENT_BINARY_DIR}/git-state-hash CACHE FILEPATH "" FORCE)
include(cmake/git_watcher.cmake)

if (NOT EXISTS ${POST_CONFIGURE_FILE})
    configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE})
endif()

# =============================================================================
# MathUtils DEPENDENCIES
# =============================================================================
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Add_Eigen3)
include(Add_Units)

#=============================================================================
# MathUtils library
#=============================================================================
add_subdirectory(src/MathUtils)

#=============================================================================
# Tests
#=============================================================================
if (${MATHUTILS_BUILD_TESTS})
    enable_testing()
    include(Add_googletest)
    add_subdirectory(tests)
endif()
