
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(MathUtils LANGAUGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Cmake scripts
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# Languages requirements
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)


# Version informations
set(MathUtils_VERSION_MAJOR   "0")
set(MathUtils_VERSION_MINOR   "1")
set(MathUtils_VERSION_PATCH   "0")
set(MathUtils_VERSION
        "${MathUtils_VERSION_MAJOR}.${MathUtils_VERSION_MINOR}.${MathUtils_VERSION_PATCH}")

include(Add_Eigen3)


























##set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#
## Adding thirdparty projects
#add_subdirectory(thirdparty)
#
## Adding both src and tests
##add_subdirectory(src)
#
## General include dirs for MathUtils (inner and thirdparty)
#list(APPEND MATHUTILS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS})
#include_directories(${MATHUTILS_INCLUDE_DIRS})
#
## Adding tests targets
#add_subdirectory(tests)
#
#
## Generating configuration files to make MathUtils easily usable in other projects
#include (CMakePackageConfigHelpers)
#configure_package_config_file(
#    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MathUtilsConfig.cmake.in
#    ${CMAKE_CURRENT_BINARY_DIR}/MathUtilsConfig.cmake
#    PATH_VARS MATHUTILS_INCLUDE_DIRS MATHUTILS_THIRDPARTY_LIBS # Here we expose these two variables as replacement into the MathUtilsConfig.cmake.in
#    INSTALL_DESTINATION ${CMAKE_BINARY_DIR}
##    NO_CHECK_REQUIRED_COMPONENTS_MACRO # Eigen does not provide components
#)
#
## Generating the version file
#write_basic_package_version_file(
#        MathUtilsConfigVersion.cmake
#        VERSION ${MATHUTILS_VERSION}
#        COMPATIBILITY SameMajorVersion
#)
