#=============================================================================
#
# CMake configuration file for the MathUtils module
#
#=============================================================================

cmake_minimum_required(VERSION 3.1)
project(MathUtils)

set(CMAKE_BUILD_TYPE Debug)  # TODO: voir pour que le reglage choisi dans le projet appelant soit utilise ici...

set(CMAKE_CXX_STANDARD 11)

# Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# VERSION INFORMATIONS
# TODO: aller chercher dans un config.h !
set(MATHUTILS_MAJOR_VERSION 0)
set(MATHUTILS_MINOR_VERSION 1)
set(MATHUTILS_PATCH_VERSION 0)
set(MATHUTILS_VERSION
        ${MATHUTILS_MAJOR_VERSION}.${MATHUTILS_MINOR_VERSION}.${MATHUTILS_PATCH_VERSION})

#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Adding thirdparty projects
add_subdirectory(thirdparty)

# Adding both src and tests
#add_subdirectory(src)

# General include dirs for MathUtils (inner and thirdparty)
list(APPEND MATHUTILS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS})
include_directories(${MATHUTILS_INCLUDE_DIRS})

# Adding tests targets
add_subdirectory(tests)


# Generating configuration files to make MathUtils easily usable in other projects
include (CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MathUtilsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MathUtilsConfig.cmake
    PATH_VARS MATHUTILS_INCLUDE_DIRS MATHUTILS_THIRDPARTY_LIBS # Here we expose these two variables as replacement into the MathUtilsConfig.cmake.in
    INSTALL_DESTINATION ${CMAKE_BINARY_DIR}
#    NO_CHECK_REQUIRED_COMPONENTS_MACRO # Eigen does not provide components
)

# Generating the version file
write_basic_package_version_file(
        MathUtilsConfigVersion.cmake
        VERSION ${MATHUTILS_VERSION}
        COMPATIBILITY SameMajorVersion
)
