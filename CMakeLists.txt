#=============================================================================
#
# CMake configuration file for the MathUtils module
#
#=============================================================================

cmake_minimum_required(VERSION 2.8)

project(MathUtils)

set(CMAKE_BUILD_TYPE Debug)  # TODO: voir pour que le reglage choisi dans le projet appelant soit utilise ici...

set(CMAKE_CXX_STANDARD 14)

# Adding both src and tests
add_subdirectory(src)
add_subdirectory(tests)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # TODO: voir si utile...

# Defining different folders for usage in other programs
set(MATHUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
