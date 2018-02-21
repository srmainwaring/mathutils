#=============================================================================
#
# CMake configuration file for the MathUtils module
#
#=============================================================================

cmake_minimum_required(VERSION 2.8)

project(MathUtils)

set(CMAKE_BUILD_TYPE Debug)  # TODO: voir pour que le reglage choisi dans le projet appelant soit utilise ici...

set(CMAKE_CXX_STANDARD 11)

#find_package(Eigen3 REQUIRED NO_MODULE)
#include_directories(${EIGEN3_INCLUDE_DIR})

# Adding thirdparty projects
add_subdirectory(thirdparty)



# Adding both src and tests
add_subdirectory(src)

#message(STATUS ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS})  # TODO: voir si utile...
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # TODO: voir si utile...

add_subdirectory(tests)

message(${MATHUTILS_THIRDPARTY_INCLUDE_DIRS})




## Defining VARIABLES for usage in other programs
#set(MATHUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}) # TODO: dans FRyDoM, ce sont ces variables qui doivent etre utilisees !!

