
set(MATHUTILS_THIRDPARTY_LIBS "")
set(MATHUTILS_THIRDPARTY_INCLUDE_DIRS "")

#=============================================================================
#
# Eigen
#
#=============================================================================

# Eigen
#message(STATUS "===================================================================")
#message(STATUS "Finding package Eigen...")
#message(STATUS "===================================================================")

#find_package(Eigen3 REQUIRED NO_MODULE)
##include_directories(${EIGEN3_INCLUDE_DIR})
#
##message(STATUS "Eigen include directory found at : " ${EIGEN3_INCLUDE_DIR})
#
#message(STATUS "INCLUDE DIRS : " ${EIGEN3_INCLUDE_DIR})
#message(STATUS "LIBRARIES : NONE" )

#=============================================================================
#
# FFTS
#
#=============================================================================

message(STATUS "===================================================================")
message(STATUS "Finding package FFTS...")
message(STATUS "===================================================================")

set(FFTS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ffts/include)
list(APPEND MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${FFTS_INCLUDE_DIR})

set(FFTS_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ffts/build/libffts.so)
list(APPEND MATHUTILS_THIRDPARTY_LIBS ${FFTS_LIB_DIR})

message(STATUS "INCLUDE DIRS : " ${FFTS_INCLUDE_DIR})
message(STATUS "LIBRARIES : " ${FFTS_LIB_DIR})

#=============================================================================
#
# Making principal variables available at parent scope (do not change !)
#
#=============================================================================

set(MATHUTILS_THIRDPARTY_LIBS ${MATHUTILS_THIRDPARTY_LIBS} PARENT_SCOPE)
set(MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS} PARENT_SCOPE)
