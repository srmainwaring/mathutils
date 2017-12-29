
set(MATHUTILS_THIRDPARTY_LIBS "")
set(MATHUTILS_THIRDPARTY_INCLUDE_DIRS "")


#=============================================================================
#
# Python 2.7 library
#
#=============================================================================
message(STATUS "===================================================================")
message(STATUS "Finding package Python...")
message(STATUS "===================================================================")
find_package(PythonLibs 2.7) # TODO: voir a terme le passage en Python 3...
list(APPEND MATHUTILS_THIRDPARTY_LIBS ${PYTHON_LIBRARIES})
list(APPEND MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})

message(STATUS "INCLUDE DIRS : " ${PYTHON_INCLUDE_DIRS})
message(STATUS "LIBRARIES : " ${PYTHON_LIBRARIES})

#=============================================================================
#
# Matplotlib-cpp
#
#=============================================================================
message(STATUS "===================================================================")
message(STATUS "Finding package Matplotlib-cpp...")
message(STATUS "===================================================================")
set(MATPLOTLIB_CPP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/matplotlib-cpp)
list(APPEND MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${MATPLOTLIB_CPP_INCLUDE_DIR})

message(STATUS "INCLUDE DIRS : " ${MATPLOTLIB_CPP_INCLUDE_DIR})
message(STATUS "LIBRARIES : NO LIB REQUIRED")


#=============================================================================
#
# Making principal variables available at parent scope (do not change !)
#
#=============================================================================

set(MATHUTILS_THIRDPARTY_LIBS ${MATHUTILS_THIRDPARTY_LIBS} PARENT_SCOPE)
set(MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS} PARENT_SCOPE)

