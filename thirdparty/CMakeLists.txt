
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
message(STATUS "LIBRARIES : NO LIB REQUIRED (HEADER ONLY)")

#=============================================================================
#
# EIGEN
#
#=============================================================================
message(STATUS "===================================================================")
message(STATUS "Finding package Eigen...")
message(STATUS "===================================================================")
#set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
set(Eigen3_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen/build)

find_package(Eigen3 REQUIRED NO_MODULE)

# TESTING
if (TARGET Eigen3::Eigen)

else ()
    message(STATUS "EIGEN Not found")
endif (TARGET Eigen3::Eigen)



list(APPEND MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
list(APPEND MATHUTILS_THIRDPARTY_LIBS Eigen3::Eigen)

message(STATUS "INCLUDE DIRS : " ${EIGEN3_INCLUDE_DIR})
message(STATUS "LIBRARIES : NO LIB REQUIRED (HEADERS ONLY)")

#=============================================================================
#
# CERES SOLVER
#
#=============================================================================
message(STATUS "===================================================================")
message(STATUS "Finding Ceres-solver module...")
message(STATUS "===================================================================")
#set(Ceres_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ceres-solver/build)

set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen/eigen)
find_package(Ceres)

# FIXME: hack qui devrait etre traite autrement !!! (et pas certain que ca fonctionne...)
# TODO: suivant ce qui a ete regle dans CMAKE_BUILD_TYPE, il faudra declarer libceres.so ou libceres-debug.so
set(CERES_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/ceres-solver/build/lib/libceres.so)

list(APPEND MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${CERES_INCLUDE_DIRS})
list(APPEND MATHUTILS_THIRDPARTY_LIBS ${CERES_LIBRARIES})

message(STATUS "INCLUDE DIRS : " ${CERES_INCLUDE_DIRS})
message(STATUS "LIBRARIES : ${CERES_LIBRARIES}")


#=============================================================================
#
# Making principal variables available at parent scope (do not change !)
#
#=============================================================================

set(MATHUTILS_THIRDPARTY_LIBS ${MATHUTILS_THIRDPARTY_LIBS} PARENT_SCOPE)
set(MATHUTILS_THIRDPARTY_INCLUDE_DIRS ${MATHUTILS_THIRDPARTY_INCLUDE_DIRS} PARENT_SCOPE)

