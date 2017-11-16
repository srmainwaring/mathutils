
message(STATUS "=========== MathUtils module =========")

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")


# Eigen
message(STATUS "Looking for Eigen...")
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

message(STATUS "Eigen include directory found at : " ${EIGEN3_INCLUDE_DIR})

# FILES

set(FILES
        Version.h
        Version.cpp
        Constants.h
        Angles.h
        Unit.h
        Interp1d.h
        LookupTable1D.h
        VectorGeneration.h
        Check.h
        Matrix.h
        )

add_library(MathUtils SHARED ${FILES})
target_link_libraries(MathUtils Eigen3::Eigen)
