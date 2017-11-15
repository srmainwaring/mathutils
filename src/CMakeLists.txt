
message(STATUS "=========== MathUtils module =========")

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")


# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

#message(${EIGEN3_INCLUDE_DIR})

# FILES

set(FILES
        Matrix.h
        Matrix.cpp

        )

add_library(MathUtils SHARED ${FILES})
target_link_libraries(MathUtils Eigen3::Eigen)
