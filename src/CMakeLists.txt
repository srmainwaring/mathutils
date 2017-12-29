
message(STATUS "=========== MathUtils module =========")

find_package(Eigen3 REQUIRED NO_MODULE)
#include_directories(${EIGEN3_INCLUDE_DIR})
#
#message(STATUS "Eigen include directory found at : " ${EIGEN3_INCLUDE_DIR})

# TODO: voir si utile !! A priori, MathUtils est header only...

# FILES
set(FILES
        Version.h
        Version.cpp
        Constants.h
        Angles.h
        Unit.h
        Interp.h
        Interp1d.h
        LookupTable1D.h
        Interp2d.h
        LookupTable2D.h
        Maths.h
        VectorGeneration.h
        Check.h
        Vector2d.h
        Vector3d.h
        Transform3d.h
        Matrix.h
        Integration1d.h
        FFT.h
        Complex.h
        StdVector.h
        Functions.h
        Convolution.h
        )

add_library(MathUtils SHARED ${FILES})
target_link_libraries(MathUtils Eigen3::Eigen)
