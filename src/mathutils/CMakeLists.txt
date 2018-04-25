
#message(STATUS "=========== MathUtils module =========")

# FIXME: this file is no more used. It coumd be used when we will work on install instructions


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
        VectorN.h
        )

#add_library(MathUtils ${FILES} IMPORTED)
#target_link_libraries(MathUtils Eigen3::Eigen)
#target_link_libraries(MathUtils)


#install(TARGETS MathUtils
#        EXPORT MathUtilsTargets
#        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)

#message(STATUS ${CMAKE_CURRENT_LIST_DIR})
#set(MATHUTILS_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} PARENT_SCOPE)