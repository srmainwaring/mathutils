find_package(Eigen3 QUIET)
if (NOT Eigen3_FOUND)
    include(FetchContent)

    set(Eigen3_URL http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz)
    FetchContent_Declare(Eigen3
            URL ${Eigen3_URL}
            )

    FetchContent_GetProperties(Eigen3)
    if(NOT Eigen3_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'Eigen3' dependency")
        FetchContent_Populate(Eigen3)

        # Eigen3 BUILD OPTIONS
        # NONE

        message("Add_Eigen3 from MathUtils")
        message("Eigen3_SOURCE_DIR = ${eigen3_SOURCE_DIR}")
        message("Eigen3_BINARY_DIR = ${eigen3_BINARY_DIR}")

        add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})
        
    endif()
endif()

if (TARGET Eigen3::Eigen)
    message("FOUND_EIGEN target in MAHTUTILS")
else()
    message("Eigen target not found in mathutils")
endif()
