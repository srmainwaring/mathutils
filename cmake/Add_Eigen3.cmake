
#message(STATUS ${eigen3_POPULATED})
#message(STATUS ${eigen3_SOURCE_DIR})
#message(STATUS ${eigen3_BINARY_DIR})
set(Eigen3_URL http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz)

if (POLICY CMP0024) # Hack to make possible to include Eigen3Targets.cmake from the Eigen3 build dir
    cmake_policy(SET CMP0024 OLD)
endif()


find_package(Eigen3 QUIET)
#message(STATUS "Eigen found by find_package : ${Eigen3_FOUND}")

if (NOT Eigen3_FOUND)

#    message("EIGEN_NOT_FOUND")

    include(FetchContent)


    FetchContent_Declare(Eigen3
            URL ${Eigen3_URL}
            )

    FetchContent_GetProperties(Eigen3)
#    message(STATUS EIGEN3)
#    message(STATUS ${eigen3_POPULATED})
#    message(STATUS ${eigen3_SOURCE_DIR})
#    message(STATUS ${eigen3_BINARY_DIR})

    if(NOT eigen3_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'Eigen3' dependency")
        FetchContent_Populate(Eigen3)

        # Eigen3 BUILD OPTIONS
        # NONE
#        message("Add_Eigen3 from project")
#        message("eigen3_SOURCE_DIR = ${eigen3_SOURCE_DIR}")
#        message("eigen3_BINARY_DIR = ${eigen3_BINARY_DIR}")

        add_subdirectory(${eigen3_SOURCE_DIR} ${eigen3_BINARY_DIR})

        find_package(Eigen3) # To make the Eigen3::Eigen target readily available at first Add

#        add_library(Eigen INTERFACE)
#        target_include_directories(Eigen INTERFACE ${eigen3_SOURCE_DIR})
    else()
        message(STATUS "Eigen3 already populated")
    endif()
endif()


if (TARGET Eigen3::Eigen)
    message(STATUS "Eigen target FOUND")
    get_target_property(INC Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS ${INC})

else()
    message(STATUS "Eigen target NOT FOUND")
endif()
