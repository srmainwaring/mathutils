
find_package(units QUIET)
if (NOT units_FOUND)
    include(FetchContent)

    set(units_URL https://github.com/nholthaus/units.git)
    FetchContent_Declare(units
            GIT_REPOSITORY ${units_URL}
            GIT_TAG v2.3.1  # For the moment (november 2018, the v3.xx does not build)
            )

    FetchContent_GetProperties(units)
    if(NOT units_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'units' dependency")
        FetchContent_Populate(units)

        # units BUILD OPTIONS
#        set(BUILD_TESTS OFF)

#        message(STATUS ${units_SOURCE_DIR})
#        message(STATUS ${units_BINARY_DIR})

        add_subdirectory(${units_SOURCE_DIR} ${units_BINARY_DIR})

    else()
        message(STATUS "units already populated")
    endif()
endif()


#add_library(units INTERFACE)
#target_include_directories(units INTERFACE
#        $<BUILD_INTERFACE:${units_SOURCE_DIR}>
##        $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>
#        $<INSTALL_INTERFACE:include>)
#target_link_libraries(units INTERFACE ${PYTHON_LIBRARIES})
#
##~ add_executable(generator generator.c)
#export(TARGETS units FILE units-exports.cmake)


