
# matplotlib-cpp depends on Python lib
find_package(PythonLibs 2.7)
#include_directories(${PYTHON_INCLUDE_DIRS})
#message(${PYTHON_LIBRARIES})

find_package(matplotlib-cpp QUIET)
if (NOT matplotlib-cpp_FOUND)
    include(FetchContent)

    set(matplotlib-cpp_URL https://github.com/lava/matplotlib-cpp.git)
    FetchContent_Declare(matplotlib-cpp
            GIT_REPOSITORY ${matplotlib-cpp_URL}
            GIT_TAG master
            )

    FetchContent_GetProperties(matplotlib-cpp)
    if(NOT matplotlib-cpp_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'matplotlib-cpp' dependency")
        FetchContent_Populate(matplotlib-cpp)

        # matplotlib-cpp BUILD OPTIONS
        # NONE

        message(matplotlib-cpp_SOURCE_DIR${matplotlib-cpp_SOURCE_DIR})
        message(matplotlib-cpp_BINARY_DIR${matplotlib-cpp_BINARY_DIR})

#        add_subdirectory(${matplotlib-cpp_SOURCE_DIR} ${matplotlib-cpp_BINARY_DIR})
#        add_subdirectory(${matplotlib-cpp_SOURCE_DIR}/contrib)
    endif()
endif()


add_library(matplotlib-cpp INTERFACE)
target_include_directories(matplotlib-cpp INTERFACE
        $<BUILD_INTERFACE:${matplotlib-cpp_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(matplotlib-cpp INTERFACE ${PYTHON_LIBRARIES})