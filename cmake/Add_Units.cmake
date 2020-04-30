include(FetchContent)

FetchContent_Declare(units
  GIT_REPOSITORY https://github.com/nholthaus/units.git
  GIT_TAG v2.3.1
)
FetchContent_GetProperties(units)
if(NOT units_POPULATED)
  set(BUILD_TESTS OFF CACHE BOOL "")
  set(BUILD_DOCS OFF CACHE BOOL "")
  FetchContent_Populate(units)
  add_subdirectory(${units_SOURCE_DIR} ${units_BINARY_DIR})
endif()

