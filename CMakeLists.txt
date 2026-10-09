cmake_minimum_required(VERSION 3.12)
project(GeDiMProject LANGUAGES CXX C VERSION 3.0)

set(ENABLE_GEDIM_UNIT_TESTS OFF CACHE BOOL "Enable Unit Tests")

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GeDiM-core)

# Include unit tests
if(${ENABLE_GEDIM_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GeDiM-core/test)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppcheck.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/clang-format.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/prepend_license.cmake)
