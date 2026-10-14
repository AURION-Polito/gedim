cmake_minimum_required(VERSION 3.12)
project(GeDiMProject LANGUAGES CXX C VERSION 3.0)

set(ENABLE_GEDIM_UNIT_TESTS OFF CACHE BOOL "Enable Unit Tests")
set(ENABLE_GEDIM_EXAMPLES OFF CACHE BOOL "Enable Examples")

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GeDiM-core)

# Include 3rd_party_libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/3rd_party_libraries)

# Include unit tests
if(${ENABLE_GEDIM_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif()

if(${ENABLE_GEDIM_EXAMPLES})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Examples)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppcheck.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/clang-format.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/prepend_license.cmake)
