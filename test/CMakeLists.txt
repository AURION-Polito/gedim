# INCLUDE SUB DIRECTORIES
add_subdirectory(Mesh)

list(APPEND test_include ${CMAKE_CURRENT_SOURCE_DIR})

set(test_sources ${test_sources} PARENT_SCOPE)
set(test_headers ${test_headers} PARENT_SCOPE)
set(test_include ${test_include} PARENT_SCOPE)
