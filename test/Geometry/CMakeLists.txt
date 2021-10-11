list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_intersection.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_point.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_polygon.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_polyhedron.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry_split.hpp)

list(APPEND test_include ${CMAKE_CURRENT_SOURCE_DIR})

set(test_sources ${test_sources} PARENT_SCOPE)
set(test_headers ${test_headers} PARENT_SCOPE)
set(test_include ${test_include} PARENT_SCOPE)
