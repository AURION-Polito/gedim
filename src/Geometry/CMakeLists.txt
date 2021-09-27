list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/GeometryUtilities.cpp)

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GeometryUtilities.hpp)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
