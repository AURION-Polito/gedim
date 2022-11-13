# INCLUDE SUB DIRECTORIES
add_subdirectory(Algebra)
add_subdirectory(Common)
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(MpiTools)
add_subdirectory(Geometry)
add_subdirectory(Quadrature)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
