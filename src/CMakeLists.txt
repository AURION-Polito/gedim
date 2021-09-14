# INCLUDE SUB DIRECTORIES
add_subdirectory(Common)
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(MpiTools)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
