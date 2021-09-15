list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/IOUtilities.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/FileTextReader.cpp)

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IOUtilities.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/FileTextReader.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IFileTextReader.hpp)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
