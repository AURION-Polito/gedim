list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MpiParallelEnvironment.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MpiProcess.cpp)

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IMpiProcess.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MpiParallelEnvironment.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MpiProcess.hpp)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
