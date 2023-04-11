list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MpiParallelEnvironment.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MpiProcess.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MetisUtilities.cpp)

list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IMpiProcess.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MpiParallelEnvironment.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MpiProcess.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MetisUtilities.hpp)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
