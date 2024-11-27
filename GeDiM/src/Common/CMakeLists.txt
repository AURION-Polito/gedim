list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/CommonUtilities.cpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/CommonUtilities.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Gedim_Macro.in)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
