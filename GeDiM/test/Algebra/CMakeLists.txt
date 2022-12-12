list(APPEND gedim_test_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_LAPACK_utilities.hpp)
list(APPEND gedim_test_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Eigen.hpp)
list(APPEND gedim_test_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_Pardiso.hpp)

list(APPEND gedim_test_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_test_src_sources ${gedim_test_src_sources} PARENT_SCOPE)
set(gedim_test_src_headers ${gedim_test_src_headers} PARENT_SCOPE)
set(gedim_test_src_includes ${gedim_test_src_includes} PARENT_SCOPE)
