list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_intersectorMesh2DSegment.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_unionMeshSegment.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_conformMesh.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_export_conformMesh.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatrices_2D_2Cells_Mock.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatrices_2D_4Cells_Mock.hpp)
list(APPEND test_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatrices_2D_26Cells_Mock.hpp)

list(APPEND test_include ${CMAKE_CURRENT_SOURCE_DIR})

set(test_sources ${test_sources} PARENT_SCOPE)
set(test_headers ${test_headers} PARENT_SCOPE)
set(test_include ${test_include} PARENT_SCOPE)
