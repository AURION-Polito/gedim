list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOExporterToCsv.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOImporterFromCsv.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOImporter2DFromCsv.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatricesDAO.cpp)

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/IMeshDAO.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOExporterToCsv.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOImporterFromCsv.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshDAOImporter2DFromCsv.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatrices.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/MeshMatricesDAO.hpp)

list(APPEND src_include ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_include ${src_include} PARENT_SCOPE)
