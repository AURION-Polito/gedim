# Install script for directory: D:/Nuova cartella/Intersection/gedim

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "D:/Nuova cartella/Intersection/gedim/Debug/library/")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/static" TYPE STATIC_LIBRARY FILES "D:/Nuova cartella/Intersection/gedim/Debug/libGedim.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "D:/Nuova cartella/Intersection/gedim/src/Common/CommonUtilities.hpp"
    "D:/Nuova cartella/Intersection/gedim/Debug/additional_include/Macro.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/IO/IOUtilities.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/IO/FileTextReader.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/IO/IFileTextReader.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/IMeshDAO.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshDAOExporterToCsv.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshDAOImporterFromCsv.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshMatrices.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshMatricesDAO.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/IntersectorMesh2DSegment.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/UnionMeshSegment.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshUtilities.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/ConformerMeshSegment.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/ConformerMeshPolygon.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Mesh/MeshFromCsvUtilities.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/MpiTools/IMpiProcess.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/MpiTools/MpiParallelEnvironment.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/MpiTools/MpiProcess.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Geometry/GeometryUtilities.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Geometry/MapTriangle.hpp"
    "D:/Nuova cartella/Intersection/gedim/src/Geometry/MapQuadrilateral.hpp"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/Nuova cartella/Intersection/gedim/Debug/external/cmake_install.cmake")
  include("D:/Nuova cartella/Intersection/gedim/Debug/src/cmake_install.cmake")
  include("D:/Nuova cartella/Intersection/gedim/Debug/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "D:/Nuova cartella/Intersection/gedim/Debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
