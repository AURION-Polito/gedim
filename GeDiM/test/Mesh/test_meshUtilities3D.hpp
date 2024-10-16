#ifndef __TEST_MESH_UTILITIES3D_H
#define __TEST_MESH_UTILITIES3D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Macro.hpp"

#include "GraphUtilities.hpp"
#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatrices_3D_1Cells_Mock.hpp"
#include "MeshMatrices_3D_68Cells_Mock.hpp"
#include "OVM_Mesh_Mock.hpp"

#include "MeshFromCsvUtilities.hpp"
#include "MeshDAOImporterFromCsv.hpp"
#include "OpenVolumeMeshInterface.hpp"
#include "PlatonicSolid.hpp"
#include "SphereMeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "test_meshUtilities2D.hpp"
#include "VTK_Unstructured_Grid_Mesh_Mock.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestMeshUtilities, TestMesh3DFromPolyhedron)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    const vector<unsigned int> vertexMarkers = { 1, 2, 3, 4, 5, 6, 7, 8 };
    const vector<unsigned int> edgeMarkers = { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
    const vector<unsigned int> faceMarkers = { 21, 22, 23, 24, 25, 26 };

    meshUtilities.Mesh3DFromPolyhedron(cube.Vertices,
                                       cube.Edges,
                                       cube.Faces,
                                       vertexMarkers,
                                       edgeMarkers,
                                       faceMarkers,
                                       meshDao);

    std::string exportFolder = "./Export/TestMesh3DFromPolyhedron/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 3);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 8);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 12);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 6);
  }

  TEST(TestMeshUtilities, TestCreateDelaunayMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;


    Eigen::MatrixXd points(3,27);

    points.col(0)<< 0.0, 0.5, 0.0;
    points.col(1)<< 0.5, 0.0, 0.0;
    points.col(2)<< 0.5, 0.5, 0.0;
    points.col(3)<< 0.0, 0.0, 0.5;
    points.col(4)<< 0.5, 0.0, 0.5;
    points.col(5)<< 0.0, 0.5, 0.5;
    points.col(6)<< 0.5, 0.5, 0.5;
    points.col(7)<< 0.5, 1.0, 0.0;
    points.col(8)<< 0.5, 1.0, 0.5;
    points.col(9)<< 0.0, 1.0, 0.5;
    points.col(10)<< 1.0, 0.0, 0.5;
    points.col(11)<< 1.0, 0.5, 0.0;
    points.col(12)<< 1.0, 0.5, 0.5;
    points.col(13)<< 1.0, 1.0, 0.5;
    points.col(14)<< 0.0, 0.5, 1.0;
    points.col(15)<< 0.5, 0.0, 1.0;
    points.col(16)<< 0.5, 0.5, 1.0;
    points.col(17)<< 1.0, 0.5, 1.0;
    points.col(18)<< 0.5, 1.0, 1.0;
    points.col(19)<< 0.0, 0.0, 0.0;
    points.col(20)<< 1.0, 0.0, 0.0;
    points.col(21)<< 1.0, 1.0, 0.0;
    points.col(22)<< 0.0, 1.0, 0.0;
    points.col(23)<< 0.0, 0.0, 1.0;
    points.col(24)<< 1.0, 0.0, 1.0;
    points.col(25)<< 1.0, 1.0, 1.0;
    points.col(26)<< 0.0, 1.0, 1.0;

    std::vector<unsigned int> points_marker(27, 0);
    points_marker[19] = 1;
    points_marker[20] = 2;
    points_marker[21] = 3;
    points_marker[22] = 4;
    points_marker[23] = 5;
    points_marker[24] = 6;
    points_marker[25] = 7;
    points_marker[26] = 8;


    meshUtilities.CreateDelaunayMesh3D(points,
                                       points_marker,
                                       meshDao);

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreateDelaunayMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    std::vector<std::vector<bool>> cell3Ds_faces_orientation(meshDao.Cell3DTotalNumber(),
                                                             std::vector<bool>(4, true));


    meshUtilities.ExportMeshToOpenVolume(meshDao,
                                         cell3Ds_faces_orientation,
                                         exportFolder + "/mesh.ovm");

    EXPECT_EQ(3,
              meshDao.Dimension());
    EXPECT_EQ(27,
              meshDao.Cell0DTotalNumber());
    EXPECT_EQ(98,
              meshDao.Cell1DTotalNumber());
    EXPECT_EQ(120,
              meshDao.Cell2DTotalNumber());
    EXPECT_EQ(48,
              meshDao.Cell3DTotalNumber());
  }

  TEST(TestMeshUtilities, TestCreateTetrahedralMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                   1.0);

    meshUtilities.CreateTetrahedralMesh(polyhedron.Vertices,
                                        polyhedron.Edges,
                                        polyhedron.Faces,
                                        0.03,
                                        meshDao,
                                        "Qpqfezna");

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreateTetrahedralMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "CreatedTetrahedralMesh");

    EXPECT_EQ(3,
              meshDao.Dimension());
    EXPECT_EQ(28,
              meshDao.Cell0DTotalNumber());
    EXPECT_EQ(103,
              meshDao.Cell1DTotalNumber());
    EXPECT_EQ(127,
              meshDao.Cell2DTotalNumber());
    EXPECT_EQ(51,
              meshDao.Cell3DTotalNumber());
  }

  TEST(TestMeshUtilities, TestCreatePolyhedralMesh)
  {
#if ENABLE_VORO == 0
    GTEST_SKIP_("Voro module not activated.");
#endif

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreatePolyhedralMesh";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                   1.0);

    const unsigned int meshGenerator = 0;
    const unsigned int numCells1D = 2;
    switch (meshGenerator)
    {
      case 0:
      {
        const Eigen::Vector3d parallelepipedsLengthTangent = polyhedron.Vertices.col(1) -
                                                             polyhedron.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsHeightTangent = polyhedron.Vertices.col(4) -
                                                             polyhedron.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsWidthTangent = polyhedron.Vertices.col(3) -
                                                            polyhedron.Vertices.col(0);

        const std::vector<double> lengthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numCells1D + 1,
                                                                                                            0.0, 1.0, 1);
        const std::vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numCells1D + 1,
                                                                                                            0.0, 1.0, 1);
        const std::vector<double> widthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(numCells1D + 1,
                                                                                                           0.0, 1.0, 1);

        const Eigen::Vector3d origin = polyhedron.Vertices.col(0);
        meshUtilities.CreateParallelepipedMesh(origin,
                                               parallelepipedsLengthTangent,
                                               parallelepipedsHeightTangent,
                                               parallelepipedsWidthTangent,
                                               lengthMeshCurvilinearCoordinates,
                                               heightMeshCurvilinearCoordinates,
                                               widthMeshCurvilinearCoordinates,
                                               meshDao);
      }
        break;
      case 1:
      {
        meshUtilities.CreatePolyhedralMesh(geometryUtilities,
                                           polyhedron.Vertices,
                                           polyhedron.Edges,
                                           polyhedron.Faces,
                                           numCells1D,
                                           10,
                                           meshDao);
      }
        break;
      default:
        ASSERT_FALSE(true);
    }

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "TestCreatePolyhedralMesh");

    {
      ExportMeshData exportData;
      exportData.Cell0Ds = meshDao.Cell0DsCoordinates();
      exportData.Cell1Ds = meshDao.Cell1DsExtremes();
      exportData.Cell2Ds = meshDao.Cell2DsExtremes();
      exportData.Cell3DsVertices = meshDao.Cell3DsVertices();
      exportData.Cell3DsEdges = meshDao.Cell3DsEdges();
      exportData.Cell3DsFaces = meshDao.Cell3DsFaces();
      exportData.Cell0DsMarker = meshDao.Cell0DsMarker();
      exportData.Cell1DsMarker = meshDao.Cell1DsMarker();
      exportData.Cell2DsMarker = meshDao.Cell2DsMarker();
      exportData.Cell3DsMarker = meshDao.Cell3DsMarker();

      ExportMeshUtilities::ExportMesh3DToText(exportData,
                                              exportFolder + "/MeshData.txt");
    }

    Gedim::MeshUtilities::MeshGeometricData3D cell3DsGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                           meshDao);

    {
      ExportMeshGeometricData3D exportData;
      exportData.PolyhedronsVertices = cell3DsGeometricData.Cell3DsVertices;
      exportData.PolyhedronsEdges = cell3DsGeometricData.Cell3DsEdges;
      exportData.PolyhedronsFaces = cell3DsGeometricData.Cell3DsFaces;
      exportData.PolyhedronsVolume = cell3DsGeometricData.Cell3DsVolumes;
      exportData.PolyhedronsDiameter = cell3DsGeometricData.Cell3DsDiameters;
      exportData.PolyhedronsCentroid = cell3DsGeometricData.Cell3DsCentroids;
      exportData.PolyhedronsTetrahedronsVertices = cell3DsGeometricData.Cell3DsTetrahedronPoints;
      exportData.PolyhedronsFacesTranslation = cell3DsGeometricData.Cell3DsFacesTranslations;
      exportData.PolyhedronsFacesRotationMatrix = cell3DsGeometricData.Cell3DsFacesRotationMatrices;
      exportData.PolyhedronsFacesNormalDirection = cell3DsGeometricData.Cell3DsFacesNormalDirections;
      exportData.PolyhedronsFacesEdgesDirection = cell3DsGeometricData.Cell3DsFacesEdgeDirections;
      exportData.PolyhedronsFaces2DVertices = cell3DsGeometricData.Cell3DsFaces2DVertices;
      exportData.PolyhedronsFacesTriangulations2DVertices = cell3DsGeometricData.Cell3DsFaces2DTriangulations;
      exportData.PolyhedronsFaces2DCentroid = cell3DsGeometricData.Cell3DsFaces2DCentroids;
      exportData.PolyhedronsFacesEdgesLength = cell3DsGeometricData.Cell3DsFacesEdgeLengths;
      exportData.PolyhedronsFacesEdges2DCentroid = cell3DsGeometricData.Cell3DsFacesEdges2DCentroid;
      exportData.PolyhedronsFacesEdges2DTangent = cell3DsGeometricData.Cell3DsFacesEdge2DTangents;
      exportData.PolyhedronsFacesEdges2DTangentNormalized = cell3DsGeometricData.Cell3DsFacesEdge2DTangents;
      exportData.PolyhedronsFacesEdges2DNormal = cell3DsGeometricData.Cell3DsFacesEdge2DNormals;

      exportData.PolyhedronsFacesArea.resize(cell3DsGeometricData.Cell3DsFacesAreas.size());
      for (unsigned int c = 0; c < cell3DsGeometricData.Cell3DsFacesAreas.size(); c++)
      {
        exportData.PolyhedronsFacesArea[c].resize(cell3DsGeometricData.Cell3DsFacesAreas[c].size());
        for (unsigned int f = 0; f < cell3DsGeometricData.Cell3DsFacesAreas[c].size(); f++)
        {
          exportData.PolyhedronsFacesArea[c][f] =
              cell3DsGeometricData.Cell3DsFacesAreas[c][f];
        }
      }

      exportData.PolyhedronsFacesNormal.resize(cell3DsGeometricData.Cell3DsFacesNormals.size());
      for (unsigned int c = 0; c < cell3DsGeometricData.Cell3DsFacesNormals.size(); c++)
      {
        exportData.PolyhedronsFacesNormal[c].resize(3, cell3DsGeometricData.Cell3DsFacesNormals[c].size());
        for (unsigned int f = 0; f < cell3DsGeometricData.Cell3DsFacesNormals[c].size(); f++)
        {
          exportData.PolyhedronsFacesNormal[c].col(f) =
              cell3DsGeometricData.Cell3DsFacesNormals[c][f];
        }
      }

      exportData.PolyhedronsFacesDiameter.resize(cell3DsGeometricData.Cell3DsFacesDiameters.size());
      for (unsigned int c = 0; c < cell3DsGeometricData.Cell3DsFacesDiameters.size(); c++)
      {
        exportData.PolyhedronsFacesDiameter[c].resize(cell3DsGeometricData.Cell3DsFacesDiameters[c].size());
        for (unsigned int f = 0; f < cell3DsGeometricData.Cell3DsFacesDiameters[c].size(); f++)
        {
          exportData.PolyhedronsFacesDiameter[c][f] =
              cell3DsGeometricData.Cell3DsFacesDiameters[c][f];
        }
      }

      for (unsigned int c = 0; c < exportData.PolyhedronsFacesEdges2DTangentNormalized.size(); c++)
      {
        for (unsigned int f = 0; f < exportData.PolyhedronsFacesEdges2DTangentNormalized[c].size(); f++)
        {
          for (unsigned int e = 0; e < exportData.PolyhedronsFacesEdges2DTangentNormalized[c][f].cols(); e++)
            exportData.PolyhedronsFacesEdges2DTangentNormalized[c][f].col(e).normalize();
        }
      }

      ExportMeshUtilities::ExportMeshGeometricData3DToText(exportData,
                                                           exportFolder + "/MeshGeometricData.txt");
    }
  }

  TEST(TestMeshUtilities, TestCheckMesh3D)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_68Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshUtilities::CheckMesh3DConfiguration config;
    ASSERT_NO_THROW(meshUtilities.CheckMesh3D(config,
                                              geometryUtilities,
                                              meshDao));
  }

  TEST(TestMeshUtilities, TestComputeCell2DCell3DNeighbours)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    mesh.Mesh.NumberCell2DNeighbourCell3D.clear();
    mesh.Mesh.Cell2DNeighbourCell3Ds.clear();
    mesh.Mesh.NumberCell2DNeighbourCell3D.resize(mesh.Mesh.NumberCell2D + 1, 0);
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDao);

    EXPECT_EQ(mesh.Mesh.NumberCell2DNeighbourCell3D,
              vector<unsigned int>({ 0,1,2,3,4,5, 6 }));
    EXPECT_EQ(mesh.Mesh.Cell2DNeighbourCell3Ds,
              vector<unsigned int>({ 0,0,0,0,0,0 }));
  }

  TEST(TestMeshUtilities, TestFillMesh3DGeometricData_Convex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData3D result = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                   meshDao);

    Gedim::MeshUtilities::MeshGeometricData3D expectedResult;
    expectedResult.Cell3DsVolumes = { 1.0 };
    expectedResult.Cell3DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.5) };
    expectedResult.Cell3DsDiameters = { sqrt(3.0) };
    expectedResult.Cell3DsVertices = { cube.Vertices };
    expectedResult.Cell3DsEdges = { cube.Edges };
    expectedResult.Cell3DsFaces = { cube.Faces };

    ASSERT_EQ(result.Cell3DsVertices, expectedResult.Cell3DsVertices);
    ASSERT_EQ(result.Cell3DsEdges, expectedResult.Cell3DsEdges);
    ASSERT_EQ(result.Cell3DsFaces, expectedResult.Cell3DsFaces);
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsVolumes[0], expectedResult.Cell3DsVolumes[0], geometryUtilities.Tolerance3D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].x(), expectedResult.Cell3DsCentroids[0].z(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].y(), expectedResult.Cell3DsCentroids[0].y(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].z(), expectedResult.Cell3DsCentroids[0].x(), geometryUtilities.Tolerance1D()));
    ASSERT_EQ(result.Cell3DsDiameters, expectedResult.Cell3DsDiameters);

    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig,
                                                           geometryUtilities,
                                                           meshDao,
                                                           result));
  }

  TEST(TestMeshUtilities, TestFillMesh3DGeometricData_Concave)
  {
    std::string exportFolder = "./Export/TestFillMesh3DGeometricData_Concave";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConcaveMesh");

    const Gedim::MeshUtilities::MeshGeometricData3D result = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                   meshDao,
                                                                                                   meshDao,
                                                                                                   std::vector<std::vector<unsigned int>> { { 0 } });

    {
      for (unsigned int c = 0; c < meshDao.Cell3DTotalNumber(); c++)
      {
        const unsigned int numFaces = result.Cell3DsFaces[c].size();

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < result.Cell3DsFaces3DTriangulations[c].size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            for (unsigned int t = 0; t < result.Cell3DsFaces3DTriangulations[c][f].size(); t++)
            {
              exporter.AddPolygon(result.Cell3DsFaces3DTriangulations[c][f][t],
                                  {
                                    {
                                      "Face",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(faceIndex.size()),
                                      faceIndex.data()
                                    }
                                  });
            }
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_Faces2DTriangulations.vtu");
        }

        std::vector<Eigen::Vector3d> faceInternalPoints(numFaces);
        for (unsigned int f = 0; f < numFaces; f++)
          faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(result.Cell3DsFaces3DTriangulations[c][f][0]);

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            exporter.AddPoint(faceInternalPoints[f],
                              {
                                {
                                  "Face",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(faceIndex.size()),
                                  faceIndex.data()
                                }
                              });
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_FacesInternalPoint.vtu");
        }

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            exporter.AddSegment(faceInternalPoints[f],
                                faceInternalPoints[f] + result.Cell3DsFacesNormals[c][f],
                                {
                                  {
                                    "Face",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(faceIndex.size()),
                                    faceIndex.data()
                                  }
                                });
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_FacesNormal.vtu");
        }
      }
    }

    Gedim::MeshUtilities::MeshGeometricData3D expectedResult;
    expectedResult.Cell3DsVolumes = { 1.0 };
    expectedResult.Cell3DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.5) };
    expectedResult.Cell3DsDiameters = { sqrt(3.0) };
    expectedResult.Cell3DsVertices = { cube.Vertices };
    expectedResult.Cell3DsEdges = { cube.Edges };
    expectedResult.Cell3DsFaces = { cube.Faces };
    expectedResult.Cell3DsFacesNormalDirections = { { false, true, false, true, true, false } };

    ASSERT_EQ(result.Cell3DsVertices, expectedResult.Cell3DsVertices);
    ASSERT_EQ(result.Cell3DsEdges, expectedResult.Cell3DsEdges);
    ASSERT_EQ(result.Cell3DsFaces, expectedResult.Cell3DsFaces);
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsVolumes[0], expectedResult.Cell3DsVolumes[0], geometryUtilities.Tolerance3D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].x(), expectedResult.Cell3DsCentroids[0].z(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].y(), expectedResult.Cell3DsCentroids[0].y(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].z(), expectedResult.Cell3DsCentroids[0].x(), geometryUtilities.Tolerance1D()));
    ASSERT_EQ(result.Cell3DsDiameters, expectedResult.Cell3DsDiameters);
    ASSERT_EQ(result.Cell3DsFacesNormalDirections,
              expectedResult.Cell3DsFacesNormalDirections);

    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig,
                                                           geometryUtilities,
                                                           meshDao,
                                                           result));
  }

  TEST(TestMeshUtilities, TestFillMesh3DGeometricData_Concave_Tetra)
  {
    std::string exportFolder = "./Export/TestFillMesh3DGeometricData_Concave_Tetra";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ConcaveMesh");

    const auto mesh_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                           meshDao);
    std::vector<std::vector<Eigen::MatrixXd>> cell3Ds_tetra_vertices(meshDao.Cell3DTotalNumber());
    std::vector<std::vector<Eigen::Matrix3d>> cell2Ds_triangles_vertices(meshDao.Cell2DTotalNumber());

    for (unsigned int c = 0; c < meshDao.Cell3DTotalNumber(); c++)
    {
      if (!meshDao.Cell3DIsActive(c))
        continue;

      cell3Ds_tetra_vertices[c] = mesh_geometric_data.Cell3DsTetrahedronPoints.at(c);

      for (unsigned int f = 0; f < meshDao.Cell3DNumberFaces(c); f++)
      {
        const unsigned int c2D_index = meshDao.Cell3DFace(c, f);
        cell2Ds_triangles_vertices[c2D_index] = mesh_geometric_data.Cell3DsFaces3DTriangulations.at(c).at(f);
      }
    }


    const Gedim::MeshUtilities::MeshGeometricData3D result = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                   meshDao,
                                                                                                   cell3Ds_tetra_vertices,
                                                                                                   cell2Ds_triangles_vertices);

    {
      for (unsigned int c = 0; c < meshDao.Cell3DTotalNumber(); c++)
      {
        const unsigned int numFaces = result.Cell3DsFaces[c].size();

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < result.Cell3DsFaces3DTriangulations[c].size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            for (unsigned int t = 0; t < result.Cell3DsFaces3DTriangulations[c][f].size(); t++)
            {
              exporter.AddPolygon(result.Cell3DsFaces3DTriangulations[c][f][t],
                                  {
                                    {
                                      "Face",
                                      Gedim::VTPProperty::Formats::Cells,
                                      static_cast<unsigned int>(faceIndex.size()),
                                      faceIndex.data()
                                    }
                                  });
            }
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_Faces2DTriangulations.vtu");
        }

        std::vector<Eigen::Vector3d> faceInternalPoints(numFaces);
        for (unsigned int f = 0; f < numFaces; f++)
          faceInternalPoints[f] = geometryUtilities.PolygonBarycenter(result.Cell3DsFaces3DTriangulations[c][f][0]);

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            exporter.AddPoint(faceInternalPoints[f],
                              {
                                {
                                  "Face",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(faceIndex.size()),
                                  faceIndex.data()
                                }
                              });
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_FacesInternalPoint.vtu");
        }

        {
          Gedim::VTKUtilities exporter;

          for (unsigned int f = 0; f < faceInternalPoints.size(); f++)
          {
            std::vector<double> faceIndex(1, f);
            exporter.AddSegment(faceInternalPoints[f],
                                faceInternalPoints[f] + result.Cell3DsFacesNormals[c][f],
                                {
                                  {
                                    "Face",
                                    Gedim::VTPProperty::Formats::Cells,
                                    static_cast<unsigned int>(faceIndex.size()),
                                    faceIndex.data()
                                  }
                                });
          }

          exporter.Export(exportFolder + "/Cell3D_" +
                          to_string(c) + "_FacesNormal.vtu");
        }
      }
    }

    Gedim::MeshUtilities::MeshGeometricData3D expectedResult;
    expectedResult.Cell3DsVolumes = { 1.0 };
    expectedResult.Cell3DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.5) };
    expectedResult.Cell3DsDiameters = { sqrt(3.0) };
    expectedResult.Cell3DsVertices = { cube.Vertices };
    expectedResult.Cell3DsEdges = { cube.Edges };
    expectedResult.Cell3DsFaces = { cube.Faces };
    expectedResult.Cell3DsFacesNormalDirections = { { false, true, false, true, true, false } };

    ASSERT_EQ(result.Cell3DsVertices, expectedResult.Cell3DsVertices);
    ASSERT_EQ(result.Cell3DsEdges, expectedResult.Cell3DsEdges);
    ASSERT_EQ(result.Cell3DsFaces, expectedResult.Cell3DsFaces);
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsVolumes[0], expectedResult.Cell3DsVolumes[0], geometryUtilities.Tolerance3D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].x(), expectedResult.Cell3DsCentroids[0].z(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].y(), expectedResult.Cell3DsCentroids[0].y(), geometryUtilities.Tolerance1D()));
    ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.Cell3DsCentroids[0].z(), expectedResult.Cell3DsCentroids[0].x(), geometryUtilities.Tolerance1D()));
    ASSERT_EQ(result.Cell3DsDiameters, expectedResult.Cell3DsDiameters);
    ASSERT_EQ(result.Cell3DsFacesNormalDirections,
              expectedResult.Cell3DsFacesNormalDirections);

    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig,
                                                           geometryUtilities,
                                                           meshDao,
                                                           result));
  }

  TEST(TestMeshUtilities, TestSetMeshMarkersOnPlane)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestSetMeshMarkersOnPlane";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Original");

    meshUtilities.SetMeshMarkersOnPlane(geometryUtilities,
                                        Eigen::Vector3d { -1.0, 0.0, 0.0 },
                                        Eigen::Vector3d { 0.0, 0.0, 0.0 },
                                        100,
                                        meshDao);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Modified");

    ASSERT_EQ(100, meshDao.Cell0DMarker(0));
    ASSERT_EQ(100, meshDao.Cell0DMarker(3));
    ASSERT_EQ(100, meshDao.Cell0DMarker(4));
    ASSERT_EQ(100, meshDao.Cell0DMarker(7));
    ASSERT_EQ(100, meshDao.Cell0DMarker(11));
    ASSERT_EQ(100, meshDao.Cell0DMarker(15));
    ASSERT_EQ(100, meshDao.Cell0DMarker(16));
    ASSERT_EQ(100, meshDao.Cell0DMarker(17));

    ASSERT_EQ(100, meshDao.Cell1DMarker(27));
    ASSERT_EQ(100, meshDao.Cell1DMarker(28));
    ASSERT_EQ(100, meshDao.Cell1DMarker(32));
    ASSERT_EQ(100, meshDao.Cell1DMarker(33));
    ASSERT_EQ(100, meshDao.Cell1DMarker(34));
    ASSERT_EQ(100, meshDao.Cell1DMarker(39));
    ASSERT_EQ(100, meshDao.Cell1DMarker(45));
    ASSERT_EQ(100, meshDao.Cell1DMarker(46));
    ASSERT_EQ(100, meshDao.Cell1DMarker(47));
    ASSERT_EQ(100, meshDao.Cell1DMarker(48));
    ASSERT_EQ(100, meshDao.Cell1DMarker(50));
    ASSERT_EQ(100, meshDao.Cell1DMarker(52));
    ASSERT_EQ(100, meshDao.Cell1DMarker(53));

    ASSERT_EQ(100, meshDao.Cell2DMarker(25));
    ASSERT_EQ(100, meshDao.Cell2DMarker(28));
    ASSERT_EQ(100, meshDao.Cell2DMarker(42));
    ASSERT_EQ(100, meshDao.Cell2DMarker(47));
    ASSERT_EQ(100, meshDao.Cell2DMarker(49));
    ASSERT_EQ(100, meshDao.Cell2DMarker(52));
  }

  TEST(TestMeshUtilities, TestSetMeshMarkersByFaceNormal)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestSetMeshMarkersByFaceNormal";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Original");

    meshUtilities.ComputeCell2DCell3DNeighbours(meshDao);

    std::vector<Eigen::Vector3d> cell2Ds_normal(meshDao.Cell2DTotalNumber());
    for (unsigned int f = 0; f < meshDao.Cell2DTotalNumber(); f++)
    {
      cell2Ds_normal[f] =
          geometryUtilities.PolygonNormal(meshDao.Cell2DVerticesCoordinates(f));
    }

    meshUtilities.SetMeshMarkersByFaceNormal(geometryUtilities,
                                             Eigen::Vector3d { -1.0, 0.0, 0.0 },
                                             cell2Ds_normal,
                                             100,
                                             meshDao);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh_Modified");

    ASSERT_EQ(100, meshDao.Cell0DMarker(0));
    ASSERT_EQ(100, meshDao.Cell0DMarker(3));
    ASSERT_EQ(100, meshDao.Cell0DMarker(4));
    ASSERT_EQ(100, meshDao.Cell0DMarker(7));
    ASSERT_EQ(100, meshDao.Cell0DMarker(11));
    ASSERT_EQ(100, meshDao.Cell0DMarker(15));
    ASSERT_EQ(100, meshDao.Cell0DMarker(16));
    ASSERT_EQ(100, meshDao.Cell0DMarker(17));
    ASSERT_EQ(100, meshDao.Cell0DMarker(5));
    ASSERT_EQ(100, meshDao.Cell0DMarker(19));
    ASSERT_EQ(100, meshDao.Cell0DMarker(1));
    ASSERT_EQ(100, meshDao.Cell0DMarker(9));
    ASSERT_EQ(100, meshDao.Cell0DMarker(13));
    ASSERT_EQ(100, meshDao.Cell0DMarker(6));
    ASSERT_EQ(100, meshDao.Cell0DMarker(18));
    ASSERT_EQ(100, meshDao.Cell0DMarker(2));

    ASSERT_EQ(100, meshDao.Cell1DMarker(27));
    ASSERT_EQ(100, meshDao.Cell1DMarker(28));
    ASSERT_EQ(100, meshDao.Cell1DMarker(32));
    ASSERT_EQ(100, meshDao.Cell1DMarker(33));
    ASSERT_EQ(100, meshDao.Cell1DMarker(34));
    ASSERT_EQ(100, meshDao.Cell1DMarker(39));
    ASSERT_EQ(100, meshDao.Cell1DMarker(45));
    ASSERT_EQ(100, meshDao.Cell1DMarker(46));
    ASSERT_EQ(100, meshDao.Cell1DMarker(47));
    ASSERT_EQ(100, meshDao.Cell1DMarker(48));
    ASSERT_EQ(100, meshDao.Cell1DMarker(50));
    ASSERT_EQ(100, meshDao.Cell1DMarker(52));
    ASSERT_EQ(100, meshDao.Cell1DMarker(53));
    ASSERT_EQ(100, meshDao.Cell1DMarker(38));
    ASSERT_EQ(100, meshDao.Cell1DMarker(36));
    ASSERT_EQ(100, meshDao.Cell1DMarker(10));
    ASSERT_EQ(100, meshDao.Cell1DMarker(22));
    ASSERT_EQ(100, meshDao.Cell1DMarker(19));
    ASSERT_EQ(100, meshDao.Cell1DMarker(21));
    ASSERT_EQ(100, meshDao.Cell1DMarker(16));
    ASSERT_EQ(100, meshDao.Cell1DMarker(41));
    ASSERT_EQ(100, meshDao.Cell1DMarker(17));
    ASSERT_EQ(100, meshDao.Cell1DMarker(40));
    ASSERT_EQ(100, meshDao.Cell1DMarker(15));
    ASSERT_EQ(100, meshDao.Cell1DMarker(55));
    ASSERT_EQ(100, meshDao.Cell1DMarker(57));

    ASSERT_EQ(100, meshDao.Cell2DMarker(25));
    ASSERT_EQ(100, meshDao.Cell2DMarker(28));
    ASSERT_EQ(100, meshDao.Cell2DMarker(42));
    ASSERT_EQ(100, meshDao.Cell2DMarker(47));
    ASSERT_EQ(100, meshDao.Cell2DMarker(49));
    ASSERT_EQ(100, meshDao.Cell2DMarker(52));
    ASSERT_EQ(100, meshDao.Cell2DMarker(33));
    ASSERT_EQ(100, meshDao.Cell2DMarker(19));
    ASSERT_EQ(100, meshDao.Cell2DMarker(15));
    ASSERT_EQ(100, meshDao.Cell2DMarker(11));
    ASSERT_EQ(100, meshDao.Cell2DMarker(40));
    ASSERT_EQ(100, meshDao.Cell2DMarker(59));
  }

  TEST(TestMeshUtilities, TestImportOVMMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    std::vector<std::vector<bool>> meshCell3DsFacesOrientation;
    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::OpenVolumeMeshInterface ovmInterface;

    const std::vector<string> lines = OVM_Mesh_Mock::FileLines();
    const Gedim::OpenVolumeMeshInterface::OVMMesh ovm_mesh = ovmInterface.StringsToOVMMesh(lines);
    ovmInterface.OVMMeshToMeshDAO(ovm_mesh,
                                  meshDao,
                                  meshCell3DsFacesOrientation);

    std::string exportFolder = "./Export/TestMeshUtilities/TestImportOVMMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ImportedOVMMesh");

    ASSERT_EQ(3,
              meshDao.Dimension());
    ASSERT_EQ(64,
              meshDao.Cell0DTotalNumber());
    ASSERT_EQ(144,
              meshDao.Cell1DTotalNumber());
    ASSERT_EQ(108,
              meshDao.Cell2DTotalNumber());
    ASSERT_EQ(27,
              meshDao.Cell3DTotalNumber());

    const Gedim::OpenVolumeMeshInterface::OVMMesh reconverted_ovm_mesh = ovmInterface.MeshDAOToOVMMesh(meshDao,
                                                                                                       meshCell3DsFacesOrientation);

    ASSERT_EQ(ovm_mesh.NumCell0Ds, reconverted_ovm_mesh.NumCell0Ds);
    ASSERT_EQ(ovm_mesh.NumCell1Ds, reconverted_ovm_mesh.NumCell1Ds);
    ASSERT_EQ(ovm_mesh.NumCell2Ds, reconverted_ovm_mesh.NumCell2Ds);
    ASSERT_EQ(ovm_mesh.NumCell3Ds, reconverted_ovm_mesh.NumCell3Ds);
    ASSERT_EQ(ovm_mesh.Cell0Ds, reconverted_ovm_mesh.Cell0Ds);
    ASSERT_EQ(ovm_mesh.Cell1Ds, reconverted_ovm_mesh.Cell1Ds);
    ASSERT_EQ(ovm_mesh.Cell2Ds, reconverted_ovm_mesh.Cell2Ds);
    for (unsigned int p = 0; p < ovm_mesh.NumCell3Ds; p++)
    {
      ASSERT_EQ(ovm_mesh.Cell3Ds[p].FacesIndex, reconverted_ovm_mesh.Cell3Ds[p].FacesIndex);
      ASSERT_EQ(ovm_mesh.Cell3Ds[p].FacesOrientation, reconverted_ovm_mesh.Cell3Ds[p].FacesOrientation);
    }

    const std::vector<std::string> reconverted_lines = ovmInterface.OVMMeshToStrings(reconverted_ovm_mesh);

    ASSERT_EQ(lines.size(), reconverted_lines.size());
    for (unsigned int l = 0; l < lines.size(); l++)
      ASSERT_EQ(lines[l], reconverted_lines[l]);

    meshUtilities.ExportMeshToOpenVolume(meshDao,
                                         meshCell3DsFacesOrientation,
                                         exportFolder + "/hex_parallel_1.ovm");
  }

  TEST(TestMeshUtilities, TestRefineMesh3D)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestRefineMesh3D/";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_R0");

    // first refine
    {
      const unsigned int newCell1DIndex = meshDAO.Cell1DAppend(2);

      meshDAO.Cell1DInsertExtremes(newCell1DIndex, 7, 5);
      meshDAO.Cell1DSetMarker(newCell1DIndex, meshDAO.Cell2DMarker(1));
      meshDAO.Cell1DSetState(newCell1DIndex, true);

      meshDAO.Cell1DInsertExtremes(newCell1DIndex + 1, 3, 1);
      meshDAO.Cell1DSetMarker(newCell1DIndex + 1, meshDAO.Cell2DMarker(0));
      meshDAO.Cell1DSetState(newCell1DIndex + 1, true);

      EXPECT_EQ(12, newCell1DIndex);

      Eigen::MatrixXi subCell_0_0(2, 3);
      subCell_0_0.row(0)<< 0, 3, 1;
      subCell_0_0.row(1)<< 3, newCell1DIndex + 1, 0;
      Eigen::MatrixXi subCell_0_1(2, 3);
      subCell_0_1.row(0)<< 1, 3, 2;
      subCell_0_1.row(1)<< newCell1DIndex + 1, 2, 1;

      const std::vector<unsigned int> newCell2DsIndex_0 = meshUtilities.SplitCell2D(0,
                                                                                    { subCell_0_0, subCell_0_1 },
                                                                                    meshDAO);
      EXPECT_EQ(std::vector<unsigned int>({ 6, 7 }), newCell2DsIndex_0);
      EXPECT_FALSE(meshDAO.Cell2DIsActive(0));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(6));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(7));
      EXPECT_EQ(meshDAO.Cell2DMarker(0), meshDAO.Cell2DMarker(6));
      EXPECT_EQ(meshDAO.Cell2DMarker(0), meshDAO.Cell2DMarker(7));
      EXPECT_EQ(0, meshDAO.Cell2DOriginalCell2D(6));
      EXPECT_EQ(0, meshDAO.Cell2DOriginalCell2D(7));
      std::list<unsigned int> updatedCell2Ds_0;
      EXPECT_FALSE(meshDAO.Cell2DUpdatedCell2Ds(0,
                                                updatedCell2Ds_0));
      updatedCell2Ds_0.sort();
      EXPECT_EQ(std::list<unsigned int>({6, 7}), updatedCell2Ds_0);

      Eigen::MatrixXi subCell_1_0(2, 3);
      subCell_1_0.row(0)<< 7, 4, 5;
      subCell_1_0.row(1)<< 7, 4, newCell1DIndex;
      Eigen::MatrixXi subCell_1_1(2, 3);
      subCell_1_1.row(0)<< 7, 5, 6;
      subCell_1_1.row(1)<< newCell1DIndex, 5, 6;

      const std::vector<unsigned int> newCell2DsIndex_1 = meshUtilities.SplitCell2D(1,
                                                                                    { subCell_1_0, subCell_1_1 },
                                                                                    meshDAO);
      EXPECT_EQ(std::vector<unsigned int>({ 8, 9 }), newCell2DsIndex_1);
      EXPECT_FALSE(meshDAO.Cell2DIsActive(1));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(8));
      EXPECT_TRUE(meshDAO.Cell2DIsActive(9));
      EXPECT_EQ(meshDAO.Cell2DMarker(1), meshDAO.Cell2DMarker(8));
      EXPECT_EQ(meshDAO.Cell2DMarker(1), meshDAO.Cell2DMarker(9));
      EXPECT_EQ(1, meshDAO.Cell2DOriginalCell2D(9));
      EXPECT_EQ(1, meshDAO.Cell2DOriginalCell2D(9));
      std::list<unsigned int> updatedCell2Ds_1;
      EXPECT_FALSE(meshDAO.Cell2DUpdatedCell2Ds(1,
                                                updatedCell2Ds_1));
      updatedCell2Ds_1.sort();
      EXPECT_EQ(std::list<unsigned int>({8, 9}), updatedCell2Ds_1);

      Eigen::MatrixXi newCell2D(2, 4);
      newCell2D.row(0)<< 7, 5, 1, 3;
      newCell2D.row(1)<< newCell1DIndex, 9, newCell1DIndex + 1, 11;
      const unsigned int newCell2DIndex = meshDAO.Cell2DAppend(1);

      meshDAO.Cell2DAddVerticesAndEdges(newCell2DIndex,
                                        newCell2D);
      meshDAO.Cell2DSetMarker(newCell2DIndex, 0);
      meshDAO.Cell2DSetState(newCell2DIndex, true);
      meshDAO.Cell2DInitializeNeighbourCell3Ds(newCell2DIndex, 2);

      const std::vector<std::vector<unsigned int>> subCell3DsVertices =
      {
        { 0,1,3,7,4,5 },
        { 1,2,3,5,6,7 }
      };
      const std::vector<std::vector<unsigned int>> subCell3DsEdges =
      {
        { 0,3,newCell1DIndex + 1,8,9,11,4,7,newCell2DIndex },
        { 1,2,newCell1DIndex + 1,10,9,11,5,6,newCell2DIndex }
      };
      const std::vector<std::vector<unsigned int>> subCell3DsFaces =
      {
        { 4,2,newCell2DsIndex_0[0], newCell2DsIndex_1[0] },
        { 3,5,newCell2DsIndex_0[1], newCell2DsIndex_1[1] }
      };

      const std::vector<unsigned int> newCell3DsIndex = meshUtilities.SplitCell3D(0,
                                                                                  subCell3DsVertices,
                                                                                  subCell3DsEdges,
                                                                                  subCell3DsFaces,
                                                                                  meshDAO);
      EXPECT_EQ(std::vector<unsigned int>({ 1, 2 }), newCell3DsIndex);
      EXPECT_FALSE(meshDAO.Cell3DIsActive(0));
      EXPECT_TRUE(meshDAO.Cell3DIsActive(1));
      EXPECT_TRUE(meshDAO.Cell3DIsActive(2));
      EXPECT_EQ(meshDAO.Cell3DMarker(0), meshDAO.Cell3DMarker(1));
      EXPECT_EQ(meshDAO.Cell3DMarker(0), meshDAO.Cell3DMarker(2));
      EXPECT_EQ(0, meshDAO.Cell3DOriginalCell3D(1));
      EXPECT_EQ(0, meshDAO.Cell3DOriginalCell3D(2));
      std::list<unsigned int> updatedCell3Ds;
      EXPECT_FALSE(meshDAO.Cell3DUpdatedCell3Ds(0,
                                                updatedCell3Ds));
      updatedCell3Ds.sort();
      EXPECT_EQ(std::list<unsigned int>({1, 2}), updatedCell3Ds);

      meshDAO.Cell2DInsertNeighbourCell3D(newCell2DIndex,
                                          0,
                                          newCell3DsIndex[1]); // right
      meshDAO.Cell2DInsertNeighbourCell3D(newCell2DIndex,
                                          1,
                                          newCell3DsIndex[0]); // left

      EXPECT_EQ(1, meshDAO.Cell2DNeighbourCell3D(2, 0));
      EXPECT_EQ(1, meshDAO.Cell2DNeighbourCell3D(4, 0));
      EXPECT_EQ(1, meshDAO.Cell2DNeighbourCell3D(6, 0));
      EXPECT_EQ(1, meshDAO.Cell2DNeighbourCell3D(8, 0));
      EXPECT_EQ(1, meshDAO.Cell2DNeighbourCell3D(10, 1));
      EXPECT_EQ(2, meshDAO.Cell2DNeighbourCell3D(3, 0));
      EXPECT_EQ(2, meshDAO.Cell2DNeighbourCell3D(5, 0));
      EXPECT_EQ(2, meshDAO.Cell2DNeighbourCell3D(7, 0));
      EXPECT_EQ(2, meshDAO.Cell2DNeighbourCell3D(9, 0));
      EXPECT_EQ(2, meshDAO.Cell2DNeighbourCell3D(10, 0));

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R1");
    }


    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Final");

  }

  TEST(TestMeshUtilities, TestFillMesh3D)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestFillMesh3D";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    Eigen::MatrixXd vertices = (Eigen::MatrixXd(3, 9)<<
                                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0).finished();
    vertices.col(4)<< 0.5 * (vertices.col(1) + vertices.col(3));
    vertices.col(5)<< 0.5 * (vertices.col(0) + vertices.col(2));
    vertices.col(6)<< 0.5 * (vertices.col(4) + vertices.col(2));
    vertices.col(7)<< 0.5 * (vertices.col(5) + vertices.col(1));

    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 21)<<
                                   0, 0, 1, 2, 0, 2, 4, 4, 4, 0, 2, 2, 6, 6, 3, 4, 5, 7, 0, 1, 2,
                                   3, 1, 2, 3, 4, 6, 6, 7, 5, 5, 5, 7, 5, 7, 4, 1, 7, 1, 8, 8, 8).finished();

    vector<Eigen::MatrixXi> triangles(19);
    triangles[0] = (Eigen::MatrixXi(2, 3)<<
                    0,1,4,
                    1,15,4).finished();
    triangles[1] = (Eigen::MatrixXi(2, 3)<<
                    0,4,3,
                    4,14,0).finished();
    triangles[2] = (Eigen::MatrixXi(2, 4)<<
                    0,5,2,3,
                    9,10,3,0).finished();
    triangles[3] = (Eigen::MatrixXi(2, 3)<<
                    0,8,1,
                    18,19,1).finished();
    triangles[4] = (Eigen::MatrixXi(2, 4)<<
                    0,8,2,5,
                    18,20,10,9).finished();
    triangles[5] = (Eigen::MatrixXi(2, 4)<<
                    0,1,7,5,
                    1,17,16,9).finished();
    triangles[6] = (Eigen::MatrixXi(2, 3)<<
                    7,1,2,
                    17,2,11).finished();
    triangles[7] = (Eigen::MatrixXi(2, 3)<<
                    5,7,2,
                    16,11,10).finished();
    triangles[8] = (Eigen::MatrixXi(2, 3)<<
                    2,1,8,
                    2,19,20).finished();
    triangles[9] = (Eigen::MatrixXi(2, 3)<<
                    0,4,5,
                    4,8,9).finished();
    triangles[10] = (Eigen::MatrixXi(2, 3)<<
                     5,4,6,
                     8,6,12).finished();
    triangles[11] = (Eigen::MatrixXi(2, 3)<<
                     5,6,2,
                     12,5,10).finished();
    triangles[12] = (Eigen::MatrixXi(2, 3)<<
                     6,7,5,
                     13,16,12).finished();
    triangles[13] = (Eigen::MatrixXi(2, 4)<<
                     3,4,6,2,
                     14,6,5,3).finished();
    triangles[14] = (Eigen::MatrixXi(2, 3)<<
                     4,7,6,
                     7,13,6).finished();
    triangles[15] = (Eigen::MatrixXi(2, 3)<<
                     6,7,2,
                     13,11,5).finished();

    triangles[16] = (Eigen::MatrixXi(2, 3)<<
                     4,1,7,
                     15,17,7).finished();
    triangles[17] = (Eigen::MatrixXi(2, 4)<<
                     4,1,2,6,
                     15,2,5,6).finished();
    triangles[18] = (Eigen::MatrixXi(2, 3)<<
                     4,7,5,
                     7,16,8).finished();

    const vector<Gedim::MeshUtilities::Mesh3DPolyhedron> polyhedrons =
    {
      {
        { 0,1,2,5,7,8 },
        { 18,19,20,1,2,9,10,11,16,17 },
        { 5,6,7,3,4,8 }
      },
      {
        { 0,2,3,4,5,6 },
        { 0,3,14,4,9,10,5,6,8,12 },
        { 1,2,9,10,11,13 }
      },
      {
        { 0,1,4,5,7 },
        { 1,4,15,9,16,17,7,8 },
        { 0,5,9,16,18 }
      },
      {
        { 4,5,6,7 },
        { 6,7,8,12,13,16 },
        { 10,12,14,18 }
      },
      {
        { 2,5,6,7 },
        { 5,10,11,12,13,16 },
        { 7,11,12,15 }
      },
      {
        { 1,2,4,6,7 },
        { 2,5,6,7,11,13,15,17 },
        { 17,6,14,15,16 }
      }
    };

    meshUtilities.FillMesh3D(vertices,
                             edges,
                             triangles,
                             polyhedrons,
                             meshDao);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    ASSERT_EQ(vertices,
              meshDao.Cell0DsCoordinates());
    ASSERT_EQ(edges,
              meshDao.Cell1DsExtremes());
    Gedim::MeshUtilities::CheckMesh3DConfiguration checkConfig;
    meshUtilities.CheckMesh3D(checkConfig,
                              geometryUtilities,
                              meshDao);
  }

  TEST(TestMeshUtilities, TestComputeMesh3DAlignedCell1Ds)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestComputeMesh3DAlignedCell1Ds";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;
    Gedim::GraphUtilities graphUtilities;

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Eigen::MatrixXd vertices = (Eigen::MatrixXd(3, 9)<<
                                0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0).finished();
    vertices.col(4)<< 0.5 * (vertices.col(1) + vertices.col(3));
    vertices.col(5)<< 0.5 * (vertices.col(0) + vertices.col(2));
    vertices.col(6)<< 0.5 * (vertices.col(4) + vertices.col(2));
    vertices.col(7)<< 0.5 * (vertices.col(5) + vertices.col(1));

    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 21)<<
                                   0, 0, 1, 2, 0, 2, 4, 4, 4, 0, 2, 2, 6, 6, 3, 4, 5, 7, 0, 1, 2,
                                   3, 1, 2, 3, 4, 6, 6, 7, 5, 5, 5, 7, 5, 7, 4, 1, 7, 1, 8, 8, 8).finished();

    vector<Eigen::MatrixXi> triangles(19);
    triangles[0] = (Eigen::MatrixXi(2, 3)<<
                    0,1,4,
                    1,15,4).finished();
    triangles[1] = (Eigen::MatrixXi(2, 3)<<
                    0,4,3,
                    4,14,0).finished();
    triangles[2] = (Eigen::MatrixXi(2, 4)<<
                    0,5,2,3,
                    9,10,3,0).finished();
    triangles[3] = (Eigen::MatrixXi(2, 3)<<
                    0,8,1,
                    18,19,1).finished();
    triangles[4] = (Eigen::MatrixXi(2, 4)<<
                    0,8,2,5,
                    18,20,10,9).finished();
    triangles[5] = (Eigen::MatrixXi(2, 4)<<
                    0,1,7,5,
                    1,17,16,9).finished();
    triangles[6] = (Eigen::MatrixXi(2, 3)<<
                    7,1,2,
                    17,2,11).finished();
    triangles[7] = (Eigen::MatrixXi(2, 3)<<
                    5,7,2,
                    16,11,10).finished();
    triangles[8] = (Eigen::MatrixXi(2, 3)<<
                    2,1,8,
                    2,19,20).finished();
    triangles[9] = (Eigen::MatrixXi(2, 3)<<
                    0,4,5,
                    4,8,9).finished();
    triangles[10] = (Eigen::MatrixXi(2, 3)<<
                     5,4,6,
                     8,6,12).finished();
    triangles[11] = (Eigen::MatrixXi(2, 3)<<
                     5,6,2,
                     12,5,10).finished();
    triangles[12] = (Eigen::MatrixXi(2, 3)<<
                     6,7,5,
                     13,16,12).finished();
    triangles[13] = (Eigen::MatrixXi(2, 4)<<
                     3,4,6,2,
                     14,6,5,3).finished();
    triangles[14] = (Eigen::MatrixXi(2, 3)<<
                     4,7,6,
                     7,13,6).finished();
    triangles[15] = (Eigen::MatrixXi(2, 3)<<
                     6,7,2,
                     13,11,5).finished();

    triangles[16] = (Eigen::MatrixXi(2, 3)<<
                     4,1,7,
                     15,17,7).finished();
    triangles[17] = (Eigen::MatrixXi(2, 4)<<
                     4,1,2,6,
                     15,2,5,6).finished();
    triangles[18] = (Eigen::MatrixXi(2, 3)<<
                     4,7,5,
                     7,16,8).finished();

    const vector<Gedim::MeshUtilities::Mesh3DPolyhedron> polyhedrons =
    {
      {
        { 0,1,2,5,7,8 },
        { 18,19,20,1,2,9,10,11,16,17 },
        { 5,6,7,3,4,8 }
      },
      {
        { 0,2,3,4,5,6 },
        { 0,3,14,4,9,10,5,6,8,12 },
        { 1,2,9,10,11,13 }
      },
      {
        { 0,1,4,5,7 },
        { 1,4,15,9,16,17,7,8 },
        { 0,5,9,16,18 }
      },
      {
        { 4,5,6,7 },
        { 6,7,8,12,13,16 },
        { 10,12,14,18 }
      },
      {
        { 2,5,6,7 },
        { 5,10,11,12,13,16 },
        { 7,11,12,15 }
      },
      {
        { 1,2,4,6,7 },
        { 2,5,6,7,11,13,15,17 },
        { 17,6,14,15,16 }
      }
    };

    const std::vector<unsigned int> cell2DsAligned =
    { 0,1,2,3,4,5,5,5,8,9,9,9,12,13,14,14,16,17,18 };

    meshUtilities.FillMesh3D(vertices,
                             edges,
                             triangles,
                             polyhedrons,
                             mesh);

    meshUtilities.ExportMeshToVTU(mesh,
                                  exportFolder,
                                  "Mesh");

    const Gedim::MeshUtilities::MeshGeometricData3D meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                              mesh);


    std::vector<std::vector<std::vector<unsigned int>>> cell3DsAlignedEdgesVertices(mesh.Cell3DTotalNumber());
    std::vector<std::vector<std::vector<unsigned int>>> cell3DsAlignedEdgesEdges(mesh.Cell3DTotalNumber());
    std::vector<Eigen::VectorXd> cell3DsAlignedEdgesLenght(mesh.Cell3DTotalNumber());

    for (unsigned int c = 0; c < mesh.Cell3DTotalNumber(); c++)
    {
      if (!mesh.Cell3DIsActive(c))
        continue;

      const Eigen::MatrixXd& cell3DVertices = meshGeometricData.Cell3DsVertices.at(c);
      const Eigen::MatrixXi& cell3DEdges = meshGeometricData.Cell3DsEdges.at(c);
      const std::vector<Eigen::MatrixXi>& cell3DFaces = meshGeometricData.Cell3DsFaces.at(c);

      const Gedim::GraphUtilities::GraphAdjacencyData cell3DAdjacency = graphUtilities.GraphConnectivityToGraphAdjacency(cell3DVertices.cols(),
                                                                                                                         cell3DEdges,
                                                                                                                         false);

      const Gedim::GeometryUtilities::AlignedPolyhedronEdgesResult alignedEdges = geometryUtilities.AlignedPolyhedronEdges(cell3DVertices,
                                                                                                                           cell3DAdjacency.GraphAdjacencyVertices,
                                                                                                                           cell3DAdjacency.GraphAdjacencyEdges,
                                                                                                                           cell3DAdjacency.GraphAdjacencyVerticesMap,
                                                                                                                           meshGeometricData.Cell3DsEdgeTangents.at(c),
                                                                                                                           meshGeometricData.Cell3DsEdgeLengths.at(c).array().square());

      cell3DsAlignedEdgesVertices[c] = alignedEdges.AlignedEdgesVertices;
      cell3DsAlignedEdgesEdges[c] = alignedEdges.AlignedEdgesEdges;
    }

    Gedim::MeshUtilities::ComputeMesh3DAlignedCell1DsResult result =
        meshUtilities.ComputeMesh3DAlignedCell1Ds(cell3DsAlignedEdgesVertices,
                                                  cell3DsAlignedEdgesEdges,
                                                  mesh);

    Eigen::MatrixXi alignedCell1Ds_expected(2, 24);
    alignedCell1Ds_expected.row(0)<< 0, 0, 0, 1, 2, 1, 1, 2, 0, 0, 2, 3, 2, 4, 5, 0, 1, 4, 4, 6, 5, 2, 2, 1;
    alignedCell1Ds_expected.row(1)<< 8, 1, 2, 8, 8, 2, 5, 7, 3, 4, 3, 4, 4, 5, 6, 5, 4, 7, 6, 7, 7, 6, 5, 7;

    std::vector<Eigen::MatrixXi> cell0DsAlignedCell1DsIndex_expected =
    {
      (Eigen::MatrixXi(2,6) << 0,  1,  2,  8,  9, 15,
      0,  0,  0,  0,  0,  0).finished(),
      (Eigen::MatrixXi(2,6) << 1,  3,  5,  6, 16, 23,
      1,  0,  0,  0,  0,  0).finished(),
      (Eigen::MatrixXi(2,8) << 2,  4,  5,  7, 10, 12, 21, 22,
      2,  0,  1,  0,  0,  0,  0,  0).finished(),
      (Eigen::MatrixXi(2,3) << 8, 10, 11,
      1,  1,  0).finished(),
      (Eigen::MatrixXi(2,7) << 9, 11, 12, 13, 16, 17, 18,
      1,  1,  2,  0,  1,  0,  0).finished(),
      (Eigen::MatrixXi(2,7) << 2,  6, 13, 14, 15, 20, 22,
      1,  2,  1,  0,  1,  0,  1).finished(),
      (Eigen::MatrixXi(2,5) << 12, 14, 18, 19, 21,
      1 , 1 , 1 , 0 , 1).finished(),
      (Eigen::MatrixXi(2,6) << 6,  7, 17, 19, 20, 23,
      1,  1,  1,  1,  1,  1).finished(),
      (Eigen::MatrixXi(2,3) << 0, 3, 4,
      1, 1, 1).finished()
    };

    std::vector<Eigen::MatrixXi> cell1DsAlignedCell1DsIndex_expected =
    {
      (Eigen::MatrixXi(2,1) << 8,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 1,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 5,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 10,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 9,
      0).finished(),
      (Eigen::MatrixXi(2,2) << 12, 21,
      0,  0).finished(),
      (Eigen::MatrixXi(2,2) << 12, 18,
      1,  0).finished(),
      (Eigen::MatrixXi(2,1) << 17,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 13,
      0).finished(),
      (Eigen::MatrixXi(2,2) << 2, 15,
      0,  0).finished(),
      (Eigen::MatrixXi(2,2) << 2, 22,
      1,  0).finished(),
      (Eigen::MatrixXi(2,1) << 7,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 14,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 19,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 11,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 16,
      0).finished(),
      (Eigen::MatrixXi(2,2) << 6, 20,
      1,  0).finished(),
      (Eigen::MatrixXi(2,2) << 6, 23,
      0,  0).finished(),
      (Eigen::MatrixXi(2,1) << 0,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 3,
      0).finished(),
      (Eigen::MatrixXi(2,1) << 4,
      0).finished()
    };

    std::vector<Eigen::MatrixXi> cell3DsAlignedCell1DsIndex_expected =
    {
      (Eigen::MatrixXi(2,8) << 0, 1, 2, 3, 4, 5, 6, 7,
      0, 0, 0, 0, 0, 0, 0, 0).finished(),
      (Eigen::MatrixXi(2,8) << 8,  9,  2, 10, 11, 12, 13, 14,
      0,  0,  1,  0,  0,  0,  0,  0).finished(),
      (Eigen::MatrixXi(2,7) << 1,  9, 15, 16,  6, 17, 13,
      1,  1,  0,  0,  1,  0,  1).finished(),
      (Eigen::MatrixXi(2,6) << 18, 17, 13, 14, 19, 20,
      0 , 1 , 2 , 1 , 0 , 0).finished(),
      (Eigen::MatrixXi(2,6) << 21, 22,  7, 14, 19, 20,
      0 , 0 , 1 , 2 , 1 , 1).finished(),
      (Eigen::MatrixXi(2,7) << 5, 16, 23, 12,  7, 17, 19,
      1,  1,  0,  1,  2,  2,  2).finished()
    };

    ASSERT_EQ(alignedCell1Ds_expected,
              result.AlignedCell1Ds);
    ASSERT_EQ(
          std::vector<std::vector<unsigned int>>({{0,8}, {0,1}, {0,5,2}, {1,8}, {2,8}, {1,2}, {1,7,5}, {2,7}, {0,3}, {0,4}, {2,3}, {3,4}, {2,6,4}, {4,5}, {5,6}, {0,5}, {1,4}, {4,7}, {4,6}, {6,7}, {5,7}, {2,6}, {2,5}, {1,7}}),
          result.AlignedCell1Ds_SubCell0Ds);
    ASSERT_EQ(
          std::vector<std::vector<unsigned int>>({{18}, {1}, {9,10}, {19}, {20}, {2}, {17,16}, {11}, {0}, {4}, {3}, {14}, {5,6}, {8}, {12}, {9}, {15}, {7}, {6}, {13}, {16}, {5}, {10}, {17}}),
          result.AlignedCell1Ds_SubCell1Ds);
    ASSERT_EQ(
          std::vector<std::vector<unsigned int>>({{0}, {0,2}, {0,1}, {0}, {0}, {0,5}, {0,2}, {0,4,5}, {1}, {1,2}, {1}, {1}, {1,5}, {1,2,3}, {1,3,4}, {2}, {2,5}, {2,3,5}, {3}, {3,4,5}, {3,4}, {4}, {4}, {5}}),
          result.AlignedCell1Ds_Cell3Ds);
    ASSERT_EQ(cell0DsAlignedCell1DsIndex_expected,
              result.Cell0DsAlignedCell1DsIndex);
    ASSERT_EQ(cell1DsAlignedCell1DsIndex_expected,
              result.Cell1DsAlignedCell1DsIndex);
    ASSERT_EQ(cell3DsAlignedCell1DsIndex_expected,
              result.Cell3DsAlignedCell1DsIndex);
  }

  TEST(TestMeshUtilities, TestAgglomerateCell3Ds_ByFace)
  {
    auto sort_result = [](const auto& array)
    {
      std::vector<unsigned int> array_copy(array.begin(),
                                           array.end());
      std::sort(array_copy.begin(),
                array_copy.end());

      return array_copy;
    };

    std::string exportFolder = "./Export/TestMeshUtilities/TestAgglomerateCell3Ds_ByFace";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    meshUtilities.ComputeCell0DCell3DNeighbours(meshDao);
    meshUtilities.ComputeCell0DCell1DNeighbours(meshDao);
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDao);
    meshUtilities.ComputeCell2DCell3DNeighbours(meshDao);

    std::vector<std::vector<unsigned int>> meshCell3DToConvexCell3DIndices(meshDao.Cell3DTotalNumber());
    for (unsigned int c3D_index = 0; c3D_index < meshDao.Cell3DTotalNumber(); c3D_index++)
    {
      meshCell3DToConvexCell3DIndices.at(c3D_index) =
          std::vector<unsigned int>({ c3D_index });
    }

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "OriginalMesh");

    {
      const auto cell3Ds = meshDao.Cell2DNeighbourCell3Ds(1);
      const auto agglomerationInfo = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                       cell3Ds.rend()),
                                                                      meshDao);

      const unsigned int agglomeratedCell3DIndex = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                                     cell3Ds.rend()),
                                                                                    agglomerationInfo.AgglomerateCell3DVertices,
                                                                                    agglomerationInfo.AgglomerateCell3DEdges,
                                                                                    agglomerationInfo.AgglomerateCell3DFaces,
                                                                                    agglomerationInfo.SubCell3DsRemovedVertices,
                                                                                    agglomerationInfo.SubCell3DsRemovedEdges,
                                                                                    agglomerationInfo.SubCell3DsRemovedFaces,
                                                                                    meshDao,
                                                                                    meshCell3DToConvexCell3DIndices);

      ASSERT_EQ(0,
                agglomeratedCell3DIndex);
    }

    {
      const auto cell3Ds = meshDao.Cell2DNeighbourCell3Ds(4);
      const auto agglomerationInfo = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                       cell3Ds.rend()),
                                                                      meshDao);


      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 5,19,12,13,8 })),
                sort_result(agglomerationInfo.AgglomerateCell3DVertices));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 38,7,36,10,6,8,37,11,9 })),
                sort_result(agglomerationInfo.AgglomerateCell3DEdges));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 34, 35, 6, 5, 33, 7 })),
                sort_result(agglomerationInfo.AgglomerateCell3DFaces));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 4 })),
                sort_result(agglomerationInfo.SubCell3DsRemovedFaces));

      const unsigned int agglomeratedCell3DIndex = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                                     cell3Ds.rend()),
                                                                                    agglomerationInfo.AgglomerateCell3DVertices,
                                                                                    agglomerationInfo.AgglomerateCell3DEdges,
                                                                                    agglomerationInfo.AgglomerateCell3DFaces,
                                                                                    agglomerationInfo.SubCell3DsRemovedVertices,
                                                                                    agglomerationInfo.SubCell3DsRemovedEdges,
                                                                                    agglomerationInfo.SubCell3DsRemovedFaces,
                                                                                    meshDao,
                                                                                    meshCell3DToConvexCell3DIndices);

      ASSERT_EQ(22,
                agglomeratedCell3DIndex);
      ASSERT_EQ(std::vector<unsigned int>({ 1, 10 }),
                meshCell3DToConvexCell3DIndices.at(22));
    }

    {
      const auto cell3Ds = meshDao.Cell2DNeighbourCell3Ds(7);
      const auto agglomerationInfo = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                       cell3Ds.rend()),
                                                                      meshDao);
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 14,8,13,12,19,5 })),
                sort_result(agglomerationInfo.AgglomerateCell3DVertices));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 2,35,9,11,37,8,6,23,10,36,7,38 })),
                sort_result(agglomerationInfo.AgglomerateCell3DEdges));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 38, 20, 33, 5, 32, 6, 35, 34 })),
                sort_result(agglomerationInfo.AgglomerateCell3DFaces));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 7 })),
                sort_result(agglomerationInfo.SubCell3DsRemovedFaces));

      const unsigned int agglomeratedCell3DIndex = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                                     cell3Ds.rend()),
                                                                                    agglomerationInfo.AgglomerateCell3DVertices,
                                                                                    agglomerationInfo.AgglomerateCell3DEdges,
                                                                                    agglomerationInfo.AgglomerateCell3DFaces,
                                                                                    agglomerationInfo.SubCell3DsRemovedVertices,
                                                                                    agglomerationInfo.SubCell3DsRemovedEdges,
                                                                                    agglomerationInfo.SubCell3DsRemovedFaces,
                                                                                    meshDao,
                                                                                    meshCell3DToConvexCell3DIndices);

      ASSERT_EQ(23,
                agglomeratedCell3DIndex);
      ASSERT_EQ(std::vector<unsigned int>({ 1, 10, 12 }),
                meshCell3DToConvexCell3DIndices.at(23));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData activeMeshData;
    meshUtilities.ExtractActiveMesh(meshDao,
                                    activeMeshData);
    std::vector<std::vector<unsigned int>> activeMeshCell3DToConvexCell3DIndices(meshDao.Cell3DTotalNumber());
    for (unsigned int c3D_index = 0; c3D_index < meshDao.Cell3DTotalNumber(); c3D_index++)
    {
      activeMeshCell3DToConvexCell3DIndices.at(c3D_index) =
          meshCell3DToConvexCell3DIndices.at(activeMeshData.NewCell3DToOldCell3D.at(c3D_index));
    }

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "AgglomeratedMesh");
  }

  TEST(TestMeshUtilities, TestAgglomerateCell3Ds_ByVertex)
  {
    auto sort_result = [](const auto& array)
    {
      std::vector<unsigned int> array_copy(array.begin(),
                                           array.end());
      std::sort(array_copy.begin(),
                array_copy.end());

      return array_copy;
    };

    std::string exportFolder = "./Export/TestMeshUtilities/TestAgglomerateCell3Ds_ByVertex";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    meshUtilities.ComputeCell0DCell3DNeighbours(meshDao);
    meshUtilities.ComputeCell0DCell1DNeighbours(meshDao);
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDao);
    meshUtilities.ComputeCell2DCell3DNeighbours(meshDao);

    std::vector<std::vector<unsigned int>> meshCell3DToConvexCell3DIndices(meshDao.Cell3DTotalNumber());
    for (unsigned int c3D_index = 0; c3D_index < meshDao.Cell3DTotalNumber(); c3D_index++)
    {
      meshCell3DToConvexCell3DIndices.at(c3D_index) =
          std::vector<unsigned int>({ c3D_index });
    }

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "OriginalMesh");

    {
      const auto cell3Ds = meshDao.Cell0DNeighbourCell3Ds(18);
      const auto agglomerationInfo = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.begin(),
                                                                                                       cell3Ds.end()),
                                                                      meshDao);

      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 6,14,13,2,18,10,9,8 })),
                sort_result(agglomerationInfo.AgglomerateCell3DVertices));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 5,40,41,42,16,3,55,18,57,15,2,58,0,56,17,35,8,1 })),
                sort_result(agglomerationInfo.AgglomerateCell3DEdges));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 3,39,41,11,1,40,14,32,59,57,60,61 })),
                sort_result(agglomerationInfo.AgglomerateCell3DFaces));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ })),
                sort_result(agglomerationInfo.SubCell3DsRemovedVertices));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 4 })),
                sort_result(agglomerationInfo.SubCell3DsRemovedEdges));
      ASSERT_EQ(sort_result(std::vector<unsigned int>({ 31,2,12,13,0,58 })),
                sort_result(agglomerationInfo.SubCell3DsRemovedFaces));

      const unsigned int agglomeratedCell3DIndex = meshUtilities.AgglomerateCell3Ds(std::unordered_set<unsigned int>(cell3Ds.rbegin(),
                                                                                                                     cell3Ds.rend()),
                                                                                    agglomerationInfo.AgglomerateCell3DVertices,
                                                                                    agglomerationInfo.AgglomerateCell3DEdges,
                                                                                    agglomerationInfo.AgglomerateCell3DFaces,
                                                                                    agglomerationInfo.SubCell3DsRemovedVertices,
                                                                                    agglomerationInfo.SubCell3DsRemovedEdges,
                                                                                    agglomerationInfo.SubCell3DsRemovedFaces,
                                                                                    meshDao,
                                                                                    meshCell3DToConvexCell3DIndices);


      ASSERT_EQ(22,
                agglomeratedCell3DIndex);
    }



    Gedim::MeshUtilities::ExtractActiveMeshData activeMeshData;
    meshUtilities.ExtractActiveMesh(meshDao,
                                    activeMeshData);
    std::vector<std::vector<unsigned int>> activeMeshCell3DToConvexCell3DIndices(meshDao.Cell3DTotalNumber());
    for (unsigned int c3D_index = 0; c3D_index < meshDao.Cell3DTotalNumber(); c3D_index++)
    {
      activeMeshCell3DToConvexCell3DIndices.at(c3D_index) =
          meshCell3DToConvexCell3DIndices.at(activeMeshData.NewCell3DToOldCell3D.at(c3D_index));
    }

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "AgglomeratedMesh");
  }

  TEST(TestMeshUtilities, TestMakeMeshTriangularFaces)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                             1.0);
    const Eigen::Vector3d parallelepipedsLengthTangent = cube.Vertices.col(1) - cube.Vertices.col(0);
    const Eigen::Vector3d parallelepipedsHeightTangent = cube.Vertices.col(3) - cube.Vertices.col(0);
    const Eigen::Vector3d parallelepipedsWidthTangent = cube.Vertices.col(4) - cube.Vertices.col(0);

    const std::vector<double> lengthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(2 + 1,
                                                                                                        0.0, 1.0, 1);
    const std::vector<double> heightMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(1 + 1,
                                                                                                        0.0, 1.0, 1);
    const std::vector<double> widthMeshCurvilinearCoordinates = geometryUtilities.EquispaceCoordinates(2 + 1,
                                                                                                       0.0, 1.0, 1);


    const Eigen::Vector3d origin = cube.Vertices.col(0);
    meshUtilities.CreateParallelepipedMesh(origin,
                                           parallelepipedsLengthTangent,
                                           parallelepipedsHeightTangent,
                                           parallelepipedsWidthTangent,
                                           lengthMeshCurvilinearCoordinates,
                                           heightMeshCurvilinearCoordinates,
                                           widthMeshCurvilinearCoordinates,
                                           meshDao);

    std::string exportFolder = "./Export/TestMakeMeshTriangularFaces/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    Gedim::MeshMatrices new_mesh = mesh;
    Gedim::MeshMatricesDAO new_meshDao(new_mesh);
    std::vector<std::vector<unsigned int>> cell2Ds_triangles(meshDao.Cell2DTotalNumber());

    for (unsigned int f = 0; f < meshDao.Cell2DTotalNumber(); f++)
    {
      const Eigen::MatrixXd face_vertices =
          meshDao.Cell2DVerticesCoordinates(f);

      const auto face_translation =
          geometryUtilities.PolygonTranslation(face_vertices);
      const auto face_normal =
          geometryUtilities.PolygonNormal(face_vertices);
      const auto face_rotation =
          geometryUtilities.PolygonRotationMatrix(face_vertices,
                                                  face_normal,
                                                  face_translation);
      const auto face_vertices_2D =
          geometryUtilities.RotatePointsFrom3DTo2D(face_vertices,
                                                   face_rotation.transpose(),
                                                   face_translation);

      cell2Ds_triangles[f] =
          geometryUtilities.PolygonTriangulationByFirstVertex(face_vertices_2D);

    }

    meshUtilities.MakeMeshTriangularFaces(cell2Ds_triangles,
                                          new_meshDao);

    Gedim::MeshUtilities::ExtractActiveMeshData extract_data;
    meshUtilities.ExtractActiveMesh(new_meshDao,
                                    extract_data);

    meshUtilities.ExportMeshToVTU(new_meshDao,
                                  exportFolder,
                                  "TriangularFaceMesh");

    ASSERT_EQ(meshDao.Cell0DTotalNumber(),
              new_meshDao.Cell0DTotalNumber());
    ASSERT_EQ(meshDao.Cell1DTotalNumber() +
              meshDao.Cell2DTotalNumber(),
              new_meshDao.Cell1DTotalNumber());
    ASSERT_EQ(2 * meshDao.Cell2DTotalNumber(),
              new_meshDao.Cell2DTotalNumber());
    ASSERT_EQ(meshDao.Cell3DTotalNumber(),
              new_meshDao.Cell3DTotalNumber());

  }

  TEST(TestMeshUtilities, TestFindPointCell3D)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_3D_68Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    std::string exportFolder = "./Export/TestFindPointCell3D";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    const Gedim::MeshUtilities::MeshGeometricData3D mesh_geometry_data =
        meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                              meshDao);

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.0, 0.0, 1.1),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox);

      ASSERT_TRUE(result.Cell3Ds_found.empty());

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.empty());
    }

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.1, 0.1, 0.1),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox,
                                                        true,
                                                        48);

      ASSERT_TRUE(result.Cell3Ds_found.size() == 1);
      ASSERT_EQ(48, result.Cell3Ds_found[0].Cell3D_index);
      ASSERT_EQ(Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Inside,
                result.Cell3Ds_found[0].Cell3D_Position.Type);

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.size() == 1);
      ASSERT_EQ(Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types::Cell3D,
                result_position.MeshPositions[0].Type);
      ASSERT_EQ(48,
                result_position.MeshPositions[0].Cell_index);
    }

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.5, 0.2, 0.2),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox);

      ASSERT_TRUE(result.Cell3Ds_found.size() == 1);
      ASSERT_EQ(1, result.Cell3Ds_found[0].Cell3D_index);
      ASSERT_EQ(Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace,
                result.Cell3Ds_found[0].Cell3D_Position.Type);
      ASSERT_EQ(2,
                result.Cell3Ds_found[0].Cell3D_Position.BorderIndex);

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.size() == 1);
      ASSERT_EQ(Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types::Cell2D,
                result_position.MeshPositions[0].Type);
      ASSERT_EQ(5,
                result_position.MeshPositions[0].Cell_index);
    }

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.5, 0.25, 0.25),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox);

      ASSERT_TRUE(result.Cell3Ds_found.size() == 1);
      ASSERT_EQ(1, result.Cell3Ds_found[0].Cell3D_index);
      ASSERT_EQ(Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge,
                result.Cell3Ds_found[0].Cell3D_Position.Type);
      ASSERT_EQ(3,
                result.Cell3Ds_found[0].Cell3D_Position.BorderIndex);

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.size() == 1);
      ASSERT_EQ(Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types::Cell1D,
                result_position.MeshPositions[0].Type);
      ASSERT_EQ(6,
                result_position.MeshPositions[0].Cell_index);
    }

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.5, 0.5, 0.0),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox);

      ASSERT_TRUE(result.Cell3Ds_found.size() == 1);
      ASSERT_EQ(1, result.Cell3Ds_found[0].Cell3D_index);
      ASSERT_EQ(Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex,
                result.Cell3Ds_found[0].Cell3D_Position.Type);
      ASSERT_EQ(2,
                result.Cell3Ds_found[0].Cell3D_Position.BorderIndex);

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.size() == 1);
      ASSERT_EQ(Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types::Cell0D,
                result_position.MeshPositions[0].Type);
      ASSERT_EQ(44,
                result_position.MeshPositions[0].Cell_index);
    }

    {
      const auto result = meshUtilities.FindPointCell3D(geometryUtilities,
                                                        Eigen::Vector3d(0.5, 0.25, 0.25),
                                                        meshDao,
                                                        mesh_geometry_data.Cell3DsFaces,
                                                        mesh_geometry_data.Cell3DsFaces3DVertices,
                                                        mesh_geometry_data.Cell3DsFaces2DVertices,
                                                        mesh_geometry_data.Cell3DsFacesNormals,
                                                        mesh_geometry_data.Cell3DsFacesNormalDirections,
                                                        mesh_geometry_data.Cell3DsFacesTranslations,
                                                        mesh_geometry_data.Cell3DsFacesRotationMatrices,
                                                        mesh_geometry_data.Cell3DsBoundingBox,
                                                        false);

      ASSERT_TRUE(result.Cell3Ds_found.size() == 6);
      ASSERT_EQ(1, result.Cell3Ds_found[0].Cell3D_index);
      ASSERT_EQ(17, result.Cell3Ds_found[1].Cell3D_index);
      ASSERT_EQ(27, result.Cell3Ds_found[2].Cell3D_index);
      ASSERT_EQ(50, result.Cell3Ds_found[3].Cell3D_index);
      ASSERT_EQ(52, result.Cell3Ds_found[4].Cell3D_index);
      ASSERT_EQ(62, result.Cell3Ds_found[5].Cell3D_index);
      for (unsigned int i = 0; i < 6; i++)
        ASSERT_EQ(Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge,
                  result.Cell3Ds_found[i].Cell3D_Position.Type);

      const auto result_position = meshUtilities.FindPointMeshPosition(result,
                                                                       meshDao);
      ASSERT_TRUE(result_position.MeshPositions.size() == 6);
      for (unsigned int i = 0; i < 6; i++)
      {
        ASSERT_EQ(Gedim::MeshUtilities::FindPointMeshPositionResult::PointMeshPosition::Types::Cell1D,
                  result_position.MeshPositions[i].Type);
        ASSERT_EQ(6,
                  result_position.MeshPositions[i].Cell_index);
      }
    }
  }

  TEST(TestMeshUtilities, TestMarkCell3Ds)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestMarkCell3Ds";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    GedimUnitTesting::MeshMatrices_3D_68Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    auto marking_function = [&geometryUtilities](const Eigen::MatrixXd& points)
    {
      Eigen::VectorXi marks(points.cols());
      const Eigen::Vector3d center(0.5, 0.5, 0.5);

      for (unsigned int p = 0; p < points.cols(); p++)
      {
        const Eigen::Vector3d point = points.col(p);

        if (geometryUtilities.IsValueGreater((center - point).norm(),
                                             0.25,
                                             geometryUtilities.Tolerance1D()))
        {
          marks[p] = 0;
          continue;
        }

        if (geometryUtilities.IsValueGreater(point.x(),
                                             0.5,
                                             geometryUtilities.Tolerance1D()))
          marks[p] = 1;
        else
          marks[p] = 2;
      }

      return marks;
    };

    const auto mesh_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                           meshDao);

    std::vector<Eigen::MatrixXd> cell3Ds_points(meshDao.Cell3DTotalNumber());
    for (unsigned int c = 0; c < meshDao.Cell3DTotalNumber(); c++)
    {
      if (!meshDao.Cell3DIsActive(c))
        continue;

      cell3Ds_points[c].resize(3, 1);
      cell3Ds_points[c].col(0)<< mesh_geometric_data.Cell3DsCentroids.at(c);
    }

    const auto cell3Ds_marked = meshUtilities.MarkCells(marking_function,
                                                        cell3Ds_points,
                                                        0);

    ASSERT_EQ(std::vector<unsigned int>({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0}),
              cell3Ds_marked);

    {
      Gedim::VTKUtilities exporter;
      std::vector<double> marks(cell3Ds_marked.begin(),
                                cell3Ds_marked.end());

      exporter.AddPolyhedrons(meshDao.Cell0DsCoordinates(),
                              meshDao.Cell3DsFacesVertices(),
                              {
                                {
                                  "Mark",
                                  Gedim::VTPProperty::Formats::Cells,
                                  static_cast<unsigned int>(marks.size()),
                                  marks.data()
                                }
                              });
      exporter.Export(exportFolder + "/cell3Ds_mark.vtu");
    }
  }

  TEST(TestMeshUtilities, TestImportVtkMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    std::string exportFolder = "./Export/TestMeshUtilities/TestImportVtkMesh";
    Gedim::Output::CreateFolder(exportFolder);

    const std::string file_test_path = exportFolder +
                                       "/test_file.vtk";

    GedimUnitTesting::VTK_Unstructured_Grid_Mesh_Mock::ExportFile(file_test_path);

    Gedim::MeshUtilities meshUtilities;
    meshUtilities.ImportVtkMesh3D(file_test_path,
                                  meshDao);

#if ENABLE_VTK == 1
    ASSERT_EQ(5,
              meshDao.Cell0DTotalNumber());
    ASSERT_EQ(9,
              meshDao.Cell1DTotalNumber());
    ASSERT_EQ(7,
              meshDao.Cell2DTotalNumber());
    ASSERT_EQ(2,
              meshDao.Cell3DTotalNumber());
#endif

    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "ImportedMesh");
  }

  TEST(TestMeshUtilities, TestIntersect_mesh_polyhedron)
  {
    std::string exportFolder = "./Export/TestMeshUtilities/TestIntersect_mesh_polyhedron";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    // Create domain 3D
    const double domain_3D_length = 1.0;
    const double domain_3D_height = 1.0;
    const double domain_3D_width = 1.0;
    const unsigned int domain_3D_length_num_cells = 4;
    const unsigned int domain_3D_height_num_cells = 4;
    const unsigned int domain_3D_width_num_cells = 4;
    const Eigen::Vector3d domain_3D_origin(0.0, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_length_direction(domain_3D_length, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_height_direction(0.0, 0.0, domain_3D_height);
    const Eigen::Vector3d domain_3D_width_direction(0.0, domain_3D_width, 0.0);


    {
      const auto domain_3D = geometryUtilities.CreateParallelepipedWithOrigin(domain_3D_origin,
                                                                              domain_3D_length_direction,
                                                                              domain_3D_height_direction,
                                                                              domain_3D_width_direction);

      Gedim::VTKUtilities exporter;
      exporter.AddPolyhedron(domain_3D.Vertices,
                             domain_3D.Edges,
                             domain_3D.Faces);
      exporter.Export(exportFolder + "/Domain_3D.vtu");
    }

    // Create mesh 3D
    const std::vector<double> domain_3D_mesh_length_coordinates = geometryUtilities.EquispaceCoordinates(domain_3D_length_num_cells + 1,
                                                                                                         0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_height_coordinates = geometryUtilities.EquispaceCoordinates(domain_3D_height_num_cells + 1,
                                                                                                         0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_width_coordinates = geometryUtilities.EquispaceCoordinates(domain_3D_width_num_cells + 1,
                                                                                                        0.0, 1.0, 1);


    Gedim::MeshMatrices mesh_3D_data;
    Gedim::MeshMatricesDAO mesh_3D(mesh_3D_data);
    meshUtilities.CreateParallelepipedMesh(domain_3D_origin,
                                           domain_3D_length_direction,
                                           domain_3D_height_direction,
                                           domain_3D_width_direction,
                                           domain_3D_mesh_length_coordinates,
                                           domain_3D_mesh_height_coordinates,
                                           domain_3D_mesh_width_coordinates,
                                           mesh_3D);

    {
      std::string mesh_3D_export_folder = exportFolder + "/Mesh3D";
      Gedim::Output::CreateFolder(mesh_3D_export_folder);
      meshUtilities.ExportMeshToVTU(mesh_3D,
                                    mesh_3D_export_folder,
                                    "Mesh3D");
    }

    std::vector<Eigen::MatrixXd> cell2Ds_vertices(mesh_3D.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> cell2Ds_2D_vertices(mesh_3D.Cell2DTotalNumber());
    std::vector<Eigen::Vector3d> cell2Ds_normal(mesh_3D.Cell2DTotalNumber());
    std::vector<Eigen::Vector3d> cell2Ds_translation(mesh_3D.Cell2DTotalNumber());
    std::vector<Eigen::Matrix3d> cell2Ds_rotation_matrix(mesh_3D.Cell2DTotalNumber());
    std::vector<Eigen::MatrixXd> cell2Ds_bounding_box(mesh_3D.Cell2DTotalNumber());
    for (unsigned int p = 0; p < mesh_3D.Cell2DTotalNumber(); p++)
    {
      cell2Ds_vertices[p] = mesh_3D.Cell2DVerticesCoordinates(p);
      cell2Ds_normal[p] = geometryUtilities.PolygonNormal(cell2Ds_vertices[p]);
      cell2Ds_bounding_box[p] = geometryUtilities.PointsBoundingBox(cell2Ds_vertices[p]);
      cell2Ds_translation[p] = geometryUtilities.PolygonTranslation(cell2Ds_vertices[p]);
      cell2Ds_rotation_matrix[p] = geometryUtilities.PolygonRotationMatrix(cell2Ds_vertices[p],
                                                                           cell2Ds_normal[p],
                                                                           cell2Ds_translation[p]);
      cell2Ds_2D_vertices[p] = geometryUtilities.RotatePointsFrom3DTo2D(cell2Ds_vertices[p],
                                                                        cell2Ds_rotation_matrix[p].transpose(),
                                                                        cell2Ds_translation[p]);
    }


    const auto cell1Ds_geometric_data = meshUtilities.FillMesh1DGeometricData(geometryUtilities,
                                                                              mesh_3D);

    const auto cell3Ds_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                              mesh_3D);


    //    const Gedim::PlatonicSolid platonicSolid = Gedim::PlatonicSolid(geometryUtilities,
    //                                                                    meshUtilities);

    //        const Gedim::GeometryUtilities::Polyhedron ico = platonicSolid.icosahedron();
    //         Gedim::GeometryUtilities::Polyhedron dual = platonicSolid.first_class_geodesic_polyhedron(ico, 4);
    //        Gedim::GeometryUtilities::Polyhedron polyhedron = platonicSolid.goldberg_polyhedron(dual);

    //    const Gedim::SphereMeshUtilities sphereMeshUtilities(geometryUtilities,
    //                                                         meshUtilities);

    //    Gedim::GeometryUtilities::Polyhedron polyhedron = sphereMeshUtilities.uv_sphere(10, 7);

    //    polyhedron.Vertices = (0.1 * polyhedron.Vertices).colwise() +
    //                          Eigen::Vector3d(0.3, 0.3, 0.3);

    const auto polyhedron = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.25, 0.25, 0.25),
                                                                             Eigen::Vector3d(0.1, 0.1, 0.0),
                                                                             Eigen::Vector3d(-0.1, 0.1, 0.0),
                                                                             Eigen::Vector3d(0.0, 0.0, 0.1));

    {
      Gedim::VTKUtilities exporter;

      exporter.AddPolyhedron(polyhedron.Vertices,
                             polyhedron.Edges,
                             polyhedron.Faces);

      exporter.Export(exportFolder + "/polyhedron.vtu");
    }

    std::vector<Eigen::MatrixXd> polyhedron_edges_vertices(polyhedron.Edges.cols());
    std::vector<Eigen::MatrixXd> polyhedron_edges_bounding_box(polyhedron.Edges.cols());
    for (unsigned int e = 0; e < polyhedron.Edges.cols(); e++)
    {
      polyhedron_edges_vertices[e].resize(3, 2);
      polyhedron_edges_vertices[e].col(0)<< polyhedron.Vertices.col(polyhedron.Edges(0, e));
      polyhedron_edges_vertices[e].col(1)<< polyhedron.Vertices.col(polyhedron.Edges(1, e));

      polyhedron_edges_bounding_box[e] = geometryUtilities.PointsBoundingBox(polyhedron_edges_vertices[e]);
    }

    const auto polyhedron_edges_tangent = geometryUtilities.PolyhedronEdgeTangents(polyhedron.Vertices,
                                                                                   polyhedron.Edges);

    const auto polyhedron_faces_vertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                    polyhedron.Faces);

    std::vector<Eigen::MatrixXd> polyhedron_faces_bounding_box(polyhedron.Faces.size());
    for (unsigned int f = 0; f < polyhedron.Faces.size(); f++)
    {
      polyhedron_faces_bounding_box[f] = geometryUtilities.PointsBoundingBox(polyhedron_faces_vertices[f]);
    }

    const auto polyhedron_internal_point = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
    const auto polyhedron_faces_translation = geometryUtilities.PolyhedronFaceTranslations(polyhedron_faces_vertices);
    const auto polyhedron_faces_normal = geometryUtilities.PolyhedronFaceNormals(polyhedron_faces_vertices);
    const auto polyhedron_faces_normal_direction = geometryUtilities.PolyhedronFaceNormalDirections(polyhedron_faces_vertices,
                                                                                                    polyhedron_internal_point,
                                                                                                    polyhedron_faces_normal);
    const auto polyhedron_faces_rotation_matrix = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedron_faces_vertices,
                                                                                                   polyhedron_faces_normal,
                                                                                                   polyhedron_faces_translation);
    const auto polyhedron_faces_rotated_vertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedron_faces_vertices,
                                                                                                   polyhedron_faces_translation,
                                                                                                   polyhedron_faces_rotation_matrix);
    const auto polyhedron_bounding_box = geometryUtilities.PointsBoundingBox(polyhedron.Vertices);

    const auto result = meshUtilities.Intersect_mesh_polyhedron(geometryUtilities,
                                                                polyhedron.Vertices,
                                                                polyhedron.Edges,
                                                                polyhedron_edges_vertices,
                                                                polyhedron_edges_tangent,
                                                                polyhedron_edges_bounding_box,
                                                                polyhedron.Faces,
                                                                polyhedron_faces_vertices,
                                                                polyhedron_faces_rotated_vertices,
                                                                polyhedron_faces_normal,
                                                                polyhedron_faces_normal_direction,
                                                                polyhedron_faces_translation,
                                                                polyhedron_faces_rotation_matrix,
                                                                polyhedron_faces_bounding_box,                                                                polyhedron_bounding_box,
                                                                mesh_3D,
                                                                cell1Ds_geometric_data.Cell1DsBoundingBox,
                                                                cell1Ds_geometric_data.Cell1DsVertices,
                                                                cell1Ds_geometric_data.Cell1DsTangents,
                                                                cell2Ds_vertices,
                                                                cell2Ds_normal,
                                                                cell2Ds_2D_vertices,
                                                                cell2Ds_translation,
                                                                cell2Ds_rotation_matrix,
                                                                cell2Ds_bounding_box,
                                                                cell3Ds_geometric_data.Cell3DsBoundingBox,
                                                                cell3Ds_geometric_data.Cell3DsFaces,
                                                                cell3Ds_geometric_data.Cell3DsFaces3DVertices,
                                                                cell3Ds_geometric_data.Cell3DsFaces2DVertices,
                                                                cell3Ds_geometric_data.Cell3DsFacesNormals,
                                                                cell3Ds_geometric_data.Cell3DsFacesNormalDirections,
                                                                cell3Ds_geometric_data.Cell3DsFacesTranslations,
                                                                cell3Ds_geometric_data.Cell3DsFacesRotationMatrices);

    if (result.Type ==
        Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Types::Vertices)
    {
      Gedim::VTKUtilities exporter;

      exporter.AddPoints(result.Intersections_Coordinates);

      exporter.Export(exportFolder + "/intersections.vtu");
    }

    ASSERT_EQ(Gedim::MeshUtilities::Intersect_mesh_polyhedron_result::Types::Vertices,
              result.Type);
    ASSERT_EQ(8,
              result.Intersections_Coordinates.cols());
    ASSERT_EQ(8,
              result.Polyhedron_intersections.size());

  }
}

#endif // __TEST_MESH_UTILITIES3D_H
