#ifndef __TEST_MESH_UTILITIES3D_H
#define __TEST_MESH_UTILITIES3D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

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
#include "VTKUtilities.hpp"

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
}

#endif // __TEST_MESH_UTILITIES3D_H
