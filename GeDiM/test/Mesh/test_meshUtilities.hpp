#ifndef __TEST_MESH_UTILITIES_H
#define __TEST_MESH_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_2D_26Cells_Mock.hpp"

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_3D_1Cells_Mock.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestMeshUtilities, TestMesh1DFromSegment)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    Eigen::MatrixXd segment(3, 4);
    segment.col(0)<< 0.0, 0.0, 0.0;
    segment.col(1)<< 1.0, 0.0, 0.0;
    vector<unsigned int> vertexMarkers = { 1, 2 };

    meshUtilities.Mesh1DFromSegment(geometryUtilities,
                                    segment,
                                    vertexMarkers,
                                    meshDao);

    std::string exportFolder = "./Export/TestMesh1DFromSegment/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 1);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 2);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 1);
  }

  TEST(TestMeshUtilities, TestMesh2DFromPolygon)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    Eigen::MatrixXd polygon(3, 4);
    polygon.col(0)<< 0.0, 0.0, 0.0;
    polygon.col(1)<< 1.0, 0.0, 0.0;
    polygon.col(2)<< 1.0, 1.0, 0.0;
    polygon.col(3)<< 0.0, 1.0, 0.0;
    vector<unsigned int> vertexMarkers = { 1, 2, 3, 4 };
    vector<unsigned int> edgeMarkers = { 5, 6, 7, 8 };

    meshUtilities.Mesh2DFromPolygon(polygon,
                                    vertexMarkers,
                                    edgeMarkers,
                                    meshDao);

    std::string exportFolder = "./Export/TestMesh2DFromPolygon/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 1);
  }

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

  TEST(TestMeshUtilities, TestCreateRectangleMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::Vector3d rectangleOrigin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBaseTangent = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeightTangent = Eigen::Vector3d(0.0, 1.0, 0.0);
    const vector<double> baseCoordinates = geometryUtilities.EquispaceCoordinates(0.25, true);
    const vector<double> heightCoordinates = geometryUtilities.EquispaceCoordinates(0.5, true);

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseCoordinates,
                                      heightCoordinates,
                                      meshDao);

    Eigen::MatrixXd cell0DCoordinates(3, 15);
    cell0DCoordinates.col(0)<< 0.0, 0.0, 0.0;
    cell0DCoordinates.col(1)<< 0.25, 0.0, 0.0;
    cell0DCoordinates.col(2)<< 0.5, 0.0, 0.0;
    cell0DCoordinates.col(3)<< 0.75, 0.0, 0.0;
    cell0DCoordinates.col(4)<< 1.0, 0.0, 0.0;
    cell0DCoordinates.col(5)<< 0.0, 0.5, 0.0;
    cell0DCoordinates.col(6)<< 0.25, 0.5, 0.0;
    cell0DCoordinates.col(7)<< 0.5, 0.5, 0.0;
    cell0DCoordinates.col(8)<< 0.75, 0.5, 0.0;
    cell0DCoordinates.col(9)<< 1.0, 0.5, 0.0;
    cell0DCoordinates.col(10)<< 0.0, 1.0, 0.0;
    cell0DCoordinates.col(11)<< 0.25, 1.0, 0.0;
    cell0DCoordinates.col(12)<< 0.5, 1.0, 0.0;
    cell0DCoordinates.col(13)<< 0.75, 1.0, 0.0;
    cell0DCoordinates.col(14)<< 1.0, 1.0, 0.0;

    Eigen::MatrixXi cell1DExtremes(2, 22);
    cell1DExtremes.col(0)<< 0, 1;
    cell1DExtremes.col(1)<< 1, 2;
    cell1DExtremes.col(2)<< 2, 3;
    cell1DExtremes.col(3)<< 3, 4;
    cell1DExtremes.col(4)<< 5, 6;
    cell1DExtremes.col(5)<< 6, 7;
    cell1DExtremes.col(6)<< 7, 8;
    cell1DExtremes.col(7)<< 8, 9;
    cell1DExtremes.col(8)<< 10, 11;
    cell1DExtremes.col(9)<< 11, 12;
    cell1DExtremes.col(10)<< 12, 13;
    cell1DExtremes.col(11)<< 13, 14;
    cell1DExtremes.col(12)<< 0, 5;
    cell1DExtremes.col(13)<< 1, 6;
    cell1DExtremes.col(14)<< 2, 7;
    cell1DExtremes.col(15)<< 3, 8;
    cell1DExtremes.col(16)<< 4, 9;
    cell1DExtremes.col(17)<< 5, 10;
    cell1DExtremes.col(18)<< 6, 11;
    cell1DExtremes.col(19)<< 7, 12;
    cell1DExtremes.col(20)<< 8, 13;
    cell1DExtremes.col(21)<< 9, 14;

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 15);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 22);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 8);
    EXPECT_EQ(meshDao.Cell0DsCoordinates(),
              cell0DCoordinates);
    EXPECT_EQ(meshDao.Cell1DsExtremes(),
              cell1DExtremes);
  }

  TEST(TestMeshUtilities, TestFillMesh2D_TriangularMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO expectedMesh(mockMesh.Mesh);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::MatrixXd vertices = (Eigen::MatrixXd(3, 5)<<
                                      0.0, 1.0, 1.0, 0.0, 0.5,
                                      0.0, 0.0, 1.0, 1.0, 0.5,
                                      0.0, 0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 8)<<
                                   1, 2, 4, 3, 0, 4, 2, 0,
                                   2, 4, 1, 0, 4, 3, 3, 1).finished();

    vector<Eigen::MatrixXi> triangles(4, Eigen::MatrixXi(2, 3));
    triangles[0] = (Eigen::MatrixXi(2, 3)<<
                    1,2,4,
                    0,1,2).finished();
    triangles[1] = (Eigen::MatrixXi(2, 3)<<
                    3,0,4,
                    3,4,5).finished();
    triangles[2] = (Eigen::MatrixXi(2, 3)<<
                    4,2,3,
                    1,6,5).finished();
    triangles[3] = (Eigen::MatrixXi(2, 3)<<
                    0,1,4,
                    7,2,4).finished();

    meshUtilities.FillMesh2D(vertices,
                             edges,
                             triangles,
                             meshDao);

    ASSERT_EQ(expectedMesh.Cell0DsCoordinates(),
              meshDao.Cell0DsCoordinates());
    ASSERT_EQ(expectedMesh.Cell1DsExtremes(),
              meshDao.Cell1DsExtremes());
    for (unsigned int c = 0; c < 4; c++)
    {
      ASSERT_EQ(expectedMesh.Cell2DVertices(c),
                meshDao.Cell2DVertices(c));
      ASSERT_EQ(expectedMesh.Cell2DEdges(c),
                meshDao.Cell2DEdges(c));
    }

  }

  TEST(TestMeshUtilities, TestCreateTriangleMesh)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO expectedMesh(mockMesh.Mesh);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Gedim::MeshUtilities meshUtilities;

    const Eigen::MatrixXd polygon = geometryUtilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                   1.0);

    meshUtilities.CreateTriangularMesh(polygon,
                                       0.033,
                                       meshDao);

    std::string exportFolder = "./Export/TestMeshUtilities/TestCreateTriangleMesh";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao,
                                  exportFolder,
                                  "CreatedTriangleMesh");
    meshUtilities.ExportMeshToVTU(expectedMesh,
                                  exportFolder,
                                  "ExpectedTriangleMesh");

    Gedim::MeshUtilities::MeshGeometricData2D cell2DsGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                           meshDao);
    const unsigned int cell2DToExportIndex = 0;
    meshUtilities.ExportCell2DToVTU(meshDao,
                                    cell2DToExportIndex,
                                    cell2DsGeometricData.Cell2DsVertices[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsTriangulations[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsAreas[cell2DToExportIndex],
                                    cell2DsGeometricData.Cell2DsCentroids[cell2DToExportIndex],
                                    exportFolder);


    EXPECT_EQ(expectedMesh.Dimension(),
              meshDao.Dimension());
    EXPECT_EQ(expectedMesh.Cell0DTotalNumber(),
              meshDao.Cell0DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell1DTotalNumber(),
              meshDao.Cell1DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell2DTotalNumber(),
              meshDao.Cell2DTotalNumber());
    EXPECT_EQ(expectedMesh.Cell0DsCoordinates(),
              meshDao.Cell0DsCoordinates());
    EXPECT_EQ(expectedMesh.Cell1DsExtremes(),
              meshDao.Cell1DsExtremes());

    const vector<unsigned int> cell0DMarkers = { 11, 12, 13, 14 };
    const vector<unsigned int> cell1DMarkers = { 15, 16, 17, 18 };
    meshUtilities.ChangePolygonMeshMarkers(polygon,
                                           cell0DMarkers,
                                           cell1DMarkers,
                                           meshDao);

    EXPECT_EQ(cell0DMarkers[0],
        meshDao.Cell0DMarker(0));
    EXPECT_EQ(cell0DMarkers[1],
        meshDao.Cell0DMarker(1));
    EXPECT_EQ(cell0DMarkers[2],
        meshDao.Cell0DMarker(2));
    EXPECT_EQ(cell0DMarkers[3],
        meshDao.Cell0DMarker(3));
    EXPECT_EQ(cell1DMarkers[1],
        meshDao.Cell1DMarker(20));
    EXPECT_EQ(cell1DMarkers[0],
        meshDao.Cell1DMarker(23));
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
                                        meshDao);

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

  TEST(TestMeshUtilities, TestCheckMesh2D)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshUtilities::CheckMesh2DConfiguration config;
    ASSERT_NO_THROW(meshUtilities.CheckMesh2D(config,
                                              geometryUtilities,
                                              meshDao));
  }

  TEST(TestMeshUtilities, TestComputeCell1DCell2DNeighbours)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    mesh.Mesh.NumberCell1DNeighbourCell2D.clear();
    mesh.Mesh.Cell1DNeighbourCell2Ds.clear();
    mesh.Mesh.NumberCell1DNeighbourCell2D.resize(mesh.Mesh.NumberCell1D + 1, 0);
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);
    Gedim::MeshUtilities meshUtilities;

    meshUtilities.ComputeCell1DCell2DNeighbours(meshDao);

    EXPECT_EQ(mesh.Mesh.NumberCell1DNeighbourCell2D,
              vector<unsigned int>({ 0,2,4,6,8 }));
    EXPECT_EQ(mesh.Mesh.Cell1DNeighbourCell2Ds,
              vector<unsigned int>({ 1,0,1,0,1,0,1,0 }));
  }

  TEST(TestMeshUtilities, TestFillMesh1DGeometricData)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData2D result = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                   meshDao);

    Gedim::MeshUtilities::MeshGeometricData2D expectedResult;
    expectedResult.Cell2DsAreas = { 1.0 };
    expectedResult.Cell2DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.0) };
    expectedResult.Cell2DsDiameters = { sqrt(2.0) };
    expectedResult.Cell2DsEdgeDirections = { { true, true, true, true } };
    Eigen::VectorXd edgeLengths(4);
    edgeLengths<< 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = { edgeLengths };
    Eigen::MatrixXd edgeNormals(3, 4);
    edgeNormals.col(0)<< -1.0, 0.0, 0.0;
    edgeNormals.col(1)<< 0.0, -1.0, 0.0;
    edgeNormals.col(2)<< 1.0, 0.0, 0.0;
    edgeNormals.col(3)<< 0.0, 1.0, 0.0;
    expectedResult.Cell2DsEdgeNormals = { edgeNormals };
    Eigen::MatrixXd edgeTangents(3, 4);
    edgeTangents.col(0)<< 0.0, -1.0, 0.0;
    edgeTangents.col(1)<< 1.0, 0.0, 0.0;
    edgeTangents.col(2)<< 0.0, 1.0, 0.0;
    edgeTangents.col(3)<< -1.0, 0.0, 0.0;
    expectedResult.Cell2DsEdgeTangents = { edgeTangents };
    Eigen::Matrix3d triangleOne;
    triangleOne.col(0)<< 0.0, 1.0, 0.0;
    triangleOne.col(1)<< 0.0, 0.0, 0.0;
    triangleOne.col(2)<< 1.0, 0.0, 0.0;
    Eigen::Matrix3d triangleTwo;
    triangleTwo.col(0)<< 0.0, 1.0, 0.0;
    triangleTwo.col(1)<< 1.0, 0.0, 0.0;
    triangleTwo.col(2)<< 1.0, 1.0, 0.0;
    expectedResult.Cell2DsTriangulations = { { triangleOne, triangleTwo } };
    Eigen::MatrixXd vertices(3, 4);
    vertices.col(0)<< 0.0, 1.0, 0.0;
    vertices.col(1)<< 0.0, 0.0, 0.0;
    vertices.col(2)<< 1.0, 0.0, 0.0;
    vertices.col(3)<< 1.0, 1.0, 0.0;
    expectedResult.Cell2DsVertices = { vertices };

    EXPECT_EQ(result.Cell2DsAreas, expectedResult.Cell2DsAreas);
    EXPECT_EQ(result.Cell2DsCentroids, expectedResult.Cell2DsCentroids);
    EXPECT_EQ(result.Cell2DsDiameters, expectedResult.Cell2DsDiameters);
    EXPECT_EQ(result.Cell2DsEdgeDirections, expectedResult.Cell2DsEdgeDirections);
    EXPECT_EQ(result.Cell2DsEdgeLengths, expectedResult.Cell2DsEdgeLengths);
    EXPECT_EQ(result.Cell2DsEdgeNormals, expectedResult.Cell2DsEdgeNormals);
    EXPECT_EQ(result.Cell2DsEdgeTangents, expectedResult.Cell2DsEdgeTangents);
    EXPECT_EQ(result.Cell2DsTriangulations, expectedResult.Cell2DsTriangulations);
    EXPECT_EQ(result.Cell2DsVertices, expectedResult.Cell2DsVertices);
  }

  TEST(TestMeshUtilities, TestFillMesh1DGeometricData_Convex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Eigen::MatrixXd segment(3, 4);
    segment.col(0)<< 0.0, 0.0, 0.0;
    segment.col(1)<< 1.0, 0.0, 0.0;
    vector<unsigned int> vertexMarkers = { 1, 2 };

    meshUtilities.Mesh1DFromSegment(geometryUtilities,
                                    segment,
                                    vertexMarkers,
                                    meshDao);

    const Gedim::MeshUtilities::MeshGeometricData1D result = meshUtilities.FillMesh1DGeometricData(geometryUtilities,
                                                                                                   meshDao);

    Gedim::MeshUtilities::MeshGeometricData1D expectedResult;
    expectedResult.Cell1DsLengths = { 1.0 };
    expectedResult.Cell1DsSquaredLengths = { 1.0 };
    expectedResult.Cell1DsCentroids = { Eigen::Vector3d(0.5, 0.0, 0.0) };
    expectedResult.Cell1DsTangents = { Eigen::Vector3d(1.0, 0.0, 0.0) };
    expectedResult.Cell1DsVertices = { (Eigen::MatrixXd(3, 2)<< 0.0, 1.0, 0.0, 0.0, 0.0, 0.0).finished() };

    EXPECT_EQ(result.Cell1DsLengths, expectedResult.Cell1DsLengths);
    EXPECT_EQ(result.Cell1DsSquaredLengths, expectedResult.Cell1DsSquaredLengths);
    EXPECT_EQ(result.Cell1DsCentroids, expectedResult.Cell1DsCentroids);
    EXPECT_EQ(result.Cell1DsTangents, expectedResult.Cell1DsTangents);
    EXPECT_EQ(result.Cell1DsVertices, expectedResult.Cell1DsVertices);
  }

  TEST(TestMeshUtilities, TestFillMesh2DGeometricData_NonConvex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    GedimUnitTesting::MeshMatrices_2D_2Cells_Mock convexMesh;
    Gedim::MeshMatricesDAO convexMeshDao(convexMesh.Mesh);

    vector<vector<unsigned int>> meshCell2DToConvexCell2DIndices = { { 0, 1 } };

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData2D result = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                   meshDao,
                                                                                                   convexMeshDao,
                                                                                                   meshCell2DToConvexCell2DIndices);

    Gedim::MeshUtilities::MeshGeometricData2D expectedResult;
    expectedResult.Cell2DsAreas = { 1.0 };
    expectedResult.Cell2DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.0) };
    expectedResult.Cell2DsDiameters = { sqrt(2.0) };
    expectedResult.Cell2DsEdgeDirections = { { true, true, true, true } };
    Eigen::VectorXd edgeLengths(4);
    edgeLengths<< 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = { edgeLengths };
    Eigen::MatrixXd edgeNormals(3, 4);
    edgeNormals.col(0)<< -1.0, 0.0, 0.0;
    edgeNormals.col(1)<< 0.0, -1.0, 0.0;
    edgeNormals.col(2)<< 1.0, 0.0, 0.0;
    edgeNormals.col(3)<< 0.0, 1.0, 0.0;
    expectedResult.Cell2DsEdgeNormals = { edgeNormals };
    Eigen::MatrixXd edgeTangents(3, 4);
    edgeTangents.col(0)<< 0.0, -1.0, 0.0;
    edgeTangents.col(1)<< 1.0, 0.0, 0.0;
    edgeTangents.col(2)<< 0.0, 1.0, 0.0;
    edgeTangents.col(3)<< -1.0, 0.0, 0.0;
    expectedResult.Cell2DsEdgeTangents = { edgeTangents };
    Eigen::Matrix3d triangleOne;
    triangleOne.col(0)<< 0.0, 1.0, 0.0;
    triangleOne.col(1)<< 0.0, 0.0, 0.0;
    triangleOne.col(2)<< 1.0, 0.0, 0.0;
    Eigen::Matrix3d triangleTwo;
    triangleTwo.col(0)<< 1.0, 0.0, 0.0;
    triangleTwo.col(1)<< 1.0, 1.0, 0.0;
    triangleTwo.col(2)<< 0.0, 1.0, 0.0;
    expectedResult.Cell2DsTriangulations = { { triangleOne, triangleTwo } };
    Eigen::MatrixXd vertices(3, 4);
    vertices.col(0)<< 0.0, 1.0, 0.0;
    vertices.col(1)<< 0.0, 0.0, 0.0;
    vertices.col(2)<< 1.0, 0.0, 0.0;
    vertices.col(3)<< 1.0, 1.0, 0.0;
    expectedResult.Cell2DsVertices = { vertices };

    EXPECT_EQ(result.Cell2DsAreas, expectedResult.Cell2DsAreas);
    EXPECT_EQ(result.Cell2DsCentroids, expectedResult.Cell2DsCentroids);
    EXPECT_EQ(result.Cell2DsDiameters, expectedResult.Cell2DsDiameters);
    EXPECT_EQ(result.Cell2DsEdgeDirections, expectedResult.Cell2DsEdgeDirections);
    EXPECT_EQ(result.Cell2DsEdgeLengths, expectedResult.Cell2DsEdgeLengths);
    EXPECT_EQ(result.Cell2DsEdgeNormals, expectedResult.Cell2DsEdgeNormals);
    EXPECT_EQ(result.Cell2DsEdgeTangents, expectedResult.Cell2DsEdgeTangents);
    EXPECT_EQ(result.Cell2DsTriangulations, expectedResult.Cell2DsTriangulations);
    EXPECT_EQ(result.Cell2DsVertices, expectedResult.Cell2DsVertices);
  }

  TEST(TestMeshUtilities, TestFillMesh3DGeometricData_Convex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-14;
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

    EXPECT_EQ(result.Cell3DsVertices, expectedResult.Cell3DsVertices);
    EXPECT_EQ(result.Cell3DsEdges, expectedResult.Cell3DsEdges);
    EXPECT_EQ(result.Cell3DsFaces, expectedResult.Cell3DsFaces);
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsVolumes[0], expectedResult.Cell3DsVolumes[0]));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].x(), expectedResult.Cell3DsCentroids[0].z()));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].y(), expectedResult.Cell3DsCentroids[0].y()));
    EXPECT_TRUE(geometryUtilities.Are1DValuesEqual(result.Cell3DsCentroids[0].z(), expectedResult.Cell3DsCentroids[0].x()));
    EXPECT_EQ(result.Cell3DsDiameters, expectedResult.Cell3DsDiameters);
  }
}

#endif // __TEST_MESH_UTILITIES_H
