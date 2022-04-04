#ifndef __TEST_MESH_UTILITIES_H
#define __TEST_MESH_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestMeshUtilities, TestMeshMatricesDAO)
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

    EXPECT_EQ(meshDao.Dimension(), 2);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 4);
    EXPECT_EQ(meshDao.Cell2DTotalNumber(), 1);
  }

  TEST(TestMeshUtilities, TestFillMesh2DGeometricData_Convex)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData result = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                 meshDao);
    Gedim::MeshUtilities::MeshGeometricData expectedResult;
    expectedResult.Cell2DsAreas = { 1.0 };
    expectedResult.Cell2DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.0) };
    expectedResult.Cell2DsDiameters = { sqrt(2.0) };
    expectedResult.Cell2DsEdgeDirections = { { true, true, true, true } };
    Eigen::VectorXd edgeLenghts(4);
    edgeLenghts<< 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = { edgeLenghts };
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

    const Gedim::MeshUtilities::MeshGeometricData result = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                 meshDao,
                                                                                                 convexMeshDao,
                                                                                                 meshCell2DToConvexCell2DIndices);
    Gedim::MeshUtilities::MeshGeometricData expectedResult;
    expectedResult.Cell2DsAreas = { 1.0 };
    expectedResult.Cell2DsCentroids = { Eigen::Vector3d(0.5, 0.5, 0.0) };
    expectedResult.Cell2DsDiameters = { sqrt(2.0) };
    expectedResult.Cell2DsEdgeDirections = { { true, true, true, true } };
    Eigen::VectorXd edgeLenghts(4);
    edgeLenghts<< 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = { edgeLenghts };
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
}

#endif // __TEST_MESH_UTILITIES_H
