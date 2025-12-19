// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __TEST_MESH_UTILITIES1D_H
#define __TEST_MESH_UTILITIES1D_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshUtilities.hpp"

#include "MeshDAOImporterFromCsv.hpp"
#include "MeshFromCsvUtilities.hpp"

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
    segment.col(0) << 0.0, 0.0, 0.0;
    segment.col(1) << 1.0, 0.0, 0.0;
    vector<unsigned int> vertexMarkers = {1, 2};

    meshUtilities.Mesh1DFromSegment(geometryUtilities, segment, vertexMarkers, meshDao);

    std::string exportFolder = "./Export/TestMesh1DFromSegment/";
    Gedim::Output::CreateFolder(exportFolder);
    meshUtilities.ExportMeshToVTU(meshDao, exportFolder, "Mesh");

    EXPECT_EQ(meshDao.Dimension(), 1);
    EXPECT_EQ(meshDao.Cell0DTotalNumber(), 2);
    EXPECT_EQ(meshDao.Cell1DTotalNumber(), 1);
}

TEST(TestMeshUtilities, TestFillMesh1DGeometricData)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    GedimUnitTesting::MeshMatrices_2D_1Cells_Mock mesh;
    Gedim::MeshMatricesDAO meshDao(mesh.Mesh);

    Gedim::MeshUtilities meshUtilities;

    const Gedim::MeshUtilities::MeshGeometricData2D result = meshUtilities.FillMesh2DGeometricData(geometryUtilities, meshDao);

    Gedim::MeshUtilities::MeshGeometricData2D expectedResult;
    expectedResult.Cell2DsAreas = {1.0};
    expectedResult.Cell2DsCentroids = {Eigen::Vector3d(0.5, 0.5, 0.0)};
    expectedResult.Cell2DsDiameters = {sqrt(2.0)};
    expectedResult.Cell2DsEdgeDirections = {{true, true, true, true}};
    Eigen::VectorXd edgeLengths(4);
    edgeLengths << 1.0, 1.0, 1.0, 1.0;
    expectedResult.Cell2DsEdgeLengths = {edgeLengths};
    Eigen::MatrixXd edgeNormals(3, 4);
    edgeNormals.col(0) << -1.0, 0.0, 0.0;
    edgeNormals.col(1) << 0.0, -1.0, 0.0;
    edgeNormals.col(2) << 1.0, 0.0, 0.0;
    edgeNormals.col(3) << 0.0, 1.0, 0.0;
    expectedResult.Cell2DsEdgeNormals = {edgeNormals};
    Eigen::MatrixXd edgeTangents(3, 4);
    edgeTangents.col(0) << 0.0, -1.0, 0.0;
    edgeTangents.col(1) << 1.0, 0.0, 0.0;
    edgeTangents.col(2) << 0.0, 1.0, 0.0;
    edgeTangents.col(3) << -1.0, 0.0, 0.0;
    expectedResult.Cell2DsEdgeTangents = {edgeTangents};
    Eigen::Matrix3d triangleOne;
    triangleOne.col(0) << 0.0, 1.0, 0.0;
    triangleOne.col(1) << 0.0, 0.0, 0.0;
    triangleOne.col(2) << 1.0, 0.0, 0.0;
    Eigen::Matrix3d triangleTwo;
    triangleTwo.col(0) << 0.0, 1.0, 0.0;
    triangleTwo.col(1) << 1.0, 0.0, 0.0;
    triangleTwo.col(2) << 1.0, 1.0, 0.0;
    expectedResult.Cell2DsTriangulations = {{triangleOne, triangleTwo}};
    Eigen::MatrixXd vertices(3, 4);
    vertices.col(0) << 0.0, 1.0, 0.0;
    vertices.col(1) << 0.0, 0.0, 0.0;
    vertices.col(2) << 1.0, 0.0, 0.0;
    vertices.col(3) << 1.0, 1.0, 0.0;
    expectedResult.Cell2DsVertices = {vertices};

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
    std::string exportFolder = "./Export/TestMeshUtilities/TestFillMesh1DGeometricData_Convex";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDao(mesh);

    Eigen::MatrixXd segment(3, 4);
    segment.col(0) << 0.0, 0.0, 0.0;
    segment.col(1) << 1.0, 0.0, 0.0;
    vector<unsigned int> vertexMarkers = {1, 2};

    meshUtilities.Mesh1DFromSegment(geometryUtilities, segment, vertexMarkers, meshDao);

    const Gedim::MeshUtilities::MeshGeometricData1D result = meshUtilities.FillMesh1DGeometricData(geometryUtilities, meshDao);

    Gedim::MeshUtilities::MeshGeometricData1D expectedResult;
    expectedResult.Cell1DsLengths = {1.0};
    expectedResult.Cell1DsSquaredLengths = {1.0};
    expectedResult.Cell1DsCentroids = {Eigen::Vector3d(0.5, 0.0, 0.0)};
    expectedResult.Cell1DsTangents = {Eigen::Vector3d(1.0, 0.0, 0.0)};
    expectedResult.Cell1DsVertices = {(Eigen::MatrixXd(3, 2) << 0.0, 1.0, 0.0, 0.0, 0.0, 0.0).finished()};

    EXPECT_EQ(result.Cell1DsLengths, expectedResult.Cell1DsLengths);
    EXPECT_EQ(result.Cell1DsSquaredLengths, expectedResult.Cell1DsSquaredLengths);
    EXPECT_EQ(result.Cell1DsCentroids, expectedResult.Cell1DsCentroids);
    EXPECT_EQ(result.Cell1DsTangents, expectedResult.Cell1DsTangents);
    EXPECT_EQ(result.Cell1DsVertices, expectedResult.Cell1DsVertices);

    meshUtilities.ExportMeshGeometricData1DToTxt(result, exportFolder + "/geometric_properties.txt");
    const auto imported_result = meshUtilities.ImportMeshGeometricData1DFromTxt(exportFolder + "/geometric_properties."
                                                                                               "txt");

    EXPECT_EQ(imported_result.Cell1DsLengths, expectedResult.Cell1DsLengths);
    EXPECT_EQ(imported_result.Cell1DsSquaredLengths, expectedResult.Cell1DsSquaredLengths);
    EXPECT_EQ(imported_result.Cell1DsCentroids, expectedResult.Cell1DsCentroids);
    EXPECT_EQ(imported_result.Cell1DsTangents, expectedResult.Cell1DsTangents);
    EXPECT_EQ(imported_result.Cell1DsVertices, expectedResult.Cell1DsVertices);
}

TEST(TestMeshUtilities, TestAgglomerateCell1Ds)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestAgglomerateCell1Ds";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO meshDao(mesh_data);

    meshUtilities.FillMesh1D(geometryUtilities, Eigen::Vector3d(1.0, 1.0, 0.0), Eigen::Vector3d(1.0, 1.0, 0.0), {0.0, 0.25, 0.75, 1.0}, meshDao);
    meshUtilities.ComputeCell0DCell1DNeighbours(meshDao);

    std::vector<std::vector<unsigned int>> meshCell1DToConvexCell1DIndices(meshDao.Cell1DTotalNumber());
    for (unsigned int c1D_index = 0; c1D_index < meshDao.Cell1DTotalNumber(); c1D_index++)
    {
        meshCell1DToConvexCell1DIndices.at(c1D_index) = std::vector<unsigned int>({c1D_index});
    }

    meshUtilities.ExportMeshToVTU(meshDao, exportFolder, "OriginalMesh");

    {
        const auto agglomerationInfo = meshUtilities.AgglomerateCell1Ds(std::unordered_set<unsigned int>({1, 0, 2}), meshDao);

        ASSERT_EQ(std::vector<unsigned int>({0, 3}), agglomerationInfo.AgglomerateCell1DVertices);
        ASSERT_EQ(std::vector<unsigned int>({1, 2}), agglomerationInfo.SubCell1DsRemovedVertices);

        const unsigned int agglomeratedCell1DIndex =
            meshUtilities.AgglomerateCell1Ds(std::unordered_set<unsigned int>({1, 0, 2}),
                                             agglomerationInfo.AgglomerateCell1DVertices,
                                             agglomerationInfo.SubCell1DsRemovedVertices,
                                             meshDao,
                                             meshCell1DToConvexCell1DIndices);

        ASSERT_EQ(3, agglomeratedCell1DIndex);
        ASSERT_EQ(std::vector<unsigned int>({2, 0, 1}), meshCell1DToConvexCell1DIndices.at(3));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData activeMeshData;
    meshUtilities.ExtractActiveMesh(meshDao, activeMeshData);
    std::vector<std::vector<unsigned int>> activeMeshCell1DToConvexCell1DIndices(meshDao.Cell1DTotalNumber());
    for (unsigned int c1D_index = 0; c1D_index < meshDao.Cell1DTotalNumber(); c1D_index++)
    {
        activeMeshCell1DToConvexCell1DIndices.at(c1D_index) =
            meshCell1DToConvexCell1DIndices.at(activeMeshData.NewCell1DToOldCell1D.at(c1D_index));
    }

    meshUtilities.ExportMeshToVTU(meshDao, exportFolder, "AgglomeratedMesh");
}

TEST(TestMeshUtilities, Test_CollapseCell1D_Mesh2D_Polygon)
{
    std::string exportFolder = "./Export/TestMeshUtilities/Test_CollapseCell1D_Mesh2D_Polygon";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.MinTolerance = 1.0e-14;
    geometry_utilities_config.Tolerance1D = 1.0e-6;
    geometry_utilities_config.Tolerance2D = 1.0e-12;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d::Zero(),
                                                          1.0);

      mesh_utilities.CreatePolygonalMesh(geometry_utilities,
                                         square,
                                         10,
                                         10,
                                         mesh,
                                         10);
    }

    mesh_utilities.ComputeCell0DCell1DNeighbours(mesh);
    mesh_utilities.ComputeCell0DCell2DNeighbours(mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Original_Mesh");

    mesh_utilities.CollapseCell1D(20,
                                  mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_1");

    mesh_utilities.CollapseCell1D(18,
                                  mesh);

   mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_2");

   mesh_utilities.CollapseCell1D(33,
                                 mesh);

   mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_3");

   mesh_utilities.CollapseCell1D(37,
                                 mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Final_Mesh");

    Gedim::MeshUtilities::ExtractActiveMeshData extraction_data;
    mesh_utilities.ExtractActiveMesh(mesh, extraction_data);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Filtered_Mesh");

    {
        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        config.Cell2D_CheckConvexity = false;
        mesh_utilities.CheckMesh2D(config,
                                   geometry_utilities,
                                   mesh);
    }

    {
        const std::vector<Gedim::GeometryUtilities::PolygonTypes> cell2Ds_types(mesh.Cell2DTotalNumber(),
                                                                          Gedim::GeometryUtilities::PolygonTypes::Generic_Concave);
        const auto mesh_geometric_data = mesh_utilities.FillMesh2DGeometricData(geometry_utilities, mesh, cell2Ds_types);

        Gedim::MeshUtilities::CheckMeshGeometricData2DConfiguration config;
        mesh_utilities.CheckMeshGeometricData2D(config, geometry_utilities, mesh, mesh_geometric_data);
    }

}

TEST(TestMeshUtilities, Test_CollapseCell1D_Mesh3D_Polyhedron)
{
    std::string exportFolder = "./Export/TestMeshUtilities/Test_CollapseCell1D_Mesh3D_Polyhedron";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.MinTolerance = 1.0e-14;
    geometry_utilities_config.Tolerance1D = 1.0e-6;
    geometry_utilities_config.Tolerance2D = 1.0e-10;
    geometry_utilities_config.Tolerance3D = 1.0e-8;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);
    Gedim::MeshUtilities mesh_utilities;

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      const auto cube = geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d::Zero(),
                                                          1.0);

      mesh_utilities.CreatePolyhedralMesh(geometry_utilities,
                                         cube.Vertices,
                                          cube.Edges,
                                          cube.Faces,
                                         10,
                                         100,
                                         mesh,
                                         10);
    }

    mesh_utilities.ComputeCell0DCell1DNeighbours(mesh);
    mesh_utilities.ComputeCell0DCell2DNeighbours(mesh);
    mesh_utilities.ComputeCell0DCell3DNeighbours(mesh);
    mesh_utilities.ComputeCell1DCell2DNeighbours(mesh);
    mesh_utilities.ComputeCell1DCell3DNeighbours(mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Original_Mesh");

    mesh_utilities.CollapseCell1D(14,
                                  mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_1");

    mesh_utilities.CollapseCell1D(9,
                                  mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_2");

    mesh_utilities.CollapseCell1D(12,
                                  mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Mesh_3");

    mesh_utilities.CollapseCell1D(17,
                                  mesh);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Final_Mesh");

    Gedim::MeshUtilities::ExtractActiveMeshData extraction_data;
    mesh_utilities.ExtractActiveMesh(mesh, extraction_data);

    mesh_utilities.ExportMeshToVTU(mesh, exportFolder, "Filtered_Mesh");

    {
      Gedim::Output::CreateFolder(exportFolder + "/Filtered_Mesh");
      mesh_utilities.ExportMeshToCsv(mesh,
                                     ';',
                                     exportFolder + "/Filtered_Mesh");
    }

    {
        Gedim::MeshUtilities::CheckMesh3DConfiguration config;
        config.Cell2D_CheckConvexity = false;
        config.Cell3D_CheckConvexity = false;
        mesh_utilities.CheckMesh3D(config,
                                   geometry_utilities,
                                   mesh);
    }

    {
        mesh_utilities.ComputeCell2DCell3DNeighbours(mesh);
        const auto mesh_geometric_data = mesh_utilities.FillMesh3DGeometricData(geometry_utilities, mesh);

        Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration config;
        config.Cell3D_CheckTetrahedra = false;
        mesh_utilities.CheckMeshGeometricData3D(config, geometry_utilities, mesh, mesh_geometric_data);
    }

}


} // namespace GedimUnitTesting

#endif // __TEST_MESH_UTILITIES1D_H
