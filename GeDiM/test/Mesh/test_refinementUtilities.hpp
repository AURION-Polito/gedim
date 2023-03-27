#ifndef __TEST_REFINEMENT_UTILITIES_H
#define __TEST_REFINEMENT_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"

#include "ConformMeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "RefinementUtilities.hpp"
#include "CommonUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestRefinementUtilities, TestRefineTriangles)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefineTriangles";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_2Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        meshDAO);

    const unsigned int cell2DToRefineIndex = 0;

    const Gedim::RefinementUtilities::MaxEdgeDirection direction = refinementUtilities.ComputeTriangleMaxEdgeDirection(meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex));
    EXPECT_EQ(2, direction.MaxEdgeIndex);
    EXPECT_EQ(1, direction.OppositeVertexIndex);

    const Gedim::RefinementUtilities::RefinePolygon_Result result = refinementUtilities.RefineTriangleCell_ByEdge(cell2DToRefineIndex,
                                                                                                                  direction.MaxEdgeIndex,
                                                                                                                  direction.OppositeVertexIndex,
                                                                                                                  meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                  meshGeometricData.Cell2DsAreas.at(cell2DToRefineIndex),
                                                                                                                  meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                  meshDAO);
    EXPECT_EQ(std::vector<unsigned int>({ 4 }),
              result.NewCell0DsIndex);
    EXPECT_EQ(2,
              result.NewCell1DsIndex.size());
    EXPECT_EQ(Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated,
              result.NewCell1DsIndex[0].Type);
    EXPECT_EQ(std::vector<unsigned int>({ 5, 6 }),
              result.NewCell1DsIndex[0].NewCell1DsIndex);
    EXPECT_EQ(2,
              result.NewCell1DsIndex[0].OriginalCell1DIndex);
    EXPECT_EQ(4,
              result.NewCell1DsIndex[0].NewCell0DIndex);
    EXPECT_EQ(Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::New,
              result.NewCell1DsIndex[1].Type);
    EXPECT_EQ(std::vector<unsigned int>({ 7 }),
              result.NewCell1DsIndex[1].NewCell1DsIndex);
    EXPECT_EQ(std::vector<unsigned int>({ 2, 3 }),
              result.NewCell2DsIndex);

    for (unsigned int e = 0; e < result.NewCell1DsIndex.size(); e++)
    {
      if (result.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefineTriangleCell_UpdateNeighbours(cell2DToRefineIndex,
                                                              result.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                              result.NewCell1DsIndex[e].NewCell0DIndex,
                                                              result.NewCell1DsIndex[e].NewCell1DsIndex,
                                                              meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(result.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                              meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    EXPECT_EQ(5, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(8, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(4, meshDAO.Cell2DTotalNumber());

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");
  }

  TEST(TestRefinementUtilities, TestRefineTriangles_ByArea)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefineTriangles_ByArea";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_2Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    const unsigned int seed = 10;
    const unsigned int maxRefinements = 5;

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);
      const std::vector<bool>& activeCell2Ds = mockMesh.Mesh.ActiveCell2D;
      std::vector<double> activeCell2DsArea(meshDAO.Cell2DTotalNumber(), 0.0);
      unsigned int numActiveCell2Ds = 0;
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      {
        if (!activeCell2Ds[c])
          continue;

        activeCell2DsArea[c] = meshGeometricData.Cell2DsAreas[c];
        numActiveCell2Ds++;
      }

      std::vector<unsigned int> cell2DsToRefineIndex = Gedim::Utilities::SortArrayIndices(activeCell2DsArea);
      std::reverse(cell2DsToRefineIndex.begin(), cell2DsToRefineIndex.end());
      cell2DsToRefineIndex.resize(0.5 * numActiveCell2Ds);

      for (unsigned int c = 0; c < cell2DsToRefineIndex.size(); c++)
      {
        const unsigned int cell2DToRefineIndex = cell2DsToRefineIndex[c];

        const Gedim::RefinementUtilities::MaxEdgeDirection direction = refinementUtilities.ComputeTriangleMaxEdgeDirection(meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex));

        const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefineTriangleCell_ByEdge(cell2DToRefineIndex,
                                                                                                                            direction.MaxEdgeIndex,
                                                                                                                            direction.OppositeVertexIndex,
                                                                                                                            meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsAreas.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

        for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
        {
          if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
            continue;

          refinementUtilities.RefineTriangleCell_UpdateNeighbours(cell2DToRefineIndex,
                                                                  refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                                  refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                                  refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                  meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                                  meshDAO);
        }
      }

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R" +
                                    to_string(r));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
    meshUtilities.CheckMesh2D(checkConfig,
                              geometryUtilities,
                              meshDAO);

    EXPECT_EQ(21, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(48, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(28, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_NoNewVertices)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_NoNewVertices";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_1Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        meshDAO);
    const std::vector<double> cell2DsQualityParameter = { meshGeometricData.Cell2DsEdgeLengths[0].minCoeff() };
    const double cell1DsQualityWeight = 1.0;
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), cell2DsQualityParameter[0]);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(1.0, 1.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.25, 0.25, 0.0);
    const unsigned int cell2DToRefineIndex = 0;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            cell1DsQualityParameter,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(4, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(5, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(2, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_NewVertexOne)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_NewVertexOne";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_2Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        meshDAO);
    const std::vector<double> cell2DsQualityParameter = { meshGeometricData.Cell2DsEdgeLengths[0].minCoeff(),
                                                          meshGeometricData.Cell2DsEdgeLengths[1].minCoeff() };
    const double cell1DsQualityWeight = 1.0;
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(-1.0, 0.5, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(1.0, 0.5, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            cell1DsQualityParameter,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(5, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(7, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_NewVertexTwo)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_NewVertexTwo";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_2Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        meshDAO);
    const std::vector<double> cell2DsQualityParameter = { meshGeometricData.Cell2DsEdgeLengths[0].minCoeff(),
                                                          meshGeometricData.Cell2DsEdgeLengths[1].minCoeff() };
    const double cell1DsQualityWeight = 1.0;
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(1.0, 1.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.25, 0.25, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            cell1DsQualityParameter,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(5, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(7, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_NewVertices)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_NewVertices";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_2Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                        meshDAO);
    const std::vector<double> cell2DsQualityParameter = { meshGeometricData.Cell2DsEdgeLengths[0].minCoeff(),
                                                          meshGeometricData.Cell2DsEdgeLengths[1].minCoeff() };
    const double cell1DsQualityWeight = 1.0;
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(0.5, 0.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.5, 0.5, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            cell1DsQualityParameter,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(6, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(8, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_CheckQuality_NoNewVertices)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_CheckQuality_NoNewVertices";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    const Eigen::Vector3d rectangleOrigin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBaseTangent = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeightTangent = Eigen::Vector3d(0.0, 1.0, 0.0);
    const vector<double> baseCoordinates = { 0.0, 1.0 };
    const vector<double> heightCoordinates = { 0.0, 0.45, 0.55, 1.0 };

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseCoordinates,
                                      heightCoordinates,
                                      meshDAO);
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

    const double cell1DsQualityWeight = 0.5;
    const Eigen::Vector3d lineTangent = Eigen::Vector3d(1.0, 0.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.5, 0.5, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2Ds.Vertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            meshGeometricData.Cell1Ds.Quality,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(8, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(10, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_CheckQuality_NewVertexOne)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_CheckQuality_NewVertexOne";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    const Eigen::Vector3d rectangleOrigin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBaseTangent = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeightTangent = Eigen::Vector3d(0.0, 1.0, 0.0);
    const vector<double> baseCoordinates = { 0.0, 1.0 };
    const vector<double> heightCoordinates = { 0.0, 0.45, 0.55, 1.0 };

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseCoordinates,
                                      heightCoordinates,
                                      meshDAO);
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

    const double cell1DsQualityWeight = 0.5;
    const Eigen::Vector3d lineTangent = Eigen::Vector3d(-0.5, -0.05, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.5, 0.55, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2Ds.Vertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            meshGeometricData.Cell1Ds.Quality,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(8, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(10, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_CheckQuality_NewVertexTwo)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_CheckQuality_NewVertexTwo";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);

    Gedim::MeshMatrices mesh;
    Gedim::MeshMatricesDAO meshDAO(mesh);

    const Eigen::Vector3d rectangleOrigin = Eigen::Vector3d(0.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleBaseTangent = Eigen::Vector3d(1.0, 0.0, 0.0);
    const Eigen::Vector3d rectangleHeightTangent = Eigen::Vector3d(0.0, 1.0, 0.0);
    const vector<double> baseCoordinates = { 0.0, 1.0 };
    const vector<double> heightCoordinates = { 0.0, 0.45, 0.55, 1.0 };

    meshUtilities.CreateRectangleMesh(rectangleOrigin,
                                      rectangleBaseTangent,
                                      rectangleHeightTangent,
                                      baseCoordinates,
                                      heightCoordinates,
                                      meshDAO);
    meshUtilities.ComputeCell1DCell2DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

    const double cell1DsQualityWeight = 0.5;
    const Eigen::Vector3d lineTangent = Eigen::Vector3d(-0.5, 0.04, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(1.0, 0.51, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                            meshGeometricData.Cell2Ds.Vertices[cell2DToRefineIndex],
                                                                                                                            lineTangent,
                                                                                                                            lineOrigin,
                                                                                                                            meshGeometricData.Cell1Ds.Quality,
                                                                                                                            cell1DsQualityWeight,
                                                                                                                            meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                            meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                            meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                             meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                             meshDAO);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    EXPECT_EQ(8, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(10, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(3, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_ByArea_MaxDiameter)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_ByArea_MaxDiameter";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_1Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    const unsigned int seed = 10;
    const unsigned int maxRefinements = 6;
    const double cell1DsQualityWeight = 0.5;

    Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      const std::vector<bool>& activeCell2Ds = mockMesh.Mesh.ActiveCell2D;

      std::vector<double> activeCell2DsArea(meshDAO.Cell2DTotalNumber(), 0.0);
      unsigned int numActiveCell2Ds = 0;
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      {
        if (!activeCell2Ds[c])
          continue;

        activeCell2DsArea[c] = meshGeometricData.Cell2Ds.Area[c];
        numActiveCell2Ds++;
      }

      std::vector<unsigned int> cell2DsToRefineIndex = Gedim::Utilities::SortArrayIndices(activeCell2DsArea);
      std::reverse(cell2DsToRefineIndex.begin(), cell2DsToRefineIndex.end());
      const unsigned int numCell2DsToRefine = numActiveCell2Ds == 1 ? 1 : 0.5 * numActiveCell2Ds;
      cell2DsToRefineIndex.resize(numCell2DsToRefine);

      for (unsigned int c = 0; c < cell2DsToRefineIndex.size(); c++)
      {
        std::list<unsigned int> cell2DsToUpdateGeometricData;

        std::list<unsigned int> updatedCell2Ds;
        meshDAO.Cell2DUpdatedCell2Ds(cell2DsToRefineIndex[c],
                                     updatedCell2Ds);

        if (updatedCell2Ds.size() > 1)
          continue;

        const unsigned int cell2DToRefineIndex = updatedCell2Ds.size() == 0 ?
                                                   cell2DsToRefineIndex[c] :
                                                   updatedCell2Ds.front();

        Gedim::RefinementUtilities::PolygonDirection direction = refinementUtilities.ComputePolygonMaxDiameterDirection(meshGeometricData.Cell2Ds.Vertices.at(cell2DToRefineIndex),
                                                                                                                        meshGeometricData.Cell2Ds.Centroid.at(cell2DToRefineIndex));

        if (r == 0)
        {
          direction.LineOrigin = Eigen::Vector3d(0.0, 0.25, 0.0);
          direction.LineTangent = Eigen::Vector3d(1.0, 0.5, 0.0);
        }

        const Gedim::RefinementUtilities::RefinePolygon_Result refineResult  = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                                 meshGeometricData.Cell2Ds.Vertices[cell2DToRefineIndex],
                                                                                                                                 direction.LineTangent,
                                                                                                                                 direction.LineOrigin,
                                                                                                                                 meshGeometricData.Cell1Ds.Quality,
                                                                                                                                 cell1DsQualityWeight,
                                                                                                                                 meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                                 meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                                 meshDAO);

        for (unsigned int rnc = 0; rnc < refineResult.NewCell2DsIndex.size(); rnc++)
          cell2DsToUpdateGeometricData.push_back(refineResult.NewCell2DsIndex[rnc]);

        for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
        {
          if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
            continue;

          const std::vector<unsigned int> newNeighboursCell2DsIndex = refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                                                                             meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                                                                                             meshDAO);
          for (unsigned int rnc = 0; rnc < newNeighboursCell2DsIndex.size(); rnc++)
            cell2DsToUpdateGeometricData.push_back(newNeighboursCell2DsIndex[rnc]);
        }

        refinementUtilities.RefinePolygonCell_UpdateGeometricData(meshDAO,
                                                                  std::vector<unsigned int>(cell2DsToUpdateGeometricData.begin(),
                                                                                            cell2DsToUpdateGeometricData.end()),
                                                                  meshGeometricData);
      }

      Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
      meshUtilities.CheckMesh2D(checkConfig,
                                geometryUtilities,
                                meshDAO);

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R" +
                                    to_string(r));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
    meshUtilities.CheckMesh2D(checkConfig,
                              geometryUtilities,
                              meshDAO);

    EXPECT_EQ(12, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(24, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(13, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefinePolygons_ByArea_MaxInertia)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_ByArea_MaxInertia";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_2D_1Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    const unsigned int seed = 10;
    const unsigned int maxRefinements = 6;
    const double cell1DsQualityWeight = 0.5;

    Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      const std::vector<bool>& activeCell2Ds = mockMesh.Mesh.ActiveCell2D;

      std::vector<double> activeCell2DsArea(meshDAO.Cell2DTotalNumber(), 0.0);
      unsigned int numActiveCell2Ds = 0;
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      {
        if (!activeCell2Ds[c])
          continue;

        activeCell2DsArea[c] = meshGeometricData.Cell2Ds.Area[c];
        numActiveCell2Ds++;
      }

      std::vector<unsigned int> cell2DsToRefineIndex = Gedim::Utilities::SortArrayIndices(activeCell2DsArea);
      std::reverse(cell2DsToRefineIndex.begin(), cell2DsToRefineIndex.end());
      const unsigned int numCell2DsToRefine = numActiveCell2Ds == 1 ? 1 : 0.5 * numActiveCell2Ds;
      cell2DsToRefineIndex.resize(numCell2DsToRefine);

      for (unsigned int c = 0; c < cell2DsToRefineIndex.size(); c++)
      {
        std::list<unsigned int> cell2DsToUpdateGeometricData;

        std::list<unsigned int> updatedCell2Ds;
        meshDAO.Cell2DUpdatedCell2Ds(cell2DsToRefineIndex[c],
                                     updatedCell2Ds);

        if (updatedCell2Ds.size() > 1)
          continue;

        const unsigned int cell2DToRefineIndex = updatedCell2Ds.size() == 0 ?
                                                   cell2DsToRefineIndex[c] :
                                                   updatedCell2Ds.front();

        Gedim::RefinementUtilities::PolygonDirection direction = refinementUtilities.ComputePolygonMaxInertiaDirection(meshGeometricData.Cell2Ds.Centroid.at(cell2DToRefineIndex),
                                                                                                                       meshGeometricData.Cell2Ds.Inertia[cell2DToRefineIndex]);

        if (r == 0)
        {
          direction.LineOrigin = Eigen::Vector3d(0.0, 0.25, 0.0);
          direction.LineTangent = Eigen::Vector3d(1.0, 0.5, 0.0);
        }

        const Gedim::RefinementUtilities::RefinePolygon_Result refineResult  = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                                 meshGeometricData.Cell2Ds.Vertices[cell2DToRefineIndex],
                                                                                                                                 direction.LineTangent,
                                                                                                                                 direction.LineOrigin,
                                                                                                                                 meshGeometricData.Cell1Ds.Quality,
                                                                                                                                 cell1DsQualityWeight,
                                                                                                                                 meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                                 meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                                 meshDAO);

        for (unsigned int rnc = 0; rnc < refineResult.NewCell2DsIndex.size(); rnc++)
          cell2DsToUpdateGeometricData.push_back(refineResult.NewCell2DsIndex[rnc]);

        for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
        {
          if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
            continue;

          const std::vector<unsigned int> newNeighboursCell2DsIndex = refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                                                                                             refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                                                                             meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                                                                                             meshDAO);
          for (unsigned int rnc = 0; rnc < newNeighboursCell2DsIndex.size(); rnc++)
            cell2DsToUpdateGeometricData.push_back(newNeighboursCell2DsIndex[rnc]);
        }

        refinementUtilities.RefinePolygonCell_UpdateGeometricData(meshDAO,
                                                                  std::vector<unsigned int>(cell2DsToUpdateGeometricData.begin(),
                                                                                            cell2DsToUpdateGeometricData.end()),
                                                                  meshGeometricData);
      }

      Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
      meshUtilities.CheckMesh2D(checkConfig,
                                geometryUtilities,
                                meshDAO);

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R" +
                                    to_string(r));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");

    Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
    meshUtilities.CheckMesh2D(checkConfig,
                              geometryUtilities,
                              meshDAO);

    EXPECT_EQ(12, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(24, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(13, meshDAO.Cell2DTotalNumber());
  }

  TEST(TestRefinementUtilities, TestRefineMeshTwoSegments_ByArea_MaxInertia)
  {
    try
    {
      string exportFolder = "./Export";
      Gedim::Output::CreateFolder(exportFolder);
      exportFolder = "./Export/TestRefinementUtilities/TestRefineMeshTwoSegments_ByArea_MaxInertia";
      Gedim::Output::CreateFolder(exportFolder);

      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-8;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      Gedim::MeshUtilities meshUtilities;

      Gedim::ConformMeshUtilities conformMeshUtilities(geometryUtilities,
                                                       meshUtilities);
      Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                     meshUtilities);

      MeshMatrices_2D_1Cells_Mock mockMesh;
      Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);

      const unsigned int numSegments = 2;

      std::vector<Eigen::MatrixXd> segmentsVertices(numSegments, Eigen::MatrixXd(3, 2));
      std::vector<Eigen::Vector3d> segmentsTangent(numSegments);
      std::vector<Eigen::Vector3d> segmentsBarycenter(numSegments);
      std::vector<double> segmentsLength(numSegments);
      std::vector<double> segmentsSquaredLength(numSegments);

      segmentsVertices[0].col(0)<< 0.1, 0.1, 0.0;
      segmentsVertices[0].col(1)<< 0.9, 0.7, 0.0;
      segmentsVertices[1].col(0)<< 0.75, 0.0, 0.0;
      segmentsVertices[1].col(1)<< 0.4, 0.8, 0.0;

      // compute segment geometric properties
      for (unsigned int s = 0; s < numSegments; s++)
      {
        segmentsTangent[s] = geometryUtilities.SegmentTangent(segmentsVertices.at(s).col(0),
                                                              segmentsVertices.at(s).col(1));
        segmentsBarycenter[s] = geometryUtilities.SegmentBarycenter(segmentsVertices.at(s).col(0),
                                                                    segmentsVertices.at(s).col(1));
        segmentsLength[s] = geometryUtilities.SegmentLength(segmentsVertices.at(s).col(0),
                                                            segmentsVertices.at(s).col(1));
        segmentsSquaredLength[s] = segmentsLength[s] *
                                   segmentsLength[s];
      }

      // compute segment intersections
      const vector<list<double>> segmentsAdditionalPoints = geometryUtilities.IntersectionsBetweenSegments(segmentsVertices,
                                                                                                           segmentsTangent,
                                                                                                           segmentsBarycenter,
                                                                                                           segmentsLength);

      // compute conformed mesh
      std::vector<Gedim::IntersectorMesh2DSegment::IntersectionMesh> segmentsIntersectionMesh(numSegments);
      std::vector<std::vector<double>> segmentsCurvilinearCoordinatesMesh(numSegments);
      std::vector<Gedim::UnionMeshSegment::UnionMesh> segmentsUnionMesh(numSegments);
      std::vector<Gedim::ConformerMeshSegment::ConformMesh> segmentsConformMesh(numSegments);
      Gedim::ConformMeshUtilities::ComputeDomainConformedMeshOptions options;
      options.PrintStatus = false;
      options.VtkExportFolder = exportFolder;

      conformMeshUtilities.ComputeConformedMeshWithSegments(segmentsAdditionalPoints,
                                                            segmentsVertices,
                                                            segmentsTangent,
                                                            segmentsBarycenter,
                                                            segmentsLength,
                                                            segmentsSquaredLength,
                                                            meshDAO,
                                                            segmentsIntersectionMesh,
                                                            segmentsCurvilinearCoordinatesMesh,
                                                            segmentsUnionMesh,
                                                            segmentsConformMesh,
                                                            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::Generalized,
                                                            options);

      segmentsUnionMesh.clear();
      segmentsIntersectionMesh.clear();
      segmentsCurvilinearCoordinatesMesh.clear();
      segmentsConformMesh.clear();

      segmentsUnionMesh.resize(numSegments);
      segmentsIntersectionMesh.resize(numSegments);
      segmentsCurvilinearCoordinatesMesh.resize(numSegments);
      segmentsConformMesh.resize(numSegments);

      conformMeshUtilities.ComputeConformedMeshWithSegments(vector<list<double>>(numSegments),
                                                            segmentsVertices,
                                                            segmentsTangent,
                                                            segmentsBarycenter,
                                                            segmentsLength,
                                                            segmentsSquaredLength,
                                                            meshDAO,
                                                            segmentsIntersectionMesh,
                                                            segmentsCurvilinearCoordinatesMesh,
                                                            segmentsUnionMesh,
                                                            segmentsConformMesh,
                                                            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges,
                                                            options);

      // refine mesh
      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_Original");

      const unsigned int seed = 10;
      const unsigned int maxRefinements = 6;
      const double cell1DsQualityWeight = 0.5;

      Gedim::RefinementUtilities::Cell2Ds_GeometricData meshGeometricData = refinementUtilities.RefinePolygonCell_InitializeGeometricData(meshDAO);

      for (unsigned int r = 0; r < maxRefinements; r++)
      {
        std::vector<double> activeCell2DsArea(meshDAO.Cell2DTotalNumber(), 0.0);
        unsigned int numActiveCell2Ds = 0;
        for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
        {
          if (!meshDAO.Cell2DIsActive(c))
            continue;

          activeCell2DsArea[c] = meshGeometricData.Cell2Ds.Area[c];
          numActiveCell2Ds++;
        }

        std::vector<unsigned int> cell2DsToRefineIndex = Gedim::Utilities::SortArrayIndices(activeCell2DsArea);
        std::reverse(cell2DsToRefineIndex.begin(), cell2DsToRefineIndex.end());
        const unsigned int numCell2DsToRefine = numActiveCell2Ds == 1 ? 1 : 0.5 * numActiveCell2Ds;
        cell2DsToRefineIndex.resize(numCell2DsToRefine);

        for (unsigned int c = 0; c < cell2DsToRefineIndex.size(); c++)
        {
          std::list<unsigned int> cell2DsToUpdateGeometricData;

          std::list<unsigned int> updatedCell2Ds;
          meshDAO.Cell2DUpdatedCell2Ds(cell2DsToRefineIndex[c],
                                       updatedCell2Ds);

          if (updatedCell2Ds.size() > 1)
            continue;

          const unsigned int cell2DToRefineIndex = updatedCell2Ds.size() == 0 ?
                                                     cell2DsToRefineIndex[c] :
                                                     updatedCell2Ds.front();

          Gedim::RefinementUtilities::PolygonDirection direction = refinementUtilities.ComputePolygonMaxInertiaDirection(meshGeometricData.Cell2Ds.Centroid.at(cell2DToRefineIndex),
                                                                                                                         meshGeometricData.Cell2Ds.Inertia.at(cell2DToRefineIndex));

          const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonCell_ByDirection(cell2DToRefineIndex,
                                                                                                                                  meshGeometricData.Cell2Ds.Vertices.at(cell2DToRefineIndex),
                                                                                                                                  direction.LineTangent,
                                                                                                                                  direction.LineOrigin,
                                                                                                                                  meshGeometricData.Cell1Ds.Quality,
                                                                                                                                  cell1DsQualityWeight,
                                                                                                                                  meshGeometricData.Cell2Ds.EdgesLength.at(cell2DToRefineIndex),
                                                                                                                                  meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex),
                                                                                                                                  meshDAO);
          for (unsigned int rnc = 0; rnc < refineResult.NewCell2DsIndex.size(); rnc++)
            cell2DsToUpdateGeometricData.push_back(refineResult.NewCell2DsIndex[rnc]);

          for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
          {
            if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
              continue;

            const std::vector<unsigned int> newNeighboursCell2DsIndex = refinementUtilities.RefinePolygonCell_UpdateNeighbours(cell2DToRefineIndex,
                                                                                                                               refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                                                                                               refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                                                                                               refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                                                                               meshGeometricData.Cell2Ds.EdgesDirection.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                                                                                               meshDAO);
            for (unsigned int rnc = 0; rnc < newNeighboursCell2DsIndex.size(); rnc++)
              cell2DsToUpdateGeometricData.push_back(newNeighboursCell2DsIndex[rnc]);
          }

          refinementUtilities.RefinePolygonCell_UpdateGeometricData(meshDAO,
                                                                    std::vector<unsigned int>(cell2DsToUpdateGeometricData.begin(),
                                                                                              cell2DsToUpdateGeometricData.end()),
                                                                    meshGeometricData);
        }

        Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
        meshUtilities.CheckMesh2D(checkConfig,
                                  geometryUtilities,
                                  meshDAO);

        meshUtilities.ExportMeshToVTU(meshDAO,
                                      exportFolder,
                                      "Mesh_R" +
                                      to_string(r));
      }

      Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
      meshUtilities.ExtractActiveMesh(meshDAO,
                                      extractionData);

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_Refined");

      Gedim::MeshUtilities::CheckMesh2DConfiguration checkConfig;
      meshUtilities.CheckMesh2D(checkConfig,
                                geometryUtilities,
                                meshDAO);

      EXPECT_EQ(66, meshDAO.Cell0DTotalNumber());
      EXPECT_EQ(107, meshDAO.Cell1DTotalNumber());
      EXPECT_EQ(42, meshDAO.Cell2DTotalNumber());
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_REFINEMENT_UTILITIES_H
