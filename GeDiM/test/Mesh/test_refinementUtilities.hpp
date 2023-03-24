#ifndef __TEST_REFINEMENT_UTILITIES_H
#define __TEST_REFINEMENT_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_2D_1Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"

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

    const Gedim::RefinementUtilities::RefinePolygon_Result result = refinementUtilities.RefineTriangleCellByEdge(cell2DToRefineIndex,
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

      refinementUtilities.RefineTriangleCellByEdge_UpdateNeighbours(cell2DToRefineIndex,
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

        const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefineTriangleCellByEdge(cell2DToRefineIndex,
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

          refinementUtilities.RefineTriangleCellByEdge_UpdateNeighbours(cell2DToRefineIndex,
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
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), cell2DsQualityParameter[0]);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(1.0, 1.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.25, 0.25, 0.0);
    const unsigned int cell2DToRefineIndex = 0;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                                                                                             meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                             lineTangent,
                                                                                                                             lineOrigin,
                                                                                                                             cell1DsQualityParameter,
                                                                                                                             meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                             meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonalCellByDirection_UpdateNeighbours(cell2DToRefineIndex,
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
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(-1.0, 0.5, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(1.0, 0.5, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                                                                                             meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                             lineTangent,
                                                                                                                             lineOrigin,
                                                                                                                             cell1DsQualityParameter,
                                                                                                                             meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                             meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonalCellByDirection_UpdateNeighbours(cell2DToRefineIndex,
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
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(1.0, 1.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.25, 0.25, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                                                                                             meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                             lineTangent,
                                                                                                                             lineOrigin,
                                                                                                                             cell1DsQualityParameter,
                                                                                                                             meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                             meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonalCellByDirection_UpdateNeighbours(cell2DToRefineIndex,
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
    const std::vector<double> cell1DsQualityParameter(meshDAO.Cell1DTotalNumber(), 0.0);

    const Eigen::Vector3d lineTangent = Eigen::Vector3d(0.5, 0.0, 0.0).normalized();
    const Eigen::Vector3d lineOrigin = Eigen::Vector3d(0.5, 0.5, 0.0);
    const unsigned int cell2DToRefineIndex = 1;

    const Gedim::RefinementUtilities::RefinePolygon_Result refineResult = refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                                                                                             meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                             lineTangent,
                                                                                                                             lineOrigin,
                                                                                                                             cell1DsQualityParameter,
                                                                                                                             meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                             meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                             meshDAO);

    for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
    {
      if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
        continue;

      refinementUtilities.RefinePolygonalCellByDirection_UpdateNeighbours(cell2DToRefineIndex,
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

  TEST(TestRefinementUtilities, TestRefinePolygons_ByArea)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefinePolygons_ByArea";
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
    const unsigned int maxRefinements = 5;

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);


      const std::vector<bool>& activeCell2Ds = mockMesh.Mesh.ActiveCell2D;

      std::vector<double> cell2DsInRadius(meshDAO.Cell2DTotalNumber(), 0.0);
      for (unsigned int c = 0; c < meshDAO.Cell2DTotalNumber(); c++)
      {
        if (!activeCell2Ds[c])
          continue;

        cell2DsInRadius[c] = geometryUtilities.PolygonInRadius(meshGeometricData.Cell2DsVertices[c],
                                                               meshGeometricData.Cell2DsCentroids[c],
                                                               meshGeometricData.Cell2DsEdgeNormals[c]);
      }

      Gedim::RefinementUtilities::MeshQuality meshQuality = refinementUtilities.ComputeMeshQualityForRefinement(meshDAO,
                                                                                                                meshGeometricData.Cell2DsEdgeLengths,
                                                                                                                cell2DsInRadius);

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
      const unsigned int numCell2DsToRefine = numActiveCell2Ds == 1 ? 1 : 0.5 * numActiveCell2Ds;
      cell2DsToRefineIndex.resize(numCell2DsToRefine);

      for (unsigned int c = 0; c < cell2DsToRefineIndex.size(); c++)
      {
        const unsigned int cell2DToRefineIndex = cell2DsToRefineIndex[c];

        const Gedim::RefinementUtilities::PolygonDirection direction = refinementUtilities.ComputePolygonMaxDiameterDirection(meshGeometricData.Cell2DsVertices.at(cell2DToRefineIndex),
                                                                                                                              meshGeometricData.Cell2DsCentroids.at(cell2DToRefineIndex));

        const Gedim::RefinementUtilities::RefinePolygon_Result refineResult  = refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                                                                                                  meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                                                                                                  direction.LineTangent,
                                                                                                                                  direction.LineOrigin,
                                                                                                                                  meshQuality.Cell1DsQuality,
                                                                                                                                  meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                                                                                                  meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                                                                                                  meshDAO);

        for (unsigned int e = 0; e < refineResult.NewCell1DsIndex.size(); e++)
        {
          if (refineResult.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolygon_Result::RefinedCell1D::Types::Updated)
            continue;

          refinementUtilities.RefinePolygonalCellByDirection_UpdateNeighbours(cell2DToRefineIndex,
                                                                              refineResult.NewCell1DsIndex[e].OriginalCell1DIndex,
                                                                              refineResult.NewCell1DsIndex[e].NewCell0DIndex,
                                                                              refineResult.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                              meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex).at(refineResult.NewCell1DsIndex[e].OriginalCell2DEdgeIndex),
                                                                              meshDAO);
        }
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

    EXPECT_EQ(10, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(18, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(9, meshDAO.Cell2DTotalNumber());
  }
}

#endif // __TEST_REFINEMENT_UTILITIES_H
