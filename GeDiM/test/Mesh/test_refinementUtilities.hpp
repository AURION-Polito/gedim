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

    refinementUtilities.RefineTriangleCellByMaxEdge(cell2DToRefineIndex,
                                                    meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                    meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                    meshDAO);

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
        refinementUtilities.RefineTriangleCellByMaxEdge(cell2DToRefineIndex,
                                                        meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                        meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                        meshDAO);
      }

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R" +
                                    to_string(r));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    EXPECT_EQ(21, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(48, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(28, meshDAO.Cell2DTotalNumber());

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");
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

    refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                       meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                       lineTangent,
                                                       lineOrigin,
                                                       cell1DsQualityParameter,
                                                       meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                       meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                       meshDAO);

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

    refinementUtilities.RefinePolygonalCellByDirection(cell2DToRefineIndex,
                                                       meshGeometricData.Cell2DsVertices[cell2DToRefineIndex],
                                                       lineTangent,
                                                       lineOrigin,
                                                       cell1DsQualityParameter,
                                                       meshGeometricData.Cell2DsEdgeLengths.at(cell2DToRefineIndex),
                                                       meshGeometricData.Cell2DsEdgeDirections.at(cell2DToRefineIndex),
                                                       meshDAO);

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
}

#endif // __TEST_REFINEMENT_UTILITIES_H
