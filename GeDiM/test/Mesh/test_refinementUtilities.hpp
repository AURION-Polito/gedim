#ifndef __TEST_REFINEMENT_UTILITIES_H
#define __TEST_REFINEMENT_UTILITIES_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

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

  TEST(TestRefinementUtilities, TestRefineTriangles_Random)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefineTriangles_Random";
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
    const unsigned int maxRefinements = 10;

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData = meshUtilities.FillMesh2DGeometricData(geometryUtilities,
                                                                                                          meshDAO);
      std::vector<unsigned int> cell2DsToRefineIndex = Gedim::Utilities::SortArrayIndices(meshGeometricData.Cell2DsAreas);
      cell2DsToRefineIndex.resize(0.5 * meshDAO.Cell2DTotalNumber());


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
  }
}

#endif // __TEST_REFINEMENT_UTILITIES_H
