#ifndef __TEST_REFINEMENT_UTILITIES_3D_H
#define __TEST_REFINEMENT_UTILITIES_3D_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "MeshMatrices_3D_22Cells_Mock.hpp"

#include "ConformMeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"
#include "RefinementUtilities.hpp"
#include "CommonUtilities.hpp"
#include "VTKUtilities.hpp"

#include <array>

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestRefinementUtilities, TestRefineTetrahedrons)
  {
    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefineTetrahedron";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);
    MeshMatrices_3D_22Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);
    meshUtilities.ComputeCell1DCell3DNeighbours(meshDAO);
    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Original");

    Gedim::MeshUtilities::MeshGeometricData3D meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                        meshDAO);

    const unsigned int cell3DToRefineIndex = 0;

    const Eigen::MatrixXd& cell3DVertices = meshGeometricData.Cell3DsVertices.at(cell3DToRefineIndex);
    const Eigen::MatrixXi& cell3DEdges = meshGeometricData.Cell3DsEdges.at(cell3DToRefineIndex);
    const std::vector<Eigen::MatrixXi>& cell3DFaces = meshGeometricData.Cell3DsFaces.at(cell3DToRefineIndex);

    const Gedim::RefinementUtilities::TetrahedronMaxEdgeDirection direction = refinementUtilities.ComputeTetrahedronMaxEdgeDirection(cell3DEdges,
                                                                                                                                     meshGeometricData.Cell3DsEdgeLengths.at(cell3DToRefineIndex));
    EXPECT_EQ(2, direction.MaxEdgeIndex);
    EXPECT_EQ(2, direction.OppositeVerticesIndex[0]);
    EXPECT_EQ(3, direction.OppositeVerticesIndex[1]);

    const Eigen::Vector3d planeOrigin = 0.5 * (cell3DVertices.col(cell3DEdges(0, direction.MaxEdgeIndex)) +
                                               cell3DVertices.col(cell3DEdges(1, direction.MaxEdgeIndex)));

    Eigen::Matrix3d planeTriangle;
    planeTriangle.col(0)<< planeOrigin;
    planeTriangle.col(1)<< cell3DVertices.col(direction.OppositeVerticesIndex[1]);
    planeTriangle.col(2)<< cell3DVertices.col(direction.OppositeVerticesIndex[0]);

    const Eigen::Vector3d planeNormal = geometryUtilities.PolygonNormal(planeTriangle);
    const Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
    const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

    const Gedim::RefinementUtilities::RefinePolyhedron_Result result =
        refinementUtilities.RefinePolyhedronCell_ByPlane(cell3DToRefineIndex,
                                                         cell3DVertices,
                                                         cell3DEdges,
                                                         meshGeometricData.Cell3DsEdgeLengths.at(cell3DToRefineIndex),
                                                         cell3DFaces,
                                                         meshGeometricData.Cell3DsFaces3DVertices.at(cell3DToRefineIndex),
                                                         meshGeometricData.Cell3DsFacesEdge3DTangents.at(cell3DToRefineIndex),
                                                         meshGeometricData.Cell3DsFacesTranslations.at(cell3DToRefineIndex),
                                                         meshGeometricData.Cell3DsFacesRotationMatrices.at(cell3DToRefineIndex),
                                                         meshGeometricData.Cell3DsVolumes.at(cell3DToRefineIndex),
                                                         planeNormal,
                                                         planeOrigin,
                                                         planeRotationMatrix,
                                                         planeTranslation,
                                                         meshDAO);
    ASSERT_EQ(std::vector<unsigned int>({ 20 }),
              result.NewCell0DsIndex);
    ASSERT_EQ(3,
              result.NewCell1DsIndex.size());
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::Updated,
              result.NewCell1DsIndex[0].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 59, 60 }),
              result.NewCell1DsIndex[0].NewCell1DsIndex);
    ASSERT_EQ(2,
              result.NewCell1DsIndex[0].OriginalCell1DIndex);
    ASSERT_EQ(20,
              result.NewCell1DsIndex[0].NewCell0DIndex);
    ASSERT_EQ(2,
              result.NewCell1DsIndex[0].OriginalCell3DEdgeIndex);
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::New,
              result.NewCell1DsIndex[1].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 61 }),
              result.NewCell1DsIndex[1].NewCell1DsIndex);
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::New,
              result.NewCell1DsIndex[2].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 62 }),
              result.NewCell1DsIndex[2].NewCell1DsIndex);
    ASSERT_EQ(3,
              result.NewCell2DsIndex.size());
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::Updated,
              result.NewCell2DsIndex[0].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 62, 63 }),
              result.NewCell2DsIndex[0].NewCell2DsIndex);
    ASSERT_EQ(2,
              result.NewCell2DsIndex[0].OriginalCell2DIndex);
    ASSERT_EQ(61,
              result.NewCell2DsIndex[0].NewCell1DIndex);
    ASSERT_EQ(3,
              result.NewCell2DsIndex[0].OriginalCell3DFaceIndex);
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::Updated,
              result.NewCell2DsIndex[1].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 64, 65 }),
              result.NewCell2DsIndex[1].NewCell2DsIndex);
    ASSERT_EQ(3,
              result.NewCell2DsIndex[1].OriginalCell2DIndex);
    ASSERT_EQ(62,
              result.NewCell2DsIndex[1].NewCell1DIndex);
    ASSERT_EQ(0,
              result.NewCell2DsIndex[1].OriginalCell3DFaceIndex);
    ASSERT_EQ(Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::New,
              result.NewCell2DsIndex[2].Type);
    ASSERT_EQ(std::vector<unsigned int>({ 66 }),
              result.NewCell2DsIndex[2].NewCell2DsIndex);
    ASSERT_EQ(std::vector<unsigned int>({ 22, 23 }),
              result.NewCell3DsIndex);

    unsigned int step = 0;

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined_Step" + std::to_string(step++));

    std::map<unsigned int, unsigned int> updatedNeighCell2Ds;

    for (unsigned int f = 0; f < result.NewCell2DsIndex.size(); f++)
    {
      if (result.NewCell2DsIndex[f].Type != Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::Updated)
        continue;

      const unsigned int cell2DIndex = result.NewCell2DsIndex[f].OriginalCell2DIndex;

      std::list<unsigned int> splitCell1DsOriginalIndex;
      std::list<unsigned int> splitCell1DsNewCell0DIndex;
      std::list<std::vector<unsigned int>> splitCell1DsUpdatedIndices;

      for (unsigned int e = 0; e < result.NewCell2DsIndex[f].NewCell1DsPosition.size(); e++)
      {
        const Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D& refinedCell1D =
            result.NewCell1DsIndex.at(result.NewCell2DsIndex[f].NewCell1DsPosition[e]);

        if (refinedCell1D.Type !=
            Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::Updated)
          continue;

        splitCell1DsOriginalIndex.push_back(refinedCell1D.OriginalCell1DIndex);
        splitCell1DsNewCell0DIndex.push_back(refinedCell1D.NewCell0DIndex);
        splitCell1DsUpdatedIndices.push_back(refinedCell1D.NewCell1DsIndex);
      }

      const Gedim::RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result updateResult = refinementUtilities.RefinePolyhedronCell_UpdateFaceNeighbours(cell3DToRefineIndex,
                                                                                                                                                             cell2DIndex,
                                                                                                                                                             result.NewCell2DsIndex[f].NewCell1DIndex,
                                                                                                                                                             std::vector<unsigned int>(splitCell1DsOriginalIndex.begin(),
                                                                                                                                                                                       splitCell1DsOriginalIndex.end()),
                                                                                                                                                             std::vector<unsigned int>(splitCell1DsNewCell0DIndex.begin(),
                                                                                                                                                                                       splitCell1DsNewCell0DIndex.end()),
                                                                                                                                                             std::vector<std::vector<unsigned int>>(splitCell1DsUpdatedIndices.begin(),
                                                                                                                                                                                                    splitCell1DsUpdatedIndices.end()),
                                                                                                                                                             result.NewCell2DsIndex[f].NewCell2DsIndex,
                                                                                                                                                             meshGeometricData.Cell3DsFacesEdgeDirections,
                                                                                                                                                             updatedNeighCell2Ds,
                                                                                                                                                             meshDAO);
      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_Refined_Step" + std::to_string(step++));
    }

    for (unsigned int e = 0; e < result.NewCell1DsIndex.size(); e++)
    {
      if (result.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::Updated)
        continue;

      const unsigned int cell1DIndex = result.NewCell1DsIndex[e].OriginalCell1DIndex;

      const Gedim::RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result updateResult = refinementUtilities.RefinePolyhedronCell_UpdateEdgeNeighbours(cell3DToRefineIndex,
                                                                                                                                                             cell1DIndex,
                                                                                                                                                             result.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                                                                                                             result.NewCell1DsIndex[e].NewCell0DIndex,
                                                                                                                                                             meshGeometricData.Cell3DsFacesEdgeDirections,
                                                                                                                                                             updatedNeighCell2Ds,
                                                                                                                                                             meshDAO);
      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_Refined_Step" + std::to_string(step++));
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    ASSERT_EQ(21, meshDAO.Cell0DTotalNumber());
    ASSERT_EQ(62, meshDAO.Cell1DTotalNumber());
    ASSERT_EQ(65, meshDAO.Cell2DTotalNumber());
    ASSERT_EQ(23, meshDAO.Cell3DTotalNumber());

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_Refined");
  }

  TEST(TestRefinementUtilities, TestRefineTetrahedrons_ByVolume)
  {
    GTEST_SKIP();

    std::string exportFolder = "./Export/TestRefinementUtilities/TestRefineTetrahedrons_ByVolume";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::RefinementUtilities refinementUtilities(geometryUtilities,
                                                   meshUtilities);

    MeshMatrices_3D_22Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO meshDAO(mockMesh.Mesh);
    meshUtilities.ComputeCell1DCell3DNeighbours(meshDAO);
    meshUtilities.ComputeCell2DCell3DNeighbours(meshDAO);

    unsigned int step = 0;
    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_R" +
                                  to_string(step++),
                                  true);

    const unsigned int seed = 10;
    const unsigned int maxRefinements = 5;

    std::vector<unsigned int> cell2DsAligned(meshDAO.Cell2DTotalNumber());
    std::iota(std::begin(cell2DsAligned),
              std::end(cell2DsAligned),
              1);
    unsigned int maxCell2DAligned = meshDAO.Cell2DTotalNumber();

    for (unsigned int r = 0; r < maxRefinements; r++)
    {
      Gedim::MeshUtilities::MeshGeometricData3D meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                          meshDAO);

      const std::vector<bool>& activeCell3Ds = mockMesh.Mesh.ActiveCell3D;
      std::vector<double> activeCell3DsVolume(meshDAO.Cell3DTotalNumber(), 0.0);
      unsigned int numActiveCell3Ds = 0;
      for (unsigned int c = 0; c < meshDAO.Cell3DTotalNumber(); c++)
      {
        if (!activeCell3Ds[c])
          continue;

        activeCell3DsVolume[c] = meshGeometricData.Cell3DsVolumes[c];
        numActiveCell3Ds++;
      }

      std::vector<unsigned int> cell3DsToRefineIndex = Gedim::Utilities::SortArrayIndices(activeCell3Ds);
      std::reverse(cell3DsToRefineIndex.begin(), cell3DsToRefineIndex.end());
      cell3DsToRefineIndex.resize(0.5 * numActiveCell3Ds);

      for (unsigned int c = 0; c < cell3DsToRefineIndex.size(); c++)
      {
        meshGeometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                  meshDAO);

        std::list<unsigned int> updatedCell3Ds;
        meshDAO.Cell3DUpdatedCell3Ds(cell3DsToRefineIndex[c],
                                     updatedCell3Ds);

        if (updatedCell3Ds.size() > 1)
          continue;

        const unsigned int cell3DToRefineIndex = updatedCell3Ds.size() == 0 ?
                                                   cell3DsToRefineIndex[c] :
                                                   updatedCell3Ds.front();

        const Eigen::MatrixXd& cell3DVertices = meshGeometricData.Cell3DsVertices.at(cell3DToRefineIndex);
        const Eigen::MatrixXi& cell3DEdges = meshGeometricData.Cell3DsEdges.at(cell3DToRefineIndex);
        const std::vector<Eigen::MatrixXi>& cell3DFaces = meshGeometricData.Cell3DsFaces.at(cell3DToRefineIndex);

        const std::vector<std::vector<unsigned int>> facesUnalignedPoints = geometryUtilities.PolyhedronFacesUnalignedVertices(meshGeometricData.Cell3DsFaces2DVertices.at(cell3DToRefineIndex));

        std::map<unsigned int, std::list<unsigned int>> polyhedronUnaligedFacesList;
        for (unsigned int f = 0; f < cell3DFaces.size(); f++)
        {
          const unsigned int cell2DIndex = meshDAO.Cell3DFace(cell3DToRefineIndex, f);
          const unsigned int cell2DAligned = cell2DsAligned.at(cell2DIndex);

          if (polyhedronUnaligedFacesList.find(cell2DAligned) ==
              polyhedronUnaligedFacesList.end())
          {
            polyhedronUnaligedFacesList.insert(std::make_pair(cell2DAligned,
                                                              std::list<unsigned int>()));
          }

          polyhedronUnaligedFacesList[cell2DAligned].push_back(f);
        }

        std::vector<std::vector<unsigned int>> polyhedronUnaligedFaces(polyhedronUnaligedFacesList.size());
        unsigned int uf = 0;
        for (const auto& polyhedronUnaligedFaceList : polyhedronUnaligedFacesList)
          polyhedronUnaligedFaces[uf++] = std::vector<unsigned int>(polyhedronUnaligedFaceList.second.begin(),
                                                                    polyhedronUnaligedFaceList.second.end());

        const std::vector<unsigned int> unalignedVertices = geometryUtilities.UnalignedPolyhedronPoints(cell3DVertices,
                                                                                                        cell3DFaces,
                                                                                                        meshGeometricData.Cell3DsFacesTranslations.at(cell3DToRefineIndex),
                                                                                                        meshGeometricData.Cell3DsFacesRotationMatrices.at(cell3DToRefineIndex),
                                                                                                        polyhedronUnaligedFaces,
                                                                                                        facesUnalignedPoints);

        // create unaligned tetrahedron
        const Gedim::GeometryUtilities::Polyhedron unalignedTetrahedron = geometryUtilities.CreateTetrahedronWithVertices(cell3DVertices.col(unalignedVertices.at(0)),
                                                                                                                          cell3DVertices.col(unalignedVertices.at(1)),
                                                                                                                          cell3DVertices.col(unalignedVertices.at(2)),
                                                                                                                          cell3DVertices.col(unalignedVertices.at(3)));

        const Eigen::VectorXd unalignedEdgesLength = geometryUtilities.PolyhedronEdgesLength(unalignedTetrahedron.Vertices,
                                                                                             unalignedTetrahedron.Edges);

        const Gedim::RefinementUtilities::TetrahedronMaxEdgeDirection direction = refinementUtilities.ComputeTetrahedronMaxEdgeDirection(unalignedTetrahedron.Edges,
                                                                                                                                         unalignedEdgesLength);

        const Eigen::Vector3d planeOrigin = 0.5 * (unalignedTetrahedron.Vertices.col(unalignedTetrahedron.Edges(0, direction.MaxEdgeIndex)) +
                                                   unalignedTetrahedron.Vertices.col(unalignedTetrahedron.Edges(1, direction.MaxEdgeIndex)));

        Eigen::Matrix3d planeTriangle;
        planeTriangle.col(0)<< planeOrigin;
        planeTriangle.col(1)<< unalignedTetrahedron.Vertices.col(direction.OppositeVerticesIndex[1]);
        planeTriangle.col(2)<< unalignedTetrahedron.Vertices.col(direction.OppositeVerticesIndex[0]);

        {
          Gedim::VTKUtilities exporter;
          exporter.AddPolyhedron(unalignedTetrahedron.Vertices,
                                 unalignedTetrahedron.Edges,
                                 unalignedTetrahedron.Faces);
          exporter.Export(exportFolder + "/CellToSPlit.vtu");
        }

        {
          Gedim::VTKUtilities exporter;
          exporter.AddPolygon(planeTriangle);
          exporter.Export(exportFolder + "/Plane.vtu");
        }

        const Eigen::Vector3d planeNormal = geometryUtilities.PolygonNormal(planeTriangle);
        const Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal);
        const Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);


        const Gedim::RefinementUtilities::RefinePolyhedron_Result result =
            refinementUtilities.RefinePolyhedronCell_ByPlane(cell3DToRefineIndex,
                                                             cell3DVertices,
                                                             cell3DEdges,
                                                             meshGeometricData.Cell3DsEdgeLengths.at(cell3DToRefineIndex),
                                                             cell3DFaces,
                                                             meshGeometricData.Cell3DsFaces3DVertices.at(cell3DToRefineIndex),
                                                             meshGeometricData.Cell3DsFacesEdge3DTangents.at(cell3DToRefineIndex),
                                                             meshGeometricData.Cell3DsFacesTranslations.at(cell3DToRefineIndex),
                                                             meshGeometricData.Cell3DsFacesRotationMatrices.at(cell3DToRefineIndex),
                                                             meshGeometricData.Cell3DsVolumes.at(cell3DToRefineIndex),
                                                             planeNormal,
                                                             planeOrigin,
                                                             planeRotationMatrix,
                                                             planeTranslation,
                                                             meshDAO);

        cell2DsAligned.resize(meshDAO.Cell2DTotalNumber());

        for (unsigned int f = 0; f < result.NewCell2DsIndex.size(); f++)
        {
          switch (result.NewCell2DsIndex[f].Type)
          {
            case Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::Updated:
            {
              const unsigned int cell2DIndex = result.NewCell2DsIndex[f].OriginalCell2DIndex;
              for (const unsigned int newCell2DIndex : result.NewCell2DsIndex[f].NewCell2DsIndex)
                cell2DsAligned[newCell2DIndex] = cell2DsAligned.at(cell2DIndex);
            }
              break;
            case Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::New:
            {
              const unsigned int newCell2DIndex = result.NewCell2DsIndex[f].NewCell2DsIndex[0];
              cell2DsAligned[newCell2DIndex] = maxCell2DAligned++;
            }
              break;
            default:
              throw std::runtime_error("Unknown face type");
          }
        }

        std::map<unsigned int, unsigned int> updatedNeighCell2Ds;
        for (unsigned int f = 0; f < result.NewCell2DsIndex.size(); f++)
        {
          if (result.NewCell2DsIndex[f].Type != Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell2D::Types::Updated)
            continue;

          const unsigned int cell2DIndex = result.NewCell2DsIndex[f].OriginalCell2DIndex;

          std::list<unsigned int> splitCell1DsOriginalIndex;
          std::list<unsigned int> splitCell1DsNewCell0DIndex;
          std::list<std::vector<unsigned int>> splitCell1DsUpdatedIndices;

          for (unsigned int e = 0; e < result.NewCell2DsIndex[f].NewCell1DsPosition.size(); e++)
          {
            const Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D& refinedCell1D =
                result.NewCell1DsIndex.at(result.NewCell2DsIndex[f].NewCell1DsPosition[e]);

            if (refinedCell1D.Type !=
                Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::Updated)
              continue;

            splitCell1DsOriginalIndex.push_back(refinedCell1D.OriginalCell1DIndex);
            splitCell1DsNewCell0DIndex.push_back(refinedCell1D.NewCell0DIndex);
            splitCell1DsUpdatedIndices.push_back(refinedCell1D.NewCell1DsIndex);
          }

          const Gedim::RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result updateResult = refinementUtilities.RefinePolyhedronCell_UpdateFaceNeighbours(cell3DToRefineIndex,
                                                                                                                                                                 cell2DIndex,
                                                                                                                                                                 result.NewCell2DsIndex[f].NewCell1DIndex,
                                                                                                                                                                 std::vector<unsigned int>(splitCell1DsOriginalIndex.begin(),
                                                                                                                                                                                           splitCell1DsOriginalIndex.end()),
                                                                                                                                                                 std::vector<unsigned int>(splitCell1DsNewCell0DIndex.begin(),
                                                                                                                                                                                           splitCell1DsNewCell0DIndex.end()),
                                                                                                                                                                 std::vector<std::vector<unsigned int>>(splitCell1DsUpdatedIndices.begin(),
                                                                                                                                                                                                        splitCell1DsUpdatedIndices.end()),
                                                                                                                                                                 result.NewCell2DsIndex[f].NewCell2DsIndex,
                                                                                                                                                                 meshGeometricData.Cell3DsFacesEdgeDirections,
                                                                                                                                                                 updatedNeighCell2Ds,
                                                                                                                                                                 meshDAO);

        }

        for (unsigned int e = 0; e < result.NewCell1DsIndex.size(); e++)
        {
          if (result.NewCell1DsIndex[e].Type != Gedim::RefinementUtilities::RefinePolyhedron_Result::RefinedCell1D::Types::Updated)
            continue;

          const unsigned int cell1DIndex = result.NewCell1DsIndex[e].OriginalCell1DIndex;

          const Gedim::RefinementUtilities::RefinePolyhedron_UpdateNeighbour_Result updateResult = refinementUtilities.RefinePolyhedronCell_UpdateEdgeNeighbours(cell3DToRefineIndex,
                                                                                                                                                                 cell1DIndex,
                                                                                                                                                                 result.NewCell1DsIndex[e].NewCell1DsIndex,
                                                                                                                                                                 result.NewCell1DsIndex[e].NewCell0DIndex,
                                                                                                                                                                 meshGeometricData.Cell3DsFacesEdgeDirections,
                                                                                                                                                                 updatedNeighCell2Ds,
                                                                                                                                                                 meshDAO);
        }

        cell2DsAligned.resize(meshDAO.Cell2DTotalNumber());

        for (const auto& updatedNeighCell2D : updatedNeighCell2Ds)
          cell2DsAligned[updatedNeighCell2D.second] = cell2DsAligned.at(updatedNeighCell2D.first);
      }

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "Mesh_R" +
                                    to_string(step++),
                                    true);
    }

    Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
    meshUtilities.ExtractActiveMesh(meshDAO,
                                    extractionData);

    meshUtilities.ExportMeshToVTU(meshDAO,
                                  exportFolder,
                                  "Mesh_R" +
                                  to_string(step++),
                                  true);

    Gedim::MeshUtilities::CheckMesh3DConfiguration checkConfig;
    meshUtilities.CheckMesh3D(checkConfig,
                              geometryUtilities,
                              meshDAO);

    EXPECT_EQ(92, meshDAO.Cell0DTotalNumber());
    EXPECT_EQ(348, meshDAO.Cell1DTotalNumber());
    EXPECT_EQ(420, meshDAO.Cell2DTotalNumber());
    EXPECT_EQ(163, meshDAO.Cell3DTotalNumber());
  }
}

#endif // __TEST_REFINEMENT_UTILITIES_3D_H
