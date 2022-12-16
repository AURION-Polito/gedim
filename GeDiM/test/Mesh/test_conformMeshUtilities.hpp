#ifndef __TEST_ConformMeshUtilities_H
#define __TEST_ConformMeshUtilities_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"

#include "MeshMatrices.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatricesDAO.hpp"

#include "VTKUtilities.hpp"

#include "ConformMeshUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestConformMeshUtilities, TestConformMeshTwoSegments)
  {
    try
    {
      string exportFolder = "./Export";
      Gedim::Output::CreateFolder(exportFolder);
      exportFolder = "./Export/TestConformMeshUtilities/TestConformMeshTwoSegments";
      Gedim::Output::CreateFolder(exportFolder);

      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1e-14;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      Gedim::MeshUtilities meshUtilities;

      GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "ConformedMesh");

      Gedim::ConformMeshUtilities conformMeshUtilities(geometryUtilities,
                                                       meshUtilities);

      const unsigned int numSegments = 2;

      std::vector<Eigen::MatrixXd> segmentsVertices(numSegments, Eigen::MatrixXd(3, 2));
      std::vector<Eigen::Vector3d> segmentsTangent(numSegments);
      std::vector<Eigen::Vector3d> segmentsBarycenter(numSegments);
      std::vector<double> segmentsLength(numSegments);
      std::vector<double> segmentsSquaredLength(numSegments);
      vector<list<double>> segmentsAdditionalPoints(numSegments);

      segmentsVertices[0].col(0)<< 0.1, 0.1, 0.0;
      segmentsVertices[0].col(1)<< 0.9, 0.7, 0.0;
      segmentsVertices[1].col(0)<< 0.75, 0.0, 0.0;
      segmentsVertices[1].col(1)<< 0.4, 0.8, 0.0;

      // export segments
      {
        Gedim::VTKUtilities exporter;
        for (unsigned int s = 0; s < numSegments; s++)
          exporter.AddSegment(segmentsVertices[s]);
        exporter.Export(exportFolder +
                        "/Segments.vtu");
      }

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
      for (unsigned int s1 = 0; s1 < numSegments; s1++)
      {
        const Eigen::MatrixXd& segmentOne = segmentsVertices[s1];
        const Eigen::Vector3d& segmentOneBarycenter = segmentsBarycenter[s1];
        const double& segmentOneLength = segmentsLength[s1];

        for (unsigned int s2 = s1 + 1; s2 < numSegments; s2++)
        {
          const Eigen::MatrixXd& segmentTwo = segmentsVertices[s2];
          const Eigen::Vector3d& segmentTwoBarycenter = segmentsBarycenter[s2];
          const double& segmentTwoLength = segmentsLength[s2];

          if (geometryUtilities.CheckNoSpheresIntersection(segmentOneBarycenter,
                                                           segmentTwoBarycenter,
                                                           segmentOneLength,
                                                           segmentTwoLength))
            continue;

          Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
              geometryUtilities.IntersectionSegmentSegment(segmentOne.col(0),
                                                           segmentOne.col(1),
                                                           segmentTwo.col(0),
                                                           segmentTwo.col(1));
          if (result.IntersectionSegmentsType !=
              Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection)
            continue;

          segmentsAdditionalPoints[s1].push_back(result.FirstSegmentIntersections[0].CurvilinearCoordinate);
          segmentsAdditionalPoints[s2].push_back(result.SecondSegmentIntersections[0].CurvilinearCoordinate);
        }
      }

      std::vector<Gedim::IntersectorMesh2DSegment::IntersectionMesh> segmentsIntersectionMesh(numSegments);
      std::vector<std::vector<double>> segmentsCurvilinearCoordinatesMesh(numSegments);
      std::vector<Gedim::UnionMeshSegment::UnionMesh> segmentsUnionMesh(numSegments);
      std::vector<Gedim::ConformerMeshSegment::ConformMesh> segmentsConformMesh(numSegments);
      std::vector<Gedim::ConformerMeshPolygon::ConformMesh> segmentsConformMeshInfo(numSegments);
      Gedim::ConformMeshUtilities::ComputeDomainConformedMeshOptions options;
      options.PrintStatus = true;
      options.VtkExportFolder = exportFolder;

      conformMeshUtilities.ComputeDomainConformedMesh(segmentsAdditionalPoints,
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
                                                      segmentsConformMeshInfo,
                                                      Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::Generalized,
                                                      options);

      segmentsUnionMesh.clear();
      segmentsIntersectionMesh.clear();
      segmentsCurvilinearCoordinatesMesh.clear();
      segmentsConformMesh.clear();
      segmentsConformMeshInfo.clear();

      segmentsUnionMesh.resize(numSegments);
      segmentsIntersectionMesh.resize(numSegments);
      segmentsCurvilinearCoordinatesMesh.resize(numSegments);
      segmentsConformMesh.resize(numSegments);
      segmentsConformMeshInfo.resize(numSegments);

      conformMeshUtilities.ComputeDomainConformedMesh(segmentsAdditionalPoints,
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
                                                      segmentsConformMeshInfo,
                                                      Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges,
                                                      options);

      meshUtilities.ExportMeshToVTU(meshDAO,
                                    exportFolder,
                                    "ConformedMesh");
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_ConformMeshUtilities_H
