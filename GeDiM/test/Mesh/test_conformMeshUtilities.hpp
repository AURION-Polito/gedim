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

#include "ConformMeshUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestConformMeshUtilities, TestConformMeshSingleSegment)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1e-14;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      Gedim::MeshUtilities meshUtilities;

      GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
      Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

      Gedim::ConformMeshUtilities conformMeshUtilities(geometryUtilities,
                                                       meshUtilities);

      const unsigned int numSegments = 2;
      std::vector<Eigen::MatrixXd> segmentsVertices(numSegments);
      std::vector<Eigen::Vector3d> segmentsTangent(numSegments);
      std::vector<Eigen::Vector3d> segmentsBarycenter(numSegments);
      std::vector<double> segmentsLength(numSegments);
      std::vector<double> segmentsSquaredLength(numSegments);
      vector<list<double>> segmentsAdditionalPoints(numSegments);

      std::vector<Gedim::IntersectorMesh2DSegment::IntersectionMesh>& segmentsIntersectionMesh,
      std::vector<std::vector<double> >& segmentsCurvilinearCoordinatesMesh,
      vector<UnionMeshSegment::UnionMesh>& segmentsUnionMesh,
      vector<ConformerMeshSegment::ConformMesh>& segmentsConformMesh,
      vector<ConformerMeshPolygon::ConformMesh>& segmentsConformMeshInfo,
      const ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types& conformDomainMeshType,
      const ComputeDomainConformedMeshOptions& options

      conformMeshUtilities.ComputeDomainConformedMesh(segmentsAdditionalPoints,
                                                      segmentsVertices,
                                                      segmentsTangent,
                                                      segmentsBarycenter,
                                                      segmentsLength,
                                                      segmentsSquaredLength,
                                                      meshDAO,
                                                      );

    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_ConformMeshUtilities_H
