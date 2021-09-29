#ifndef __TEST_ExportConformedMesh_H
#define __TEST_ExportConformedMesh_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"
#include "ConformerMeshSegment.hpp"
#include "ConformerMeshPolygon.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestExportConformedMesh, TestExportConformedMesh2D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // conform of simple 2 points mesh with one segment
      {
        GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockOriginalMesh;
        Gedim::MeshMatricesDAO domainOriginalMesh(mockOriginalMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.25, 0.0);
        Vector3d segmentEnd(0.35, 0.35, 0.0);

        Gedim::IntersectorMesh2DSegment intersectorMeshSegment(domainOriginalMesh,
                                                               geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMesh;
        intersectorMeshSegment.CreateIntersectionMesh(segmentOrigin,
                                                      segmentEnd,
                                                      intersectionMesh);

        // Create mesh union
        vector<double> curvilinearCoordinatesMesh;
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMesh,
                                                                  curvilinearCoordinatesMesh);
        Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

        Gedim::UnionMeshSegment::UnionMesh unionMesh;
        unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMesh,
                                         curvilinearCoordinatesMesh,
                                         unionMesh);

        Gedim::ConformerMeshSegment conformMeshSegment(geometryUtilities);

        Gedim::ConformerMeshSegment::ConformMesh conformMesh;
        ASSERT_NO_THROW(conformMeshSegment.CreateConformMesh(intersectionMesh,
                                                             unionMesh,
                                                             0,
                                                             conformMesh));

        Gedim::ConformerMeshPolygon conformMeshPolygon(geometryUtilities);
        Gedim::ConformerMeshPolygon::ConformMesh domainConformedMesh;

        ASSERT_NO_THROW(conformMeshPolygon.CreateConformMesh(segmentOrigin,
                                                             segmentEnd,
                                                             conformMesh,
                                                             domainOriginalMesh,
                                                             domainConformedMesh));

        Gedim::MeshUtilities meshUtilities;

        Gedim::MeshUtilities::ExtractActiveMeshData extractionData;
        ASSERT_NO_THROW(meshUtilities.ExtractActiveMesh(domainOriginalMesh,
                                                        extractionData));

        {
          using namespace Gedim;
          cerr<< extractionData.NewCell0DToOldCell0D<< endl;
          cerr<< extractionData.NewCell1DToOldCell1D<< endl;
          cerr<< extractionData.NewCell2DToOldCell2D<< endl;
          cerr<< extractionData.NewCell3DToOldCell3D<< endl;
          cerr<< extractionData.OldCell0DToNewCell0D<< endl;
          cerr<< extractionData.OldCell1DToNewCell1D<< endl;
          cerr<< extractionData.OldCell2DToNewCell2D<< endl;
          cerr<< extractionData.OldCell3DToNewCell3D<< endl;
        }

        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell0D, 7);
        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell1D, 9);
        EXPECT_EQ(mockOriginalMesh.Mesh.NumberCell2D, 3);

        Gedim::ConformerMeshSegment::ToString(conformMesh);
        for (const auto& mesh1Dpoint : conformMesh.Points)
          EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
        for (const auto& mesh1Dsegment : conformMesh.Segments)
          EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 2);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_ExportConformedMesh_H
