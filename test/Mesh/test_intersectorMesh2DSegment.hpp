#ifndef __TEST_IntersectorMesh2DSegment_H
#define __TEST_IntersectorMesh2DSegment_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "IntersectorMesh2DSegment.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
  TEST(TestIntersectorMesh2DSegment, TestIntersectMesh)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1e-15;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // segment inside mesh2D
      {
        GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMesh;
        Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.5, 0.0);
        Vector3d segmentEnd(0.5, 0.25, 0.0);

        Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh,
                                                                 geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
        ASSERT_NO_THROW(IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin,
                                                                        segmentEnd,
                                                                        result));

        EXPECT_EQ(result.Points.size(), 2);
        EXPECT_EQ(result.Points[0.0].Cell2DIds.size(), 1);
        EXPECT_EQ(result.Points[1.0].Cell2DIds.size(), 1);
        EXPECT_EQ(result.Points[1.0].Cell2DIds, result.Points[1.0].Cell2DIds);
        EXPECT_EQ(result.Segments.size(), 1);
        EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 0);
        EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 1);
      }

      // segmentalong single edge
      {
        GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMesh;
        Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.75, 0.0);
        Vector3d segmentEnd(   0.75, 0.25, 0.0);

        Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh,
                                                                 geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
        ASSERT_NO_THROW(IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin,
                                                                        segmentEnd,
                                                                        result));

        EXPECT_EQ(result.Points.size(), 2);
        EXPECT_EQ(result.Points[0.0].Vertex2DIds.size(), 0);
        EXPECT_EQ(result.Points[0.0].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Points[0.0].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Points[1.0].Vertex2DIds.size(), 0);
        EXPECT_EQ(result.Points[1.0].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Points[1.0].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Points[0.0].Edge2DIds, result.Points[1.0].Edge2DIds);
        EXPECT_EQ(result.Points[0.0].Cell2DIds, result.Points[1.0].Cell2DIds);
        EXPECT_EQ(result.Segments.size(), 1);
        EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 2);
      }

      // segmentalong two edges
      {
        GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMesh;
        Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.75, 0.0);
        Vector3d segmentEnd(   0.75, 0.25, 0.0);

        Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh,
                                                                 geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
        ASSERT_NO_THROW(IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin,
                                                                        segmentEnd,
                                                                        result));

        EXPECT_EQ(result.Points.size(), 3);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Points[4.9999999999999994e-01].Vertex2DIds.size(), 1);
        EXPECT_EQ(result.Points[4.9999999999999994e-01].Edge2DIds.size(), 4);
        EXPECT_EQ(result.Points[4.9999999999999994e-01].Cell2DIds.size(), 4);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Segments.size(), 2);
        EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Segments[1].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Segments[1].Cell2DIds.size(), 2);
      }

      // segmentalong generic mesh
      {
        GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
        Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

        Vector3d segmentOrigin(0.25, 0.75, 0.0);
        Vector3d segmentEnd(   0.75, 0.25, 0.0);

        Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh,
                                                                 geometryUtilities);

        Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
        ASSERT_NO_THROW(IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin,
                                                                        segmentEnd,
                                                                        result));

        EXPECT_EQ(result.Points.size(), 4);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Vertex2DIds.size(), 1);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Edge2DIds.size(), 7);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Cell2DIds.size(), 7);
        EXPECT_EQ(result.Points[2.5000000000000006e-01].Vertex2DIds.size(), 0);
        EXPECT_EQ(result.Points[2.5000000000000006e-01].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Points[2.5000000000000006e-01].Cell2DIds.size(), 2);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Vertex2DIds.size(), 1);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Edge2DIds.size(), 7);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Cell2DIds.size(), 7);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Vertex2DIds.size(), 1);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Edge2DIds.size(), 8);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Cell2DIds.size(), 8);
        EXPECT_EQ(result.Segments.size(), 3);
        EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 0);
        EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 1);
        EXPECT_EQ(result.Segments[1].Edge2DIds.size(), 0);
        EXPECT_EQ(result.Segments[1].Cell2DIds.size(), 1);
        EXPECT_EQ(result.Segments[2].Edge2DIds.size(), 1);
        EXPECT_EQ(result.Segments[2].Cell2DIds.size(), 2);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }


}

#endif // __TEST_IntersectorMesh2DSegment_H
