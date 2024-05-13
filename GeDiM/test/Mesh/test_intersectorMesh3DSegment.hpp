#ifndef __TEST_IntersectorMesh3DSegment_H
#define __TEST_IntersectorMesh3DSegment_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_3D_1Cells_Mock.hpp"
#include "IntersectorMesh3DSegment.hpp"
#include "VTKUtilities.hpp"

using namespace testing;

namespace GedimUnitTesting
{
  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentInside)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    // segment inside mesh3D
    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh3D(mockMesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData3D mesh3D_geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                                 mesh3D);

    const Eigen::Vector3d segmentOrigin(0.25, 0.5, 0.0);
    const Eigen::Vector3d segmentEnd(0.5, 0.25, 0.0);
    const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                            segmentEnd);
    const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                  segmentEnd);
    const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                 segmentEnd);

    const Gedim::IntersectorMesh3DSegment intersectorMesh3DSegment(geometryUtilities,
                                                                   meshUtilities);

    const Gedim::IntersectorMesh3DSegment::IntersectionMesh result =
        intersectorMesh3DSegment.CreateIntersectionMesh(segmentOrigin,
                                                        segmentEnd,
                                                        segmentTangent,
                                                        segmentBarycenter,
                                                        segmentLength,
                                                        mesh3D,
                                                        mesh3D_geometricData);



    EXPECT_EQ(result.Points.size(), 2);
    EXPECT_EQ(result.Points[0].CurvilinearCoordinate, 0.0);
    EXPECT_EQ(result.Points[1].CurvilinearCoordinate, 1.0);
    EXPECT_EQ(result.Points[0].Cell3DIds.size(), 1);
    EXPECT_EQ(result.Points[1].Cell3DIds.size(), 1);
    EXPECT_EQ(result.Points[0].Cell3DIds[0], 0);
    EXPECT_EQ(result.Points[1].Cell3DIds[0], 0);
    EXPECT_EQ(result.Segments.size(), 1);
    EXPECT_EQ(result.Segments[0].Cell3DIds.size(), 1);
    EXPECT_EQ(result.Segments[0].Cell3DIds[0], 0);
  }
}

#endif // __TEST_IntersectorMesh3DSegment_H
