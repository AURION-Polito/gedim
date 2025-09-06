// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __TEST_GEOMETRY_INTERSECTION_H
#define __TEST_GEOMETRY_INTERSECTION_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestGeometryUtilities, TestIntersectionSegmentSegment)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check simply no intersection
        {
            const Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.0);
            const Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.0);
            const Eigen::Vector3d segmentTwoOrigin(0.0, 0.0, 2.0);
            const Eigen::Vector3d segmentTwoEnd(1.0, 0.0, 2.0);

            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);
            ASSERT_TRUE(geometryUtilities.CheckNoSpheresIntersection(segmentOneBarycenter, segmentTwoBarycenter, segmentOneLength, segmentTwoLength));

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
        }

        // check no coplanarity
        {
            const Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.0);
            const Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.0);
            const Eigen::Vector3d segmentTwoOrigin(0.0, 0.0, 0.25);
            const Eigen::Vector3d segmentTwoEnd(0.0, 1.0, 0.25);

            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);
            ASSERT_TRUE(!geometryUtilities.CheckNoSpheresIntersection(segmentOneBarycenter, segmentTwoBarycenter, segmentOneLength, segmentTwoLength));

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
        }

        // check MultipleIntersections, no inclusion
        {
            const Eigen::Vector3d segmentOneOrigin(0.25, 0.75, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.25, 0.0);
            const Eigen::Vector3d segmentTwoOrigin(0.5, 0.5, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.0, 1.0, 0.0);

            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);
            ASSERT_FALSE(geometryUtilities.CheckNoSpheresIntersection(segmentOneBarycenter, segmentTwoBarycenter, segmentOneLength, segmentTwoLength));

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);

            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
            ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
            ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 0.5) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 0.5) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        }

        // check MultipleIntersections, no inclusion
        {
            Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25);
            Eigen::Vector3d segmentOneEnd(0.75, 0.0, 0.25);
            Eigen::Vector3d segmentTwoOrigin(2.0, 0.0, 0.25);
            Eigen::Vector3d segmentTwoEnd(0.5, 0.0, 0.25);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);

            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
            ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
            ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 2.0 / 3.0) <
                        geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 5.0 / 6.0) <
                        geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        }

        // check MultipleIntersections, with inclusion
        {
            Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25);
            Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.25);
            Eigen::Vector3d segmentTwoOrigin(0.5, 0.0, 0.25);
            Eigen::Vector3d segmentTwoEnd(0.75, 0.0, 0.25);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);

            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
            ASSERT_EQ(result.SecondIntersectionRelation[0], 0);
            ASSERT_EQ(result.SecondIntersectionRelation[1], 1);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 0.75) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        }

        // check OneIntersection
        {
            Eigen::Vector3d segmentOneOrigin(-1.0, 0.0, 0.25);
            Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.25);
            Eigen::Vector3d segmentTwoOrigin(0.0, -1.0, 0.25);
            Eigen::Vector3d segmentTwoEnd(0.0, 1.0, 0.25);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        }

        // check OneIntersection origin
        {
            Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25);
            Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.25);
            Eigen::Vector3d segmentTwoOrigin(0.0, -1.0, 0.25);
            Eigen::Vector3d segmentTwoEnd(0.0, 0.0, 0.25);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        }

        // check parallel segments
        {
            Eigen::Vector3d segmentOneOrigin(5.0000000000000056e-01, 1.1458333333333328e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(4.9999999999999994e-01, 6.2500000000000022e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(9.9999999999999978e-01, 1.9999999999999996e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoEnd(9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
        }

        // check parallel segments
        {
            Eigen::Vector3d segmentOneOrigin(1.1249999999999996e+00, 2.0000000000000000e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(1.1249999999999996e+00, 1.8602389358491287e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
            Eigen::Vector3d segmentTwoEnd(1.0000000000000000e+00, 2.2204460492503131e-16, 2.4825341532472739e-17);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
        }

        // check intersection in first segment end and second segment origin
        {
            Eigen::Vector3d segmentOneOrigin(1.2499999999999996e+00, 1.7105836549996734e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(9.9999999999999956e-01, 2.0000000000000000e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
            Eigen::Vector3d segmentTwoEnd(1.0000000000000000e+00, 2.2204460492503131e-16, 2.4825341532472739e-17);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
        }

        // check intersection inside segments
        {
            Eigen::Vector3d segmentOneOrigin(9.3749999999999956e-01, 1.9234735079187608e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(1.0409363447223732e+00, 1.8120621438385331e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoEnd(1.0000000000000000e+00, 2.2204460492503131e-16, 0.0000000000000000e+00);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
        }

        // check parallel intersection
        {
            Eigen::Vector3d segmentOneOrigin(6.9388939039072284e-18, 9.9999999999999989e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(5.2041704279304221e-18, 7.5000000000000000e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoEnd(6.9388939039072284e-18, 9.9999999999999989e-01, 0.0000000000000000e+00);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionLinesType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarParallel);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
            ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
            ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.75) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        }

        // check parallel intersection
        {
            Eigen::Vector3d segmentOneOrigin(2.3409475899275122e+02, 2.1032887616417568e+02, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(2.5749609565741929e+02, 1.8464458777130039e+02, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(2.5749609556478464e+02, 1.8464458790048394e+02, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoEnd(1.9112000445373033e+02, 2.5749609588054147e+02, 0.0000000000000000e+00);

            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
            Gedim::GeometryUtilities geometryUtilitiesLocal(geometryUtilitiesConfig);
            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilitiesLocal.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionLinesType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarParallel);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
            ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
            ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilitiesConfig.Tolerance1D);
            ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 0.35255671401423067) <
                        geometryUtilitiesConfig.Tolerance1D);
            ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        }

        // check parallel intersection
        {
            Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
            geometryUtilitiesConfig.Tolerance1D = 1.0e-13;
            Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

            Eigen::Vector3d segmentOneOrigin(5.5647058823529416e-01, 4.4235294117647062e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentOneEnd(5.5714285714285716e-01, 4.4285714285714284e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoOrigin(1.0000000000000001e-01, 1.0000000000000001e-01, 0.0000000000000000e+00);
            Eigen::Vector3d segmentTwoEnd(9.0000000000000002e-01, 6.9999999999999996e-01, 0.0000000000000000e+00);

            Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result =
                geometryUtilities.IntersectionSegmentSegment(segmentOneOrigin, segmentOneEnd, segmentTwoOrigin, segmentTwoEnd);
            ASSERT_EQ(result.IntersectionLinesType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarParallel);
            ASSERT_EQ(result.IntersectionSegmentsType,
                      Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionSegmentPlane)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check no intersection
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(1.0, 0.0, 3.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        }

        // check single intersection before origin
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 2.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
            ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, -1.1000000000000000e+01);
        }

        // check single intersection in origin
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 2.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(1.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 0.0000000000000000e+00);
        }

        // check single intersection inside segment
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 2.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(1.5, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 5.0000000000000000e-01);
        }

        // check single intersection in end
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 2.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(2.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 1.0000000000000000e+00);
        }

        // check single intersection after end
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 2.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(10.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
            ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 9.0000000000000000e+00);
        }

        // check multiple intersection
        {
            Eigen::Vector3d segmentOrigin(-10.0, 1.0, 13.0);
            Eigen::Vector3d segmentEnd(-10.0, 2.0, 13.0);

            Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
            Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
        }

        // check no intersection parallel case
        {
            Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
            Eigen::Vector3d segmentEnd(0.0, 1.0, 0.0);

            Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
            Eigen::Vector3d planeOrigin(0.0, 0.0, 2.0);

            Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result =
                geometryUtilities.IntersectionSegmentPlane(segmentOrigin, segmentEnd, planeNormal, planeOrigin);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronPlane)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-15;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // test no intersection with reference tetrahedron
        {
            Eigen::Vector3d v1(+0.0, +0.0, +0.0);
            Eigen::Vector3d v2(+1.0, +0.0, +0.0);
            Eigen::Vector3d v3(+0.0, +1.0, +0.0);
            Eigen::Vector3d v4(+0.0, +0.0, +1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
            Eigen::Vector3d planeOrigin(0.0, 0.0, 2.0);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::None);
        }

        // test intersection with vertex of reference tetrahedron
        {
            Eigen::Vector3d v1(+0.0, +0.0, +0.0);
            Eigen::Vector3d v2(+1.0, +0.0, +0.0);
            Eigen::Vector3d v3(+0.0, +1.0, +0.0);
            Eigen::Vector3d v4(+0.0, +0.0, +1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
            Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnVertex);
            ASSERT_EQ(result.IntersectionId, 3);
            ASSERT_EQ(result.VertexIntersections.size(), 4);
            ASSERT_EQ(result.VertexIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.EdgeIntersections.size(), 6);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.FaceIntersections.size(), 4);
            ASSERT_EQ(result.FaceIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        }

        // test intersection with edge of reference tetrahedron
        {
            Eigen::Vector3d v1(+0.0, +0.0, +0.0);
            Eigen::Vector3d v2(+1.0, +0.0, +0.0);
            Eigen::Vector3d v3(+0.0, +1.0, +0.0);
            Eigen::Vector3d v4(+0.0, +0.0, +1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            Eigen::Vector3d planeNormal(-1.0 / sqrt(2), -1.0 / sqrt(2), 0.0);
            Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnEdge);
            ASSERT_EQ(result.IntersectionId, 3);
            ASSERT_EQ(result.VertexIntersections.size(), 4);
            ASSERT_EQ(result.VertexIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.VertexIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.EdgeIntersections.size(), 6);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
            ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.FaceIntersections.size(), 4);
            ASSERT_EQ(result.FaceIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        }

        // test intersection with face of reference tetrahedron
        {
            Eigen::Vector3d v1(+0.0, +0.0, +0.0);
            Eigen::Vector3d v2(+1.0, +0.0, +0.0);
            Eigen::Vector3d v3(+0.0, +1.0, +0.0);
            Eigen::Vector3d v4(+0.0, +0.0, +1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            Eigen::Vector3d planeNormal(-1.0 / sqrt(3), -1.0 / sqrt(3), -1.0 / sqrt(3));
            Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnFace);
            ASSERT_EQ(result.IntersectionId, 3);
            ASSERT_EQ(result.VertexIntersections.size(), 4);
            ASSERT_EQ(result.VertexIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.VertexIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.VertexIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.EdgeIntersections.size(), 6);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
            ASSERT_EQ(result.FaceIntersections.size(), 4);
            ASSERT_EQ(result.FaceIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::Intersection);
        }

        // test triangular intersection with reference tetrahedron
        {
            Eigen::Vector3d v1(+0.0, +0.0, +0.0);
            Eigen::Vector3d v2(+1.0, +0.0, +0.0);
            Eigen::Vector3d v3(+0.0, +1.0, +0.0);
            Eigen::Vector3d v4(+0.0, +0.0, +1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            const string exportFolder = "./Export/TestIntersectionPolyhedronPlane/"
                                        "TriangularIntersectionWithReferenceTetra";
            Gedim::Output::CreateFolder(exportFolder);

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPoints(polyhedron.Vertices);
                exporter.Export(exportFolder + "/Vertices.vtu");
            }

            {
                Gedim::VTKUtilities exporter;
                exporter.AddSegments(polyhedron.Vertices, polyhedron.Edges);
                exporter.Export(exportFolder + "/Edges.vtu");
            }

            {
                vector<vector<unsigned int>> faces(polyhedron.Faces.size());
                for (unsigned int f = 0; f < faces.size(); f++)
                {
                    faces[f].resize(polyhedron.Faces[f].cols());
                    for (unsigned int v = 0; v < polyhedron.Faces[f].cols(); v++)
                        faces[f][v] = polyhedron.Faces[f](0, v);
                }
                Gedim::VTKUtilities exporter;
                exporter.AddPolygons(polyhedron.Vertices, faces);
                exporter.Export(exportFolder + "/Faces.vtu");
            }

            {
                Gedim::VTKUtilities exporter;
                exporter.AddPolyhedron(polyhedron.Vertices, polyhedron.Edges, polyhedron.Faces);
                exporter.Export(exportFolder + "/Polyhedron.vtu");
            }

            Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
            Eigen::Vector3d planeOrigin(0.0, 0.0, 0.25);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::NewPolygon);
            ASSERT_EQ(result.Intersections.size(), 3);
            ASSERT_EQ(result.IntersectionCoordinates.rows(), 3);
            ASSERT_EQ(result.IntersectionCoordinates.cols(), 3);
            ASSERT_EQ(result.Intersections.at(0).EdgeId, 4);
            ASSERT_EQ(result.Intersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
            ASSERT_EQ(result.Intersections.at(1).EdgeId, 5);
            ASSERT_EQ(result.Intersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
            ASSERT_EQ(result.Intersections.at(2).EdgeId, 3);
            ASSERT_EQ(result.Intersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
            ASSERT_EQ(result.IntersectionCoordinates.col(0), Eigen::Vector3d(0.75, 0.0, 0.25));
            ASSERT_EQ(result.IntersectionCoordinates.col(1), Eigen::Vector3d(0.0, 0.75, 0.25));
            ASSERT_EQ(result.IntersectionCoordinates.col(2), Eigen::Vector3d(0.0, 0.0, 0.25));
        }

        // test intersection with generic tetrahedron and generic plane
        {
            Eigen::Vector3d v1(+1.0, +1.0, +1.0);
            Eigen::Vector3d v2(-1.0, -1.0, +1.0);
            Eigen::Vector3d v3(-1.0, +1.0, -1.0);
            Eigen::Vector3d v4(+1.0, -1.0, -1.0);

            Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithVertices(v1, v2, v3, v4);

            Eigen::Vector3d planeNormal(-8.9087080637474766e-02, -4.4543540318737396e-01, 8.9087080637474791e-01);
            Eigen::Vector3d planeOrigin(1.4761904761904763e+00, 1.3809523809523809e+00, 1.2380952380952377e+00);
            Eigen::Matrix3d planeRotationMatrix = geometryUtilities.PlaneRotationMatrix(planeNormal).transpose();
            Eigen::Vector3d planeTranslation = geometryUtilities.PlaneTranslation(planeOrigin);

            Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result =
                geometryUtilities.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                              polyhedron.Edges,
                                                              polyhedron.Faces,
                                                              planeNormal,
                                                              planeOrigin,
                                                              planeRotationMatrix,
                                                              planeTranslation);

            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::NewPolygon);
            ASSERT_EQ(result.VertexIntersections.size(), 4);
            ASSERT_EQ(result.VertexIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
            ASSERT_EQ(result.VertexIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.VertexIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
            ASSERT_EQ(result.EdgeIntersections.size(), 6);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
            ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(2).Intersection.SingleIntersection.CurvilinearCoordinate, 0.40000000000000013);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.CurvilinearCoordinate, 0.54545454545454564);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type,
                      Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
            ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.Type,
                      Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
            ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.CurvilinearCoordinate, 2.2499999999999996);
            ASSERT_EQ(result.FaceIntersections.size(), 4);
            ASSERT_EQ(result.FaceIntersections.at(0).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(1).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(2).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
            ASSERT_EQ(result.FaceIntersections.at(3).Type,
                      Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionSegmentCircle)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check simply no intersection
        {
            Eigen::Vector3d segmentOrigin(1.0, 0.0, 0.0);
            Eigen::Vector3d segmentEnd(2.0, 0.0, 0.0);
            Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
            double circleRadius = 2.0;

            Gedim::GeometryUtilities::IntersectionSegmentCircleResult result =
                geometryUtilities.IntersectionSegmentCircle(segmentOrigin, segmentEnd, circleCenter, circleRadius);
            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentCircleResult::Types::NoIntersection);
        }

        // check simply one intersection
        {
            Eigen::Vector3d segmentOrigin(1.0, 5.0, 0.0);
            Eigen::Vector3d segmentEnd(2.0, 5.0, 0.0);
            Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
            double circleRadius = 2.0;

            Gedim::GeometryUtilities::IntersectionSegmentCircleResult result =
                geometryUtilities.IntersectionSegmentCircle(segmentOrigin, segmentEnd, circleCenter, circleRadius);
            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentCircleResult::Types::TangentIntersection);
            ASSERT_EQ(result.SegmentIntersections.size(), 1);
            ASSERT_DOUBLE_EQ(result.SegmentIntersections[0].CurvilinearCoordinate, -1.0);
            ASSERT_EQ(result.SegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
        }

        // check simply two intersections
        {
            Eigen::Vector3d segmentOrigin(1.0, 3.0, 0.0);
            Eigen::Vector3d segmentEnd(2.0, 3.0, 0.0);
            Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
            double circleRadius = 2.0;

            Gedim::GeometryUtilities::IntersectionSegmentCircleResult result =
                geometryUtilities.IntersectionSegmentCircle(segmentOrigin, segmentEnd, circleCenter, circleRadius);
            ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentCircleResult::Types::TwoIntersections);
            ASSERT_EQ(result.SegmentIntersections.size(), 2);
            ASSERT_DOUBLE_EQ(result.SegmentIntersections[0].CurvilinearCoordinate, -3.0);
            ASSERT_EQ(result.SegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
            ASSERT_DOUBLE_EQ(result.SegmentIntersections[1].CurvilinearCoordinate, 1.0);
            ASSERT_EQ(result.SegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolygonCircle)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // check simply no intersection
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
            double circleRadius = 1.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult result =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(result.Intersections.size(), 0);
        }

        // check Polygon Inside Circle with center outside
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
            double circleRadius = 10.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 0);
        }

        // Polygon Inside Circle with center outside one intersection
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(sqrt(2.0), sqrt(2.0), 0.0);
            double circleRadius = 2.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 1);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 0);
        }

        // check circle inside polygon
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
            double circleRadius = 0.125;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 0);
        }

        // check Polygon Inside Circle with center inside
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
            double circleRadius = 10.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 0);
        }

        // check Polygon Intersects Circle on border
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
            double circleRadius = sqrt(2) / 2;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 3);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 0);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Index, 1);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Index, 2);
        }

        // check Circle Intersect Polygon in one vertex
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.0, 2.0, 0.0);
            double circleRadius = 1.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 1);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 2);
        }

        // check Circle Intersect Polygon tangent to edge
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.5, -1.0, 0.0);
            double circleRadius = 1.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 1);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Tangent);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 0);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[0].CurvilinearCoordinate, 0.5);
        }

        // check generic intersections with triangle
        {
            Eigen::Matrix3d polygonVertices;
            polygonVertices.col(0) << 0.0, 0.0, 0.0;
            polygonVertices.col(1) << 1.0, 0.0, 0.0;
            polygonVertices.col(2) << 0.0, 1.0, 0.0;
            Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
            double circleRadius = 0.5;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 4);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Tangent);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 0);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[0].CurvilinearCoordinate, 0.5);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Index, 1);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[1].CurvilinearCoordinate, 0.5 * (1.0 - 1.0 / sqrt(2)));
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Index, 1);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[2].CurvilinearCoordinate, 1.0 - 0.5 * (1.0 - 1.0 / sqrt(2)));
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Tangent);
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].Index, 2);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[3].CurvilinearCoordinate, 0.5);
        }

        // check generic intersections with quadrilateral
        {
            Eigen::MatrixXd polygonVertices(3, 4);
            polygonVertices.col(0) << 1.0, 3.0, 0.0;
            polygonVertices.col(1) << 3.0, 3.0 - 2.0 / sqrt(3.0), 0.0;
            polygonVertices.col(2) << 4.0 + 1.0 / 10.0, 3.0, 0.0;
            polygonVertices.col(3) << 3.0, 4.0, 0.0;
            Eigen::Vector3d circleCenter(3.0, 3.0, 0.0);
            double circleRadius = 1.0;

            Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections =
                geometryUtilities.IntersectionPolygonCircle(polygonVertices, circleCenter, circleRadius);
            ASSERT_EQ(polygonCircleIntersections.Intersections.size(), 6);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Tangent);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[0].Index, 0);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[0].CurvilinearCoordinate, 0.74999999999999989);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[1].Index, 1);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[1].CurvilinearCoordinate, 0.14507270926633217);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[2].Index, 1);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[2].CurvilinearCoordinate, 0.90342008234572591);
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[3].Index, 2);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[3].CurvilinearCoordinate, 0.095022624434388858);
            ASSERT_EQ(polygonCircleIntersections.Intersections[4].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[4].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Vertex);
            ASSERT_EQ(polygonCircleIntersections.Intersections[4].Index, 3);
            ASSERT_EQ(polygonCircleIntersections.Intersections[5].Type,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::Types::Secant);
            ASSERT_EQ(polygonCircleIntersections.Intersections[5].IndexType,
                      Gedim::GeometryUtilities::IntersectionPolygonCircleResult::Intersection::IndexTypes::Edge);
            ASSERT_EQ(polygonCircleIntersections.Intersections[5].Index, 3);
            ASSERT_DOUBLE_EQ(polygonCircleIntersections.Intersections[5].CurvilinearCoordinate, 0.40000000000000002);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronLine_Cube_NoIntersection)
{
    // test no intersection with reference cube
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);

        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineTangent(1.0, 0.0, 0.0);
        Eigen::Vector3d lineOrigin(0.0, 0.0, 2.0);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::None);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronLine_Cube_OneIntersection_Vertex)
{
    // test single intersection with reference cube
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineTangent(1.0, 1.0, 1.0);
        Eigen::Vector3d lineOrigin(0.0, 0.0, 1.0);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::OneIntersection);
        ASSERT_EQ(result.LineIntersections.size(), 1);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[0].CurvilinearCoordinate,
                                                     0.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[0].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnVertex);
        ASSERT_EQ(result.LineIntersections[0].PolyhedronIndex, 4);
        ASSERT_EQ(result.PolyhedronVertexIntersections[4].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronVertexIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronVertexIntersections[4].LineIntersectionIndex, 0);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronLine_Cube_TwoIntersections_Faces)
{
    // test two intersections with reference cube
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineTangent(0.0, 0.0, 1.0);
        Eigen::Vector3d lineOrigin(0.5, 0.5, 0.0);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::TwoIntersections);
        ASSERT_EQ(result.LineIntersections.size(), 2);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[0].CurvilinearCoordinate,
                                                     0.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[0].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnFace);
        ASSERT_EQ(result.LineIntersections[0].PolyhedronIndex, 0);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[1].CurvilinearCoordinate,
                                                     1.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[1].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnFace);
        ASSERT_EQ(result.LineIntersections[1].PolyhedronIndex, 1);
        ASSERT_EQ(result.PolyhedronFaceIntersections[0].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronFaceIntersections[0].LineIntersectionIndex, 0);
        ASSERT_EQ(result.PolyhedronFaceIntersections[1].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronFaceIntersections[1].LineIntersectionIndex, 1);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}
TEST(TestGeometryUtilities, TestIntersectionPolyhedronLine_Cube_TwoIntersections_EdgeFace)
{
    // test two intersections with reference cube
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineTangent(1.0, 0.5, 0.0);
        Eigen::Vector3d lineOrigin(-1.0, -0.5, 0.5);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::TwoIntersections);
        ASSERT_EQ(result.LineIntersections.size(), 2);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[0].CurvilinearCoordinate,
                                                     1.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[0].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnEdge);
        ASSERT_EQ(result.LineIntersections[0].PolyhedronIndex, 8);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[1].CurvilinearCoordinate,
                                                     2.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[1].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnFace);
        ASSERT_EQ(result.LineIntersections[1].PolyhedronIndex, 3);
        ASSERT_EQ(result.PolyhedronEdgeIntersections[8].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronEdgeIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronEdgeIntersections[8].LineIntersectionIndex, 0);
        ASSERT_EQ(result.PolyhedronFaceIntersections[3].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronFaceIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronFaceIntersections[3].LineIntersectionIndex, 1);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronSegment_Cube_NoIntersection)
{
    // test no intersection with reference cube
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineOrigin(0.0, 0.0, 2.0);
        Eigen::Vector3d lineTangent(1.0, 0.0, 0.0);

        Eigen::Vector3d segmentOrigin = lineOrigin;
        Eigen::Vector3d segmentEnd = lineOrigin + lineTangent;
        Eigen::Vector3d segmentTangent = lineTangent;

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult polyhedronLineIntersections =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronSegment(polyhedron.Vertices,
                                                            polyhedron.Edges,
                                                            polyhedron.Faces,
                                                            segmentOrigin,
                                                            segmentEnd,
                                                            segmentTangent,
                                                            polyhedronLineIntersections);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::None);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronSegment_Cube_OneIntersection_Vertex)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineOrigin(0.0, 0.0, 1.0);
        Eigen::Vector3d lineTangent(1.0, 1.0, 1.0);

        Eigen::Vector3d segmentOrigin = lineOrigin;
        Eigen::Vector3d segmentEnd = lineOrigin + lineTangent;
        Eigen::Vector3d segmentTangent = lineTangent;

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult polyhedronLineIntersections =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronSegment(polyhedron.Vertices,
                                                            polyhedron.Edges,
                                                            polyhedron.Faces,
                                                            segmentOrigin,
                                                            segmentEnd,
                                                            segmentTangent,
                                                            polyhedronLineIntersections);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::OneIntersection);
        ASSERT_EQ(result.LineIntersections.size(), 1);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[0].CurvilinearCoordinate,
                                                     0.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[0].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnVertex);
        ASSERT_EQ(result.LineIntersections[0].PolyhedronIndex, 4);
        ASSERT_EQ(result.PolyhedronVertexIntersections[4].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronVertexIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronVertexIntersections[4].LineIntersectionIndex, 0);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronSegment_Cube_OneIntersection_Edge)
{
    // corrisponde al test "twointersections_edgeFace" ma solo un'intersezione è dentro il segmento considerato

    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        Eigen::Vector3d origin(+0.0, +0.0, +0.0);
        Eigen::Vector3d length(+1.0, +0.0, +0.0);
        Eigen::Vector3d width(+0.0, +1.0, +0.0);
        Eigen::Vector3d height(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron =
            geometryUtilities.CreateParallelepipedWithOrigin(origin, length, height, width);
        const vector<Eigen::MatrixXd> ployhedronFaceVertices =
            geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices, polyhedron.Faces);
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
        const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(ployhedronFaceVertices);
        const vector<bool> polyhedronFaceNormalDirections =
            geometryUtilities.PolyhedronFaceNormalDirections(ployhedronFaceVertices, polyhedronBarycenter, polyhedronFaceNormals);

        Eigen::Vector3d lineOrigin(-1.0, -0.5, 0.5);
        Eigen::Vector3d lineTangent(1.0, 0.5, 0.0);

        Eigen::Vector3d segmentOrigin = lineOrigin;
        Eigen::Vector3d segmentEnd = lineOrigin + lineTangent;
        Eigen::Vector3d segmentTangent = lineTangent;

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult polyhedronLineIntersections =
            geometryUtilities.IntersectionPolyhedronLine(polyhedron.Vertices,
                                                         polyhedron.Edges,
                                                         polyhedron.Faces,
                                                         polyhedronFaceNormals,
                                                         polyhedronFaceNormalDirections,
                                                         lineTangent,
                                                         lineOrigin);

        Gedim::GeometryUtilities::IntersectionPolyhedronLineResult result =
            geometryUtilities.IntersectionPolyhedronSegment(polyhedron.Vertices,
                                                            polyhedron.Edges,
                                                            polyhedron.Faces,
                                                            segmentOrigin,
                                                            segmentEnd,
                                                            segmentTangent,
                                                            polyhedronLineIntersections);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::Types::OneIntersection);
        ASSERT_EQ(result.LineIntersections.size(), 1);
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(result.LineIntersections[0].CurvilinearCoordinate,
                                                     1.0,
                                                     geometryUtilities.Tolerance1D()));
        ASSERT_EQ(result.LineIntersections[0].PolyhedronType,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::LineIntersection::Types::OnEdge);
        ASSERT_EQ(result.LineIntersections[0].PolyhedronIndex, 8);
        ASSERT_EQ(result.PolyhedronEdgeIntersections[8].Type,
                  Gedim::GeometryUtilities::IntersectionPolyhedronLineResult::PolyhedronEdgeIntersection::Types::Intersection);
        ASSERT_EQ(result.PolyhedronEdgeIntersections[8].LineIntersectionIndex, 0);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestIntersectionPolyhedronSegment_SimpleHexahedronMesh)
{
    // test no intersection with simple hexahedron mesh
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // create cell3Ds
        const unsigned int numCell3Ds = 4;
        vector<Gedim::GeometryUtilities::Polyhedron> cell3Ds(numCell3Ds);

        cell3Ds[0] = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                      Eigen::Vector3d(0.0, 0.5, 0.0));
        cell3Ds[1] = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                      Eigen::Vector3d(0.0, 0.5, 0.0));
        cell3Ds[2] = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0, 0.5, 0.0),
                                                                      Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                      Eigen::Vector3d(0.0, 0.5, 0.0));
        cell3Ds[3] = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                                      Eigen::Vector3d(0.5, 0.0, 0.0),
                                                                      Eigen::Vector3d(0.0, 0.0, 1.0),
                                                                      Eigen::Vector3d(0.0, 0.5, 0.0));

        vector<vector<Eigen::MatrixXd>> cell3DsFaceVertices(numCell3Ds);
        vector<Eigen::Vector3d> cell3DsBarycenters(numCell3Ds);
        vector<vector<Eigen::Vector3d>> cell3DsFaceNormals(numCell3Ds);
        vector<vector<bool>> cell3DsFaceNormalDirections(numCell3Ds);

        for (unsigned int c = 0; c < numCell3Ds; c++)
        {
            cell3DsFaceVertices[c] = geometryUtilities.PolyhedronFaceVertices(cell3Ds[c].Vertices, cell3Ds[c].Faces);
            cell3DsBarycenters[c] = geometryUtilities.PolyhedronBarycenter(cell3Ds[c].Vertices);
            cell3DsFaceNormals[c] = geometryUtilities.PolyhedronFaceNormals(cell3DsFaceVertices[c]);
            cell3DsFaceNormalDirections[c] = geometryUtilities.PolyhedronFaceNormalDirections(cell3DsFaceVertices[c],
                                                                                              cell3DsBarycenters[c],
                                                                                              cell3DsFaceNormals[c]);
        }

        // create segments
        const unsigned int numSegments = 4;
        Eigen::MatrixXd segmentOrigins(3, numSegments);
        Eigen::MatrixXd segmentEnds(3, numSegments);
        Eigen::MatrixXd segmentTagents(3, numSegments);

        segmentOrigins.col(0) << 0.25, 0.25, 1.0;
        segmentEnds.col(0) << 0.75, 0.75, 0.0;

        segmentOrigins.col(1) << 0.1, 0.1, 0.5;
        segmentEnds.col(1) << 0.0, 0.5, 0.0;

        segmentOrigins.col(2) << 0.9, 0.9, 0.95;
        segmentEnds.col(2) << 0.95, 0.95, 0.65;

        segmentOrigins.col(3) << 0.8, 0.0, 0.6;
        segmentEnds.col(3) << 0.8, 1.0, 0.6;

        for (unsigned int s = 0; s < numSegments; s++)
        {
            segmentTagents.col(s) = geometryUtilities.SegmentTangent(segmentOrigins.col(s), segmentEnds.col(s));
        }

        // intersects
        vector<Gedim::GeometryUtilities::IntersectionPolyhedronsSegmentResult> result(numSegments);
        for (unsigned int s = 0; s < numSegments; s++)
        {
            result[s] = geometryUtilities.IntersectionPolyhedronsSegment(cell3Ds,
                                                                         cell3DsFaceNormals,
                                                                         cell3DsFaceNormalDirections,
                                                                         segmentOrigins.col(s),
                                                                         segmentEnds.col(s),
                                                                         segmentTagents.col(s));
        }

        // check result
        ASSERT_EQ(result.size(), numSegments);
        ASSERT_EQ(result[0].Points.size(), 3);
        ASSERT_EQ(result[0].Segments.size(), 2);
        ASSERT_NE(result[0].Points.find(0.0), result[0].Points.end());
        ASSERT_NE(result[0].Points.find(0.5), result[0].Points.end());
        ASSERT_NE(result[0].Points.find(1.0), result[0].Points.end());
        ASSERT_EQ(result[0].Points.at(0.0).Cell3DIndices, vector<unsigned int>({0}));
        ASSERT_EQ(result[0].Points.at(0.5).Cell3DIndices, vector<unsigned int>({0, 1, 2, 3}));
        ASSERT_EQ(result[0].Points.at(1.0).Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[0].Segments[0].Cell3DIndices, vector<unsigned int>({0}));
        ASSERT_EQ(result[0].Segments[0].Points, vector<double>({0.0, 0.5}));
        ASSERT_EQ(result[0].Segments[1].Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[0].Segments[1].Points, vector<double>({0.5, 1.0}));

        ASSERT_EQ(result[1].Points.size(), 2);
        ASSERT_EQ(result[1].Segments.size(), 1);
        ASSERT_NE(result[1].Points.find(0.0), result[0].Points.end());
        ASSERT_NE(result[1].Points.find(1.0), result[0].Points.end());
        ASSERT_EQ(result[1].Points.at(0.0).Cell3DIndices, vector<unsigned int>({0}));
        ASSERT_EQ(result[1].Points.at(1.0).Cell3DIndices, vector<unsigned int>({0, 2}));
        ASSERT_EQ(result[1].Segments[0].Cell3DIndices, vector<unsigned int>({0}));
        ASSERT_EQ(result[1].Segments[0].Points, vector<double>({0.0, 1.0}));

        ASSERT_EQ(result[2].Points.size(), 2);
        ASSERT_EQ(result[2].Segments.size(), 1);
        ASSERT_NE(result[2].Points.find(0.0), result[0].Points.end());
        ASSERT_NE(result[2].Points.find(1.0), result[0].Points.end());
        ASSERT_EQ(result[2].Points.at(0.0).Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[2].Points.at(1.0).Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[2].Segments[0].Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[2].Segments[0].Points, vector<double>({0.0, 1.0}));

        ASSERT_EQ(result[3].Points.size(), 3);
        ASSERT_EQ(result[3].Segments.size(), 2);
        ASSERT_NE(result[3].Points.find(0.0), result[0].Points.end());
        ASSERT_NE(result[3].Points.find(0.5), result[0].Points.end());
        ASSERT_NE(result[3].Points.find(1.0), result[0].Points.end());
        ASSERT_EQ(result[3].Points.at(0.0).Cell3DIndices, vector<unsigned int>({1}));
        ASSERT_EQ(result[3].Points.at(0.5).Cell3DIndices, vector<unsigned int>({1, 3}));
        ASSERT_EQ(result[3].Points.at(1.0).Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[3].Segments[0].Cell3DIndices, vector<unsigned int>({1}));
        ASSERT_EQ(result[3].Segments[0].Points, vector<double>({0.0, 0.5}));
        ASSERT_EQ(result[3].Segments[1].Cell3DIndices, vector<unsigned int>({3}));
        ASSERT_EQ(result[3].Segments[1].Points, vector<double>({0.5, 1.0}));
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestGeometryUtilities, TestCheckTrianglesIntersection)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 5.0, 0.0, 0.0;
        triangle_one.col(2) << 0.0, 5.0, 0.0;

        triangle_two.col(0) << 0.0, 0.0, 0.0;
        triangle_two.col(1) << 5.0, 0.0, 0.0;
        triangle_two.col(2) << 0.0, 6.0, 0.0;

        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 5.0, 0.0, 0.0;
        triangle_one.col(2) << 0.0, 5.0, 0.0;

        triangle_two.col(0) << 0.0, 0.0, 0.0;
        triangle_two.col(1) << 5.0, 0.0, 0.0;
        triangle_two.col(2) << 0.0, 5.0, 0.0;

        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 5.0, 0.0, 0.0;
        triangle_one.col(2) << 0.0, 5.0, 0.0;

        triangle_two.col(0) << -10.0, 0.0, 0.0;
        triangle_two.col(1) << -5.0, 0.0, 0.0;
        triangle_two.col(2) << -1.0, 6.0, 0.0;

        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 5.0, 0.0, 0.0;
        triangle_one.col(2) << 2.5, 5.0, 0.0;

        triangle_two.col(0) << 0.0, 4.0, 0.0;
        triangle_two.col(1) << 2.5, -1.0, 0.0;
        triangle_two.col(2) << 5.0, 4.0, 0.0;

        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 1.0, 1.0, 0.0;
        triangle_one.col(2) << 0.0, 2.0, 0.0;

        triangle_two.col(0) << 2.0, 1.0, 0.0;
        triangle_two.col(1) << 3.0, 0.0, 0.0;
        triangle_two.col(2) << 3.0, 2.0, 0.0;

        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 1.0, 1.0, 0.0;
        triangle_one.col(2) << 0.0, 2.0, 0.0;

        triangle_two.col(0) << 2.0, 1.0, 0.0;
        triangle_two.col(1) << 3.0, -2.0, 0.0;
        triangle_two.col(2) << 3.0, 4.0, 0.0;

        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }

    {
        Eigen::Matrix3d triangle_one, triangle_two;

        triangle_one.col(0) << 0.0, 0.0, 0.0;
        triangle_one.col(1) << 1.0, 0.0, 0.0;
        triangle_one.col(2) << 0.0, 1.0, 0.0;

        triangle_two.col(0) << 1.0, 0.0, 0.0;
        triangle_two.col(1) << 2.0, 0.0, 0.0;
        triangle_two.col(2) << 1.0, 1.0, 0.0;

        ASSERT_TRUE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, true));
        ASSERT_FALSE(geometryUtilities.CheckTrianglesIntersection(triangle_one, triangle_two, false));
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_GEOMETRY_INTERSECTION_H
