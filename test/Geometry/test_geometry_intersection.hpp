#ifndef __TEST_GEOMETRY_INTERSECTION_H
#define __TEST_GEOMETRY_INTERSECTION_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestIntersectionSegmentSegment)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check simply no intersection
      {
        Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentOneEnd(   1.0, 0.0, 0.0);
        Eigen::Vector3d segmentTwoOrigin(0.0, 0.0, 2.0);
        Eigen::Vector3d segmentTwoEnd(   1.0, 0.0, 2.0);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
      }

      // check no coplanarity
      {
        Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.0 );
        Eigen::Vector3d segmentOneEnd(   1.0, 0.0, 0.0 );
        Eigen::Vector3d segmentTwoOrigin(0.0, 0.0, 0.25);
        Eigen::Vector3d segmentTwoEnd(   0.0, 1.0, 0.25);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
      }

      // check MultipleIntersections, no inclusion
      {
        Eigen::Vector3d segmentOneOrigin(0.25, 0.75, 0.0);
        Eigen::Vector3d segmentOneEnd(   0.75, 0.25, 0.0);
        Eigen::Vector3d segmentTwoOrigin(0.5, 0.5, 0.0  );
        Eigen::Vector3d segmentTwoEnd(   0.0, 1.0, 0.0  );

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);

        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
        ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
        ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 0.5) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 0.5) < geometryUtilityConfig.Tolerance);
        ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
      }

      // check MultipleIntersections, no inclusion
      {
        Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25 );
        Eigen::Vector3d segmentOneEnd(   0.75, 0.0, 0.25);
        Eigen::Vector3d segmentTwoOrigin(2.0, 0.0, 0.25 );
        Eigen::Vector3d segmentTwoEnd(   0.5, 0.0, 0.25 );

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);

        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
        ASSERT_EQ(result.SecondIntersectionRelation[0], 1);
        ASSERT_EQ(result.SecondIntersectionRelation[1], 0);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 2.0 / 3.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 5.0 / 6.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilityConfig.Tolerance);
        ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
      }

      // check MultipleIntersections, with inclusion
      {
        Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25 );
        Eigen::Vector3d segmentOneEnd(   1.0, 0.0, 0.25 );
        Eigen::Vector3d segmentTwoOrigin(0.5, 0.0, 0.25 );
        Eigen::Vector3d segmentTwoEnd(   0.75, 0.0, 0.25);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);

        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::MultipleIntersections);
        ASSERT_EQ(result.SecondIntersectionRelation[0], 0);
        ASSERT_EQ(result.SecondIntersectionRelation[1], 1);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[1].CurvilinearCoordinate - 0.75) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[1].CurvilinearCoordinate - 1.0) < geometryUtilityConfig.Tolerance);
        ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.FirstSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_EQ(result.SecondSegmentIntersections[1].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
      }


      // check OneIntersection
      {
        Eigen::Vector3d segmentOneOrigin(-1.0, 0.0, 0.25);
        Eigen::Vector3d segmentOneEnd(   1.0, 0.0, 0.25 );
        Eigen::Vector3d segmentTwoOrigin(0.0, -1.0, 0.25);
        Eigen::Vector3d segmentTwoEnd(   0.0, 1.0, 0.25 );

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 0.5) < geometryUtilityConfig.Tolerance);
        ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
      }

      // check OneIntersection origin
      {
        Eigen::Vector3d segmentOneOrigin(0.0, 0.0, 0.25 );
        Eigen::Vector3d segmentOneEnd(   1.0, 0.0, 0.25 );
        Eigen::Vector3d segmentTwoOrigin(0.0, -1.0, 0.25);
        Eigen::Vector3d segmentTwoEnd(   0.0, 0.0, 0.25 );

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
        ASSERT_TRUE(abs(result.FirstSegmentIntersections[0].CurvilinearCoordinate - 0.0) < geometryUtilityConfig.Tolerance);
        ASSERT_TRUE(abs(result.SecondSegmentIntersections[0].CurvilinearCoordinate - 1.0) < geometryUtilityConfig.Tolerance);
        ASSERT_EQ(result.FirstSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_EQ(result.SecondSegmentIntersections[0].Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
      }

      // check parallel segments
      {
        Eigen::Vector3d segmentOneOrigin(5.0000000000000056e-01, 1.1458333333333328e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentOneEnd(   4.9999999999999994e-01, 6.2500000000000022e-01, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoOrigin(9.9999999999999978e-01, 1.9999999999999996e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoEnd(   9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
      }

      // check parallel segments
      {
        Eigen::Vector3d segmentOneOrigin(1.1249999999999996e+00, 2.0000000000000000e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentOneEnd(   1.1249999999999996e+00, 1.8602389358491287e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
        Eigen::Vector3d segmentTwoEnd(   1.0000000000000000e+00, 2.2204460492503131e-16, 2.4825341532472739e-17);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
      }

      // check intersection in first segment end and second segment origin
      {
        Eigen::Vector3d segmentOneOrigin(1.2499999999999996e+00, 1.7105836549996734e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentOneEnd(   9.9999999999999956e-01, 2.0000000000000000e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
        Eigen::Vector3d segmentTwoEnd(   1.0000000000000000e+00, 2.2204460492503131e-16, 2.4825341532472739e-17);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
      }

      // check intersection inside segments
      {
        Eigen::Vector3d segmentOneOrigin(9.3749999999999956e-01, 1.9234735079187608e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentOneEnd(   1.0409363447223732e+00, 1.8120621438385331e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoOrigin(9.9999999999999956e-01, 2.0000000000000000e+00, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoEnd(   1.0000000000000000e+00, 2.2204460492503131e-16, 0.0000000000000000e+00);

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::SingleIntersection);
      }

      // check no intersection
      {
        Eigen::Vector3d segmentOneOrigin(3.6586807759999989e+00, 1.6586807760000004e+00, 0.0000000000000000e+00 );
        Eigen::Vector3d segmentOneEnd(   3.6586807759999997e+00, 1.6586807760000004e+00, 0.0000000000000000e+00 );
        Eigen::Vector3d segmentTwoOrigin(1.1586807759999997e+00, -2.3999999858728053e-08, 0.0000000000000000e+00);
        Eigen::Vector3d segmentTwoEnd(   1.1586807759999993e+00, 3.3173615759999997e+00, 0.0000000000000000e+00 );

        Gedim::GeometryUtilities::IntersectionSegmentSegmentResult result = geometryUtility.IntersectionSegmentSegment(segmentOneOrigin,
                                                                                                                       segmentOneEnd,
                                                                                                                       segmentTwoOrigin,
                                                                                                                       segmentTwoEnd);
        ASSERT_EQ(result.IntersectionLinesType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionLineTypes::CoPlanarIntersecting);
        ASSERT_EQ(result.IntersectionSegmentsType, Gedim::GeometryUtilities::IntersectionSegmentSegmentResult::IntersectionSegmentTypes::NoIntersection);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestIntersectionSegmentPlane)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check no intersection
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (1.0, 0.0, 3.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
      }

      // check single intersection before origin
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (2.0, 0.0, 2.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
        ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, -1.1000000000000000e+01);
      }

      // check single intersection in origin
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (2.0, 0.0, 2.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(1.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 0.0000000000000000e+00);
      }

      // check single intersection inside segment
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (2.0, 0.0, 2.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(1.5, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 5.0000000000000000e-01);
      }

      // check single intersection in end
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (2.0, 0.0, 2.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(2.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 1.0000000000000000e+00);
      }

      // check single intersection after end
      {
        Eigen::Vector3d segmentOrigin(1.0, 0.0, 2.0);
        Eigen::Vector3d segmentEnd   (2.0, 0.0, 2.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(10.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
        ASSERT_DOUBLE_EQ(result.SingleIntersection.CurvilinearCoordinate, 9.0000000000000000e+00);
      }

      // check multiple intersection
      {
        Eigen::Vector3d segmentOrigin(-10.0, 1.0, 13.0);
        Eigen::Vector3d segmentEnd   (-10.0, 2.0, 13.0);

        Eigen::Vector3d planeNormal(1.0, 0.0, 0.0);
        Eigen::Vector3d planeOrigin(-10.0, 11.0, 13.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
      }

      // check no intersection parallel case
      {
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd   (0.0, 1.0, 0.0);

        Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        Eigen::Vector3d planeOrigin(0.0, 0.0, 2.0);

        Gedim::GeometryUtilities::IntersectionSegmentPlaneResult result = geometryUtility.IntersectionSegmentPlane(segmentOrigin,
                                                                                                                   segmentEnd,
                                                                                                                   planeNormal,
                                                                                                                   planeOrigin);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }


  TEST(TestGeometryUtilities, TestIntersectionPolyhedronPlane)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // test no intersection with reference tetrahedron
      {
        Eigen::Vector3d v1(+0.0, +0.0, +0.0);
        Eigen::Vector3d v2(+1.0, +0.0, +0.0);
        Eigen::Vector3d v3(+0.0, +1.0, +0.0);
        Eigen::Vector3d v4(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtility.CreateTetrahedronWithVertices(v1,v2,v3,v4);

        Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        Eigen::Vector3d planeOrigin(0.0, 0.0, 2.0);
        Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal).transpose();

        Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result = geometryUtility.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::None);
      }

      // test intersection with vertex of reference tetrahedron
      {
        Eigen::Vector3d v1(+0.0, +0.0, +0.0);
        Eigen::Vector3d v2(+1.0, +0.0, +0.0);
        Eigen::Vector3d v3(+0.0, +1.0, +0.0);
        Eigen::Vector3d v4(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtility.CreateTetrahedronWithVertices(v1,v2,v3,v4);

        Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
        Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal).transpose();

        Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result = geometryUtility.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnVertex);
        ASSERT_EQ(result.IntersectionId, 3);
        ASSERT_EQ(result.VertexIntersections.size(), 4);
        ASSERT_EQ(result.VertexIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.EdgeIntersections.size(), 6);
        ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.FaceIntersections.size(), 4);
        ASSERT_EQ(result.FaceIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
      }

      // test intersection with edge of reference tetrahedron
      {
        Eigen::Vector3d v1(+0.0, +0.0, +0.0);
        Eigen::Vector3d v2(+1.0, +0.0, +0.0);
        Eigen::Vector3d v3(+0.0, +1.0, +0.0);
        Eigen::Vector3d v4(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtility.CreateTetrahedronWithVertices(v1,v2,v3,v4);

        Eigen::Vector3d planeNormal(-1.0 / sqrt(2), -1.0 / sqrt(2), 0.0);
        Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
        Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal).transpose();

        Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result = geometryUtility.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnEdge);
        ASSERT_EQ(result.IntersectionId, 3);
        ASSERT_EQ(result.VertexIntersections.size(), 4);
        ASSERT_EQ(result.VertexIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.VertexIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.EdgeIntersections.size(), 6);
        ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
        ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.CurvilinearCoordinate, 0.0);
        ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::NoIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
        ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(4).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(5).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.FaceIntersections.size(), 4);
        ASSERT_EQ(result.FaceIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
      }

      // test intersection with face of reference tetrahedron
      {
        Eigen::Vector3d v1(+0.0, +0.0, +0.0);
        Eigen::Vector3d v2(+1.0, +0.0, +0.0);
        Eigen::Vector3d v3(+0.0, +1.0, +0.0);
        Eigen::Vector3d v4(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtility.CreateTetrahedronWithVertices(v1,v2,v3,v4);

        Eigen::Vector3d planeNormal(-1.0 / sqrt(3), -1.0 / sqrt(3), -1.0 / sqrt(3));
        Eigen::Vector3d planeOrigin(0.0, 0.0, 1.0);
        Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal).transpose();

        Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result = geometryUtility.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::OnFace);
        ASSERT_EQ(result.IntersectionId, 3);
        ASSERT_EQ(result.VertexIntersections.size(), 4);
        ASSERT_EQ(result.VertexIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::NoIntersection);
        ASSERT_EQ(result.VertexIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.VertexIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.VertexIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::VertexIntersection::Types::Intersection);
        ASSERT_EQ(result.EdgeIntersections.size(), 6);
        ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(0).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(1).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(2).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
        ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::SingleIntersection);
        ASSERT_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.Type, Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(result.EdgeIntersections.at(3).Intersection.SingleIntersection.CurvilinearCoordinate, 1.0);
        ASSERT_EQ(result.EdgeIntersections.at(4).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
        ASSERT_EQ(result.EdgeIntersections.at(5).Intersection.Type, Gedim::GeometryUtilities::IntersectionSegmentPlaneResult::Types::MultipleIntersections);
        ASSERT_EQ(result.FaceIntersections.size(), 4);
        ASSERT_EQ(result.FaceIntersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::NoIntersection);
        ASSERT_EQ(result.FaceIntersections.at(3).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::FaceIntersection::Types::Intersection);
      }

      // test triangular intersection with reference tetrahedron
      {
        Eigen::Vector3d v1(+0.0, +0.0, +0.0);
        Eigen::Vector3d v2(+1.0, +0.0, +0.0);
        Eigen::Vector3d v3(+0.0, +1.0, +0.0);
        Eigen::Vector3d v4(+0.0, +0.0, +1.0);

        Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtility.CreateTetrahedronWithVertices(v1,v2,v3,v4);

        Eigen::Vector3d planeNormal(0.0, 0.0, 1.0);
        Eigen::Vector3d planeOrigin(0.0, 0.0, 0.25);
        Eigen::Matrix3d planeRotationMatrix = geometryUtility.PlaneRotationMatrix(planeNormal).transpose();

        Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult result = geometryUtility.IntersectionPolyhedronPlane(polyhedron.Vertices,
                                                                                                                         polyhedron.Edges,
                                                                                                                         polyhedron.Faces,
                                                                                                                         planeNormal,
                                                                                                                         planeOrigin,
                                                                                                                         planeRotationMatrix);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Types::NewPolygon);
        ASSERT_EQ(result.Intersections.size(), 3);
        ASSERT_EQ(result.IntersectionCoordinates.rows(), 3);
        ASSERT_EQ(result.IntersectionCoordinates.cols(), 3);
        ASSERT_EQ(result.Intersections.at(0).EdgeId, 3);
        ASSERT_EQ(result.Intersections.at(0).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
        ASSERT_EQ(result.Intersections.at(1).EdgeId, 4);
        ASSERT_EQ(result.Intersections.at(1).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
        ASSERT_EQ(result.Intersections.at(2).EdgeId, 5);
        ASSERT_EQ(result.Intersections.at(2).Type, Gedim::GeometryUtilities::IntersectionPolyhedronPlaneResult::Intersection::Types::Edge);
        ASSERT_EQ(result.IntersectionCoordinates.col(0), Eigen::Vector3d(0.0, 0.0, 0.25));
        ASSERT_EQ(result.IntersectionCoordinates.col(1), Eigen::Vector3d(0.75, 0.0, 0.25));
        ASSERT_EQ(result.IntersectionCoordinates.col(2), Eigen::Vector3d(0.0, 0.75, 0.25));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_INTERSECTION_H
