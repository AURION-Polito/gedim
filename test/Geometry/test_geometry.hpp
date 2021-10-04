#ifndef __TEST_GEOMETRY_H
#define __TEST_GEOMETRY_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

  TEST(TestGeometryUtilities, TestCompareValues)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check Compare1DValues
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare1DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare2DValues
      {
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare2DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare2DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check Compare3DValues
      {
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, 1.0), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare3DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
        ASSERT_EQ(geometryUtility.Compare1DValues(0.0, geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance), Gedim::GeometryUtilities::CompareTypes::Coincident);
        ASSERT_EQ(geometryUtility.Compare3DValues(1.0, -1.0), Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(0.0));
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue1DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(2 * geometryUtilityConfig.Tolerance));
        ASSERT_TRUE(geometryUtility.IsValue1DPositive(10.0));
      }

      // check IsAreaPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(0.0));
        ASSERT_TRUE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue2DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue2DPositive(10.0));
      }

      // check IsLenghtPositive
      {
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(0.0));
        ASSERT_TRUE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance * geometryUtilityConfig.Tolerance));
        ASSERT_FALSE(geometryUtility.IsValue3DPositive(-1.0));
        ASSERT_TRUE(geometryUtility.IsValue3DPositive(10.0));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }


  TEST(TestGeometryUtilities, TestPointPointLinePosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check curvilinear coordinate before
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(-0.5, 1.5),
                  Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
      }

      // check curvilinear coordinate coincident
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(0.5, 0.5 + geometryUtilityConfig.Tolerance / 2.0),
                  Gedim::GeometryUtilities::CompareTypes::Coincident);
      }

      // check curvilinear coordinate before
      {
        ASSERT_EQ(geometryUtility.Compare1DValues(0.5, -1.5),
                  Gedim::GeometryUtilities::CompareTypes::SecondBeforeFirst);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPointCurvilinearCoordinate)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check curvilinear coordinate origin
      {
        Eigen::Vector3d point(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);

        ASSERT_TRUE(abs(geometryUtility.PointCurvilinearCoordinate(point,
                                                                   segmentOrigin,
                                                                   segmentEnd) - 0.0) < geometryUtilityConfig.Tolerance);
      }

      // check curvilinear coordinate end
      {
        Eigen::Vector3d point(10.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);

        ASSERT_TRUE(abs(geometryUtility.PointCurvilinearCoordinate(point,
                                                                   segmentOrigin,
                                                                   segmentEnd) - 1.0) < geometryUtilityConfig.Tolerance);
      }

      // check curvilinear coordinate inside
      {
        Eigen::Vector3d point(5.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);

        ASSERT_TRUE(abs(geometryUtility.PointCurvilinearCoordinate(point,
                                                                   segmentOrigin,
                                                                   segmentEnd) - 0.5) < geometryUtilityConfig.Tolerance);
      }

      // check curvilinear coordinate before origin
      {
        Eigen::Vector3d point(-5.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);

        ASSERT_TRUE(abs(geometryUtility.PointCurvilinearCoordinate(point,
                                                                   segmentOrigin,
                                                                   segmentEnd) - -0.5) < geometryUtilityConfig.Tolerance);
      }

      // check curvilinear coordinate after end
      {
        Eigen::Vector3d point(15.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);

        ASSERT_TRUE(abs(geometryUtility.PointCurvilinearCoordinate(point,
                                                                   segmentOrigin,
                                                                   segmentEnd) - 1.5) < geometryUtilityConfig.Tolerance);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, IntersectionSegmentSegment)
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

  TEST(TestGeometryUtilities, PointSegmentPosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check right
      {
        Eigen::Vector3d point(0.5, -1.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::RightTheSegment);
      }

      // check left
      {
        Eigen::Vector3d point(0.5, 1.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment);
      }

      // check on segment line before origin
      {
        Eigen::Vector3d point(-10.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
      }

      // check on segment line after end
      {
        Eigen::Vector3d point(10.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
      }

      // check on segment origin
      {
        Eigen::Vector3d point(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
      }

      // check on segment end
      {
        Eigen::Vector3d point(1.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
      }

      // check inside segment
      {
        Eigen::Vector3d point(0.5, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, PointPolygonPosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check outside
      {
        Eigen::Vector3d point(0.5, -1.0, 0.0);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtility.PointPolygonPosition(point,
                                                                                                           polygonVertices);
        ASSERT_EQ(result.PositionType, Gedim::GeometryUtilities::PointPolygonPositionResult::PositionTypes::Outside);
      }

      // check border
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        // border edge
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderEdge= geometryUtility.PointPolygonPosition(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                                                                                    polygonVertices);
        ASSERT_EQ(resultBorderEdge.PositionType, Gedim::GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderEdge);
        ASSERT_EQ(resultBorderEdge.BorderIndex, 1);

        // border vertex
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderVertex = geometryUtility.PointPolygonPosition(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                                       polygonVertices);
        ASSERT_EQ(resultBorderVertex.PositionType, Gedim::GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderVertex);
        ASSERT_EQ(resultBorderVertex.BorderIndex, 1);
      }

      // check inside
      {
        Eigen::Vector3d point(0.50, 0.25, 0.0);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtility.PointPolygonPosition(point,
                                                                                                           polygonVertices);
        ASSERT_EQ(result.PositionType, Gedim::GeometryUtilities::PointPolygonPositionResult::PositionTypes::Inside);
      }

      // check on vertex
      {
        Eigen::Vector3d point(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.row(0)<< 9.9999999999999956e-01, 9.3749999999999956e-01, 1.0409363447223732e+00;
        polygonVertices.row(1)<< 2.0000000000000000e+00, 1.9234735079187608e+00, 1.8120621438385331e+00;
        polygonVertices.row(2)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtility.PointPolygonPosition(point,
                                                                                                           polygonVertices);
        ASSERT_EQ(result.PositionType, Gedim::GeometryUtilities::PointPolygonPositionResult::PositionTypes::BorderVertex);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, SplitPolygon)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // no action
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::NoAction);
        ASSERT_EQ(result.NewVertices.size(), 0);
        ASSERT_EQ(result.NewEdges.size(), 0);
        ASSERT_EQ(result.NewPolygons.size(), 0);
      }

      // update square with one point
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 2);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 5);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 5);
      }

      // update square with two points
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 6);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 6);
      }

      // update square with two points with aligned edges
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 5;
        input.AlignedEdges.resize(1);
        input.AlignedEdges[0].OriginVertexIndex = 1;
        input.AlignedEdges[0].EndVertexIndex = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 2;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonUpdate);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 4);
        ASSERT_EQ(result.NewPolygons.size(), 1);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 7);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 7);
      }

      // split square in two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 0);
        ASSERT_EQ(result.NewEdges.size(), 1);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
      }

      // split square in a triangle and a trapezioid
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 4);
      }

      // split square in a two trapezioids
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 2);
        ASSERT_EQ(result.NewEdges.size(), 5);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 4);
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 1;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 0;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 0);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 1);
        ASSERT_EQ(newEdgeIter->OriginId, 1);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 1);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 2);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 3, 5, 2 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 3, 0 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 4, 3, 0 }));
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 2);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 0);
        ASSERT_EQ(newEdgeIter->OriginId, 0);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 0);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 1);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 4, 3, 2 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 1, 3, 5 }));
      }

      // split triangle in a two triangles
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 3;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 2;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 1;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewVertices.front().Type, Gedim::GeometryUtilities::SplitPolygonResult::NewVertex::Types::SegmentOrigin);
        ASSERT_EQ(result.NewEdges.size(), 3);
        auto newEdgeIter = result.NewEdges.begin();
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeNew);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 1);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 2);
        ASSERT_EQ(newEdgeIter->OriginId, 2);
        ASSERT_EQ(newEdgeIter->EndId, 3);
        newEdgeIter++;
        ASSERT_EQ(newEdgeIter->Type, Gedim::GeometryUtilities::SplitPolygonResult::NewEdge::Types::EdgeUpdate);
        ASSERT_EQ(newEdgeIter->OldEdgeId, 2);
        ASSERT_EQ(newEdgeIter->OriginId, 3);
        ASSERT_EQ(newEdgeIter->EndId, 0);

        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 1, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 0, 3, 5 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 2, 3, 1 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 4, 3, 1 }));
      }

      // split square in a a triangle and a trapezioid
      {
        Gedim::GeometryUtilities::SplitPolygonInput input;
        input.NumberPolygonVertices = 4;
        input.Segment.Origin.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Edge;
        input.Segment.Origin.Index = 0;
        input.Segment.End.Type = Gedim::GeometryUtilities::SplitPolygonInput::SplitSegment::Vertex::Types::Vertex;
        input.Segment.End.Index = 2;

        Gedim::GeometryUtilities::SplitPolygonResult result = geometryUtility.SplitPolygon(input);

        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::SplitPolygonResult::Types::PolygonCreation);
        ASSERT_EQ(result.NewVertices.size(), 1);
        ASSERT_EQ(result.NewEdges.size(), 3);
        ASSERT_EQ(result.NewPolygons.size(), 2);
        ASSERT_EQ(result.NewPolygons[0].Vertices.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Vertices, list<unsigned int>({ 0, 4, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[0].Edges.size(), 4);
        ASSERT_EQ(result.NewPolygons[0].Edges, list<unsigned int>({ 5, 4, 2, 3 }));
        ASSERT_EQ(result.NewPolygons[1].Vertices.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Vertices, list<unsigned int>({ 1, 2, 4 }));
        ASSERT_EQ(result.NewPolygons[1].Edges.size(), 3);
        ASSERT_EQ(result.NewPolygons[1].Edges, list<unsigned int>({ 1, 4, 6 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonNormal)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check normal of polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        ASSERT_DOUBLE_EQ(normal[0], 0.0);
        ASSERT_DOUBLE_EQ(normal[1], 0.0);
        ASSERT_DOUBLE_EQ(normal[2], 1.0);
      }

      // check normal of polygon 3D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 1.0, 0.0, 0.0;
        polygonVertices.col(1)<< 0.0, 1.0, 0.0;
        polygonVertices.col(2)<< 0.0, 0.0, 1.0;

        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        ASSERT_DOUBLE_EQ(normal[0], 1.0 / sqrt(3));
        ASSERT_DOUBLE_EQ(normal[1], 1.0 / sqrt(3));
        ASSERT_DOUBLE_EQ(normal[2], 1.0 / sqrt(3));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonRotationMatrix)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check rotation matrix of polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Matrix3d rotationMatrix;
        Eigen::Vector3d translation;
        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        ASSERT_NO_THROW(geometryUtility.PolygonRotation(polygonVertices,
                                                        normal,
                                                        rotationMatrix,
                                                        translation));

        ASSERT_DOUBLE_EQ(translation[0], polygonVertices(0, 0));
        ASSERT_DOUBLE_EQ(translation[1], polygonVertices(1, 0));
        ASSERT_DOUBLE_EQ(translation[2], polygonVertices(2, 0));
        ASSERT_DOUBLE_EQ(rotationMatrix(0, 0), 1.0);
        ASSERT_DOUBLE_EQ(rotationMatrix(1, 1), 1.0);
        ASSERT_DOUBLE_EQ(rotationMatrix(2, 2), 1.0);
      }

      // check normal of polygon 3D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 1.0, 0.0, 0.0;
        polygonVertices.col(1)<< 0.0, 1.0, 0.0;
        polygonVertices.col(2)<< 0.0, 0.0, 1.0;

        Eigen::Matrix3d rotationMatrix;
        Eigen::Vector3d translation;
        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        ASSERT_NO_THROW(geometryUtility.PolygonRotation(polygonVertices,
                                                        normal,
                                                        rotationMatrix,
                                                        translation));

        ASSERT_DOUBLE_EQ(translation[0], polygonVertices(0, 0));
        ASSERT_DOUBLE_EQ(translation[1], polygonVertices(1, 0));
        ASSERT_DOUBLE_EQ(translation[2], polygonVertices(2, 0));
        ASSERT_DOUBLE_EQ(rotationMatrix(0, 0), -7.0710678118654724e-01);
        ASSERT_DOUBLE_EQ(rotationMatrix(1, 1), -4.0824829046386313e-01);
        ASSERT_DOUBLE_EQ(rotationMatrix(2, 2), 5.7735026918962651e-01);

        // export rotated polygons
        Eigen::MatrixXd rotatedPoints = geometryUtility.RotatePointsFrom3DTo2D(polygonVertices,
                                                                               rotationMatrix.transpose(),
                                                                               translation);
        Eigen::MatrixXd rotatedBackPoints = geometryUtility.RotatePointsFrom2DTo3D(rotatedPoints,
                                                                                   rotationMatrix,
                                                                                   translation);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, IntersectionSegmentPlane)
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
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_H
