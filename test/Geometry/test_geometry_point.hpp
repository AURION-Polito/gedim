#ifndef __TEST_GEOMETRY_POINT_H
#define __TEST_GEOMETRY_POINT_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

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

  TEST(TestGeometryUtilities, TestPointSegmentPosition)
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
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Eigen::Vector3d(0.5, 0.0, 0.0));
      }

      // check left
      {
        Eigen::Vector3d point(0.5, 1.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment);
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Eigen::Vector3d(0.5, 0.0, 0.0));
      }

      // check on segment line before origin
      {
        Eigen::Vector3d point(-10.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  point);
      }

      // check on segment line after end
      {
        Eigen::Vector3d point(10.0, 0.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  point);
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
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  point);
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
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  point);
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
        ASSERT_EQ(geometryUtility.PointSegmentProjection(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  point);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPointPolygonPosition)
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
}

#endif // __TEST_GEOMETRY_POINT_H
