#ifndef __TEST_GEOMETRY_POINT_H
#define __TEST_GEOMETRY_POINT_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {
  
  TEST(TestGeometryUtilities, TestPointsAreCoincident)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);
      
      // check coincident points
      {
        Eigen::Vector3d firstPoint(1.0, 2.0, 3.0);
        Eigen::Vector3d secondPoint(1.0, 2.0, 3.0);
        ASSERT_TRUE(geometryUtility.PointsAreCoincident(firstPoint,
                                                        secondPoint));
      }
      
      // check not coincident points
      {
        Eigen::Vector3d firstPoint(1.0, 2.0, 3.0);
        Eigen::Vector3d secondPoint(2.0, 2.0, 3.0);
        ASSERT_FALSE(geometryUtility.PointsAreCoincident(firstPoint,
                                                         secondPoint));
      }
      
      // check Find Point In Points
      {
        Eigen::MatrixXd points(3, 5);
        points.col(0)<< 0.0, 0.0, 0.0;
        points.col(1)<< 1.0, 2.0, 3.0;
        points.col(2)<< 4.0, 5.0, 6.0;
        points.col(3)<< 1.0, 2.0, 3.0;
        points.col(4)<< 1.1, 1.2, 1.3;
        
        ASSERT_EQ(geometryUtility.FindPointInPoints(points,
                                                    Eigen::Vector3d(6.0, 6.0, 6.0)),
                  vector<unsigned int>({ }));
        ASSERT_EQ(geometryUtility.FindPointInPoints(points,
                                                    Eigen::Vector3d(0.0, 0.0, 0.0)),
                  vector<unsigned int>({ 0 }));
        ASSERT_EQ(geometryUtility.FindPointInPoints(points,
                                                    Eigen::Vector3d(1.0, 2.0, 3.0)),
                  vector<unsigned int>({ 1, 3 }));
        ASSERT_EQ(geometryUtility.FindPointInPoints(points,
                                                    Eigen::Vector3d(4.0, 5.0, 6.0)),
                  vector<unsigned int>({ 2 }));
        ASSERT_EQ(geometryUtility.FindPointInPoints(points,
                                                    Eigen::Vector3d(1.1, 1.2, 1.3)),
                  vector<unsigned int>({ 4 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
  
  TEST(TestGeometryUtilities, TestPointsAre2D)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);
      
      // check 2D points
      {
        Eigen::MatrixXd points;
        points.setZero(3,2);
        points.row(0).setConstant(1.0);
        points.row(1).setConstant(2.0);
        ASSERT_TRUE(geometryUtility.PointsAre2D(points));
      }
      
      // check 2D points
      {
        Eigen::MatrixXd points;
        points.setZero(3,2);
        points.row(0).setConstant(1.0);
        points.row(1).setConstant(2.0);
        points.row(2).setConstant(3.0);
        ASSERT_FALSE(geometryUtility.PointsAre2D(points));
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
  
  TEST(TestGeometryUtilities, TestIsPointOnLine)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check not in line
      {
        const Eigen::Vector3d point(0.5, -1.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtility.IsPointOnLine(point,
                                                   lineOrigin,
                                                   lineTangent,
                                                   lineTangentSquaredLength));
      }

      // check not in line
      {
        const Eigen::Vector3d point(0.5, 1.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtility.IsPointOnLine(point,
                                                   lineOrigin,
                                                   lineTangent,
                                                   lineTangentSquaredLength));
      }

      // check on segment line before origin
      {
        const Eigen::Vector3d point(-10.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
      }

      // check on segment line after end
      {
        const Eigen::Vector3d point(10.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
      }

      // check on segment origin
      {
        const Eigen::Vector3d point(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
      }

      // check on segment end
      {
        const Eigen::Vector3d point(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
      }

      // check inside segment
      {
        const Eigen::Vector3d point(0.5, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
      }

      // check not in line 3D
      {
        const Eigen::Vector3d point(0.0, -1.0, -7.0);
        const Eigen::Vector3d lineOrigin(1.0, 2.0, 3.0);
        const Eigen::Vector3d lineEnd(5.0, 7.0, 9.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.norm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtility.IsPointOnLine(point,
                                                   lineOrigin,
                                                   lineTangent,
                                                   lineTangentSquaredLength));
      }

      // check in line 3D
      {
        const Eigen::Vector3d point(3.0, 4.5, 6.0);
        const Eigen::Vector3d lineOrigin(1.0, 2.0, 3.0);
        const Eigen::Vector3d lineEnd(5.0, 7.0, 9.0);
        const Eigen::Vector3d lineTangent = geometryUtility.SegmentTangent(lineOrigin,
                                                                           lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtility.IsPointOnLine(point,
                                                  lineOrigin,
                                                  lineTangent,
                                                  lineTangentSquaredLength));
        ASSERT_TRUE(geometryUtility.Are1DValuesEqual(geometryUtility.PointLineCurvilinearCoordinate(point,
                                                                                                    lineOrigin,
                                                                                                    lineTangent,
                                                                                                    lineTangentSquaredLength),
                                                     0.5));
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
      
      // check 2D right
      {
        Eigen::Vector3d point(0.5, -1.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);
        
        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtility.PointSegmentPosition(point,
                                                       segmentOrigin,
                                                       segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::RightTheSegment);
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         0.5);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         0.5);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         -10.0);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         10.0);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         0.0);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         1.0);
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
        ASSERT_DOUBLE_EQ(geometryUtility.PointSegmentProjection(point,
                                                                segmentOrigin,
                                                                segmentEnd),
                         0.5);
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
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside);
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
        ASSERT_EQ(resultBorderEdge.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge);
        ASSERT_EQ(resultBorderEdge.BorderIndex, 1);
        
        // border vertex
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderVertex = geometryUtility.PointPolygonPosition(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                                       polygonVertices);
        ASSERT_EQ(resultBorderVertex.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
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
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside);
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
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
      }

      // check on vertex aligned edge
      {
        Eigen::Vector3d point(7.5000000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00);
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.row(0)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 3.7500000000000000e-01, 7.5000000000000000e-01, 1.0000000000000000e+00;
        polygonVertices.row(1)<< 1.0000000000000000e+00, 7.5000000000000000e-01, 3.7500000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonVertices.row(2)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtility.PointPolygonPosition(point,
                                                                                                           polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPointPointPlanePosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      const Eigen::Vector3d planeNormal = Eigen::Vector3d::Constant(1.0).normalized();
      const Eigen::Vector3d planeOrigin = Eigen::Vector3d(0.0, 0.0, 1.0);

      // check point on plane
      {
        const Eigen::Vector3d point = Eigen::Vector3d(5.0, -5.0, 1.0);
        ASSERT_EQ(geometryUtility.PointPlanePosition(point,
                                                     planeNormal,
                                                     planeOrigin),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::OnPlane);
      }

      // check curvilinear coordinate coincident
      {
        const Eigen::Vector3d point = Eigen::Vector3d(-1.0, -1.0, -2.0);
        ASSERT_EQ(geometryUtility.PointPlanePosition(point,
                                                     planeNormal,
                                                     planeOrigin),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::Negative);
      }

      // check curvilinear coordinate before
      {
        const Eigen::Vector3d point = Eigen::Vector3d(0.0, 1.0, 2.0);
        ASSERT_EQ(geometryUtility.PointPlanePosition(point,
                                                     planeNormal,
                                                     planeOrigin),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::Positive);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPointCirclePosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);
      
      // check point outside
      {
        Eigen::Vector3d point(6.0, 4.0, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtility.PointCirclePosition(point,
                                                                                                         circleCenter,
                                                                                                         circleRadius);
        ASSERT_EQ(result, Gedim::GeometryUtilities::PointCirclePositionResult::Outside);
      }
      
      // check point on border
      {
        Eigen::Vector3d point(0.0, 4.0, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtility.PointCirclePosition(point,
                                                                                                         circleCenter,
                                                                                                         circleRadius);
        ASSERT_EQ(result, Gedim::GeometryUtilities::PointCirclePositionResult::OnBorder);
      }
      
      // check point inside
      {
        Eigen::Vector3d point(0.1, 2.5, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtility.PointCirclePosition(point,
                                                                                                         circleCenter,
                                                                                                         circleRadius);
        ASSERT_EQ(result, Gedim::GeometryUtilities::PointCirclePositionResult::Inside);
      }
      
      // check generic points
      {
        Eigen::MatrixXd points(3, 4);
        points.col(0)<< 1.0, 3.0, 0.0;
        points.col(1)<< 3.0, 3.0 - 2.0 / sqrt(3.0), 0.0;
        points.col(2)<< 4.0 + 1.0 / 10.0, 3.0, 0.0;
        points.col(3)<< 3.0, 4.0, 0.0;
        Eigen::Vector3d circleCenter(3.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> result = geometryUtility.PointCirclePositions(points,
                                                                                                                  circleCenter,
                                                                                                                  circleRadius);
        ASSERT_EQ(result.size(), 4);
        ASSERT_EQ(result[0], Gedim::GeometryUtilities::PointCirclePositionResult::Outside);
        ASSERT_EQ(result[1], Gedim::GeometryUtilities::PointCirclePositionResult::Outside);
        ASSERT_EQ(result[2], Gedim::GeometryUtilities::PointCirclePositionResult::Outside);
        ASSERT_EQ(result[3], Gedim::GeometryUtilities::PointCirclePositionResult::OnBorder);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPoint_PointPolyhedronPosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtility.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                    1.0);
      const Eigen::Vector3d polyhedronBarycenter = geometryUtility.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFaceVertices = geometryUtility.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                    polyhedron.Faces);
      const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtility.PolyhedronFaceNormals(polyhedronFaceVertices);
      const vector<bool> polyhedronFaceNormalDirections = geometryUtility.PolyhedronFaceNormalDirections(polyhedronFaceVertices,
                                                                                                         polyhedronBarycenter,
                                                                                                         polyhedronFaceNormals);
      const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtility.PolyhedronFaceTranslations(polyhedronFaceVertices);
      const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtility.PolyhedronFaceRotationMatrices(polyhedronFaceVertices,
                                                                                                                    polyhedronFaceNormals,
                                                                                                                    polyhedronFaceTranslations);

      // check point outside
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtility.PointPolyhedronPosition(Eigen::Vector3d(1.2, 1.7, -15.0),
                                                    polyhedron.Faces,
                                                    polyhedronFaceVertices,
                                                    polyhedronFaceNormals,
                                                    polyhedronFaceNormalDirections,
                                                    polyhedronFaceTranslations,
                                                    polyhedronFaceRotationMatrices);
        ASSERT_EQ(result.Type,
                  Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Outside);
      }

      // check point inside
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtility.PointPolyhedronPosition(Eigen::Vector3d(0.5, 0.75, 0.25),
                                                    polyhedron.Faces,
                                                    polyhedronFaceVertices,
                                                    polyhedronFaceNormals,
                                                    polyhedronFaceNormalDirections,
                                                    polyhedronFaceTranslations,
                                                    polyhedronFaceRotationMatrices);
        ASSERT_EQ(result.Type,
                  Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::Inside);
      }

      // check point on face
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtility.PointPolyhedronPosition(Eigen::Vector3d(0.5, 0.75, 1.0),
                                                    polyhedron.Faces,
                                                    polyhedronFaceVertices,
                                                    polyhedronFaceNormals,
                                                    polyhedronFaceNormalDirections,
                                                    polyhedronFaceTranslations,
                                                    polyhedronFaceRotationMatrices);
        ASSERT_EQ(result.Type,
                  Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderFace);
        ASSERT_EQ(result.BorderIndex,
                  1);
      }

      // check point on edge
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtility.PointPolyhedronPosition(Eigen::Vector3d(0.0, 0.5, 1.0),
                                                    polyhedron.Faces,
                                                    polyhedronFaceVertices,
                                                    polyhedronFaceNormals,
                                                    polyhedronFaceNormalDirections,
                                                    polyhedronFaceTranslations,
                                                    polyhedronFaceRotationMatrices);
        ASSERT_EQ(result.Type,
                  Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderEdge);
        ASSERT_EQ(result.BorderIndex,
                  7);
      }

      // check point on vertex
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtility.PointPolyhedronPosition(Eigen::Vector3d(1.0, 1.0, 1.0),
                                                    polyhedron.Faces,
                                                    polyhedronFaceVertices,
                                                    polyhedronFaceNormals,
                                                    polyhedronFaceNormalDirections,
                                                    polyhedronFaceTranslations,
                                                    polyhedronFaceRotationMatrices);
        ASSERT_EQ(result.Type,
                  Gedim::GeometryUtilities::PointPolyhedronPositionResult::Types::BorderVertex);
        ASSERT_EQ(result.BorderIndex,
                  6);
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
