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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check coincident points
      {
        Eigen::Vector3d firstPoint(1.0, 2.0, 3.0);
        Eigen::Vector3d secondPoint(1.0, 2.0, 3.0);
        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(firstPoint,
                                                          secondPoint));
      }
      
      // check not coincident points
      {
        Eigen::Vector3d firstPoint(1.0, 2.0, 3.0);
        Eigen::Vector3d secondPoint(2.0, 2.0, 3.0);
        ASSERT_FALSE(geometryUtilities.PointsAreCoincident(firstPoint,
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
        
        ASSERT_EQ(geometryUtilities.FindPointInPoints(points,
                                                      Eigen::Vector3d(6.0, 6.0, 6.0)),
                  vector<unsigned int>({ }));
        ASSERT_EQ(geometryUtilities.FindPointInPoints(points,
                                                      Eigen::Vector3d(0.0, 0.0, 0.0)),
                  vector<unsigned int>({ 0 }));
        ASSERT_EQ(geometryUtilities.FindPointInPoints(points,
                                                      Eigen::Vector3d(1.0, 2.0, 3.0)),
                  vector<unsigned int>({ 1, 3 }));
        ASSERT_EQ(geometryUtilities.FindPointInPoints(points,
                                                      Eigen::Vector3d(4.0, 5.0, 6.0)),
                  vector<unsigned int>({ 2 }));
        ASSERT_EQ(geometryUtilities.FindPointInPoints(points,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check 2D points
      {
        Eigen::MatrixXd points;
        points.setZero(3,2);
        points.row(0).setConstant(1.0);
        points.row(1).setConstant(2.0);
        ASSERT_TRUE(geometryUtilities.PointsAre2D(points));
      }
      
      // check 2D points
      {
        Eigen::MatrixXd points;
        points.setZero(3,2);
        points.row(0).setConstant(1.0);
        points.row(1).setConstant(2.0);
        points.row(2).setConstant(3.0);
        ASSERT_FALSE(geometryUtilities.PointsAre2D(points));
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check curvilinear coordinate before
      {
        ASSERT_EQ(geometryUtilities.Compare1DValues(-0.5, 1.5),
                  Gedim::GeometryUtilities::CompareTypes::FirstBeforeSecond);
      }
      
      // check curvilinear coordinate coincident
      {
        ASSERT_EQ(geometryUtilities.Compare1DValues(0.5, 0.5 + geometryUtilitiesConfig.Tolerance / 2.0),
                  Gedim::GeometryUtilities::CompareTypes::Coincident);
      }
      
      // check curvilinear coordinate before
      {
        ASSERT_EQ(geometryUtilities.Compare1DValues(0.5, -1.5),
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check curvilinear coordinate origin
      {
        Eigen::Vector3d point(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);
        
        ASSERT_TRUE(abs(geometryUtilities.PointCurvilinearCoordinate(point,
                                                                     segmentOrigin,
                                                                     segmentEnd) - 0.0) < geometryUtilitiesConfig.Tolerance);
      }
      
      // check curvilinear coordinate end
      {
        Eigen::Vector3d point(10.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);
        
        ASSERT_TRUE(abs(geometryUtilities.PointCurvilinearCoordinate(point,
                                                                     segmentOrigin,
                                                                     segmentEnd) - 1.0) < geometryUtilitiesConfig.Tolerance);
      }
      
      // check curvilinear coordinate inside
      {
        Eigen::Vector3d point(5.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);
        
        ASSERT_TRUE(abs(geometryUtilities.PointCurvilinearCoordinate(point,
                                                                     segmentOrigin,
                                                                     segmentEnd) - 0.5) < geometryUtilitiesConfig.Tolerance);
      }
      
      // check curvilinear coordinate before origin
      {
        Eigen::Vector3d point(-5.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);
        
        ASSERT_TRUE(abs(geometryUtilities.PointCurvilinearCoordinate(point,
                                                                     segmentOrigin,
                                                                     segmentEnd) - -0.5) < geometryUtilitiesConfig.Tolerance);
      }
      
      // check curvilinear coordinate after end
      {
        Eigen::Vector3d point(15.0, 0.0, 4.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 4.0);
        Eigen::Vector3d segmentEnd(10.0, 0.0, 4.0);
        
        ASSERT_TRUE(abs(geometryUtilities.PointCurvilinearCoordinate(point,
                                                                     segmentOrigin,
                                                                     segmentEnd) - 1.5) < geometryUtilitiesConfig.Tolerance);
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check not in line
      {
        const Eigen::Vector3d point(0.5, -1.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(1.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtilities.IsPointOnLine(point,
                                                     lineOrigin,
                                                     lineTangent,
                                                     lineTangentSquaredLength));
      }

      // check not in line
      {
        const Eigen::Vector3d point(0.5, 1.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(1.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtilities.IsPointOnLine(point,
                                                     lineOrigin,
                                                     lineTangent,
                                                     lineTangentSquaredLength));
      }

      // check on segment line before origin
      {
        const Eigen::Vector3d point(-10.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(0.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
      }

      // check on segment line after end
      {
        const Eigen::Vector3d point(10.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(0.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
      }

      // check on segment origin
      {
        const Eigen::Vector3d point(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(0.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
      }

      // check on segment end
      {
        const Eigen::Vector3d point(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(0.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
      }

      // check inside segment
      {
        const Eigen::Vector3d point(0.5, 0.0, 0.0);
        const Eigen::Vector3d lineOrigin(0.0, 0.0, 0.0);
        const Eigen::Vector3d lineEnd(1.0, 0.0, 0.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = Eigen::Vector3d(lineTangent.y(), -lineTangent.x(), 0.0);

        ASSERT_DOUBLE_EQ(0.0, geometryUtilities.PointLineDistance(point,
                                                                  lineOrigin,
                                                                  normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
      }

      // check not in line 3D
      {
        const Eigen::Vector3d point(0.0, -1.0, -7.0);
        const Eigen::Vector3d lineOrigin(1.0, 2.0, 3.0);
        const Eigen::Vector3d lineEnd(5.0, 7.0, 9.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();
        const Eigen::Vector3d normalToLine = lineTangent.cross(point - lineOrigin).cross(lineTangent);

        ASSERT_DOUBLE_EQ(5.3803393896716170e+00, geometryUtilities.PointLineDistance(point,
                                                                                     lineOrigin,
                                                                                     normalToLine));

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_FALSE(geometryUtilities.IsPointOnLine(point,
                                                     lineOrigin,
                                                     lineTangent,
                                                     lineTangentSquaredLength));
      }

      // check in line 3D
      {
        const Eigen::Vector3d point(3.0, 4.5, 6.0);
        const Eigen::Vector3d lineOrigin(1.0, 2.0, 3.0);
        const Eigen::Vector3d lineEnd(5.0, 7.0, 9.0);
        const Eigen::Vector3d lineTangent = geometryUtilities.SegmentTangent(lineOrigin,
                                                                             lineEnd);
        const double lineTangentSquaredLength = lineTangent.squaredNorm();

        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_TRUE(geometryUtilities.IsPointOnLine(point,
                                                    lineOrigin,
                                                    lineTangent,
                                                    lineTangentSquaredLength));
        ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(geometryUtilities.PointLineCurvilinearCoordinate(point,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check 2D right
      {
        Eigen::Vector3d point(0.5, -1.0, 0.0);
        Eigen::Vector3d segmentOrigin(0.0, 0.0, 0.0);
        Eigen::Vector3d segmentEnd(   1.0, 0.0, 0.0);
        
        Gedim::GeometryUtilities::PointSegmentPositionTypes result;
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::RightTheSegment);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::LeftTheSegment);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineBeforeOrigin);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd),
                  Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentLineAfterEnd);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentOrigin);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::OnSegmentEnd);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
        ASSERT_EQ(geometryUtilities.PointSegmentPosition(point,
                                                         segmentOrigin,
                                                         segmentEnd), Gedim::GeometryUtilities::PointSegmentPositionTypes::InsideSegment);
        ASSERT_DOUBLE_EQ(geometryUtilities.PointSegmentProjection(point,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check outside
      {
        Eigen::Vector3d point(0.5, -1.0, 0.0);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        
        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtilities.PointPolygonPosition(point,
                                                                                                             polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside);
        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside);
        ASSERT_FALSE(geometryUtilities.IsPointInsidePolygon(point,
                                                            polygonVertices));
        ASSERT_FALSE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                       polygonVertices));
      }
      
      // check border
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        
        // border edge
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderEdge= geometryUtilities.PointPolygonPosition(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                                                                                      polygonVertices);
        ASSERT_EQ(resultBorderEdge.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge);
        ASSERT_EQ(resultBorderEdge.BorderIndex, 1);
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderEdge_ray = geometryUtilities.PointPolygonPosition_RayCasting(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                                                                                                      polygonVertices);
        ASSERT_EQ(resultBorderEdge_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderEdge);
        ASSERT_EQ(resultBorderEdge_ray.BorderIndex, 1);

        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                           polygonVertices));
        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(Eigen::Vector3d(0.5, 0.5, 0.0),
                                                                      polygonVertices));
        
        // border vertex
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderVertex = geometryUtilities.PointPolygonPosition(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                                         polygonVertices);
        ASSERT_EQ(resultBorderVertex.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
        ASSERT_EQ(resultBorderVertex.BorderIndex, 1);
        Gedim::GeometryUtilities::PointPolygonPositionResult resultBorderVertex_ray = geometryUtilities.PointPolygonPosition_RayCasting(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                                                                                        polygonVertices);
        ASSERT_EQ(resultBorderVertex_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
        ASSERT_EQ(resultBorderVertex_ray.BorderIndex, 1);

        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                           polygonVertices));
        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(Eigen::Vector3d(1.0, 0.0, 0.0),
                                                                      polygonVertices));
      }
      
      // check inside
      {
        Eigen::Vector3d point(0.50, 0.25, 0.0);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        
        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtilities.PointPolygonPosition(point,
                                                                                                             polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside);
        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside);


        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon(point,
                                                           polygonVertices));
        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                      polygonVertices));
      }

      // check inside concave
      {
        Eigen::Vector3d point(-0.5, 0.0, 0.0);
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        polygonVertices.col(3)<< -1.0, 0.0, 0.0;
        polygonVertices.col(4)<< 0.0, -1.0, 0.0;

        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside);


        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                      polygonVertices));
      }

      // check inside on rombo
      {
        Eigen::Vector3d point(-0.5, 0.0, 0.0);
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 1.0, 0.0, 0.0;
        polygonVertices.col(1)<< 0.0, 1.0, 0.0;
        polygonVertices.col(2)<< -1.0, 0.0, 0.0;
        polygonVertices.col(3)<< 0.0, -1.0, 0.0;

        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Inside);


        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                      polygonVertices));
      }
      
      // check on vertex
      {
        Eigen::Vector3d point(9.9999999999999956e-01, 2.0000000000000000e+00, 1.2412670766236366e-16);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.row(0)<< 9.9999999999999956e-01, 9.3749999999999956e-01, 1.0409363447223732e+00;
        polygonVertices.row(1)<< 2.0000000000000000e+00, 1.9234735079187608e+00, 1.8120621438385331e+00;
        polygonVertices.row(2)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        
        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtilities.PointPolygonPosition(point,
                                                                                                             polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);


        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon(point,
                                                           polygonVertices));
        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                      polygonVertices));
      }

      // check outside aligned on edge
      {
        Eigen::Vector3d point(-0.5, 0.0, 0.0);
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtilities.PointPolygonPosition(point,
                                                                                                             polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside);
        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::Outside);


        ASSERT_FALSE(geometryUtilities.IsPointInsidePolygon(point,
                                                            polygonVertices));
        ASSERT_FALSE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                       polygonVertices));
      }

      // check on vertex aligned edge
      {
        Eigen::Vector3d point(7.5000000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00);
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.row(0)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 3.7500000000000000e-01, 7.5000000000000000e-01, 1.0000000000000000e+00;
        polygonVertices.row(1)<< 1.0000000000000000e+00, 7.5000000000000000e-01, 3.7500000000000000e-01, 0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonVertices.row(2)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Gedim::GeometryUtilities::PointPolygonPositionResult result = geometryUtilities.PointPolygonPosition(point,
                                                                                                             polygonVertices);
        ASSERT_EQ(result.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);
        Gedim::GeometryUtilities::PointPolygonPositionResult result_ray = geometryUtilities.PointPolygonPosition_RayCasting(point,
                                                                                                                            polygonVertices);
        ASSERT_EQ(result_ray.Type, Gedim::GeometryUtilities::PointPolygonPositionResult::Types::BorderVertex);


        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon(point,
                                                           polygonVertices));
        ASSERT_TRUE(geometryUtilities.IsPointInsidePolygon_RayCasting(point,
                                                                      polygonVertices));
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      const Eigen::Vector3d planeNormal = Eigen::Vector3d::Constant(1.0).normalized();
      const Eigen::Vector3d planeOrigin = Eigen::Vector3d(0.0, 0.0, 1.0);
      const std::array<Eigen::Vector3d, 3> planePoints =
      {
        Eigen::Vector3d(0.0, 0.0, 1.0),
        Eigen::Vector3d(2.0, -1.0, 0.0),
        Eigen::Vector3d(3.0, -2.0, 0.0)
      };

      // check point on plane
      {
        const Eigen::Vector3d point = Eigen::Vector3d(5.0, -5.0, 1.0);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planeNormal,
                                                         planeOrigin)),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::OnPlane);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planePoints)),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::OnPlane);
      }

      // check curvilinear coordinate coincident
      {
        const Eigen::Vector3d point = Eigen::Vector3d(-1.0, -1.0, -2.0);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planeNormal,
                                                         planeOrigin)),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::Negative);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planePoints)),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::Negative);
      }

      // check curvilinear coordinate before
      {
        const Eigen::Vector3d point = Eigen::Vector3d(0.0, 1.0, 2.0);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planeNormal,
                                                         planeOrigin)),
                  Gedim::GeometryUtilities::PointPlanePositionTypes::Positive);
        ASSERT_EQ(geometryUtilities.PointPlanePosition(
                    geometryUtilities.PointPlaneDistance(point,
                                                         planePoints)),
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
      
      // check point outside
      {
        Eigen::Vector3d point(6.0, 4.0, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtilities.PointCirclePosition(point,
                                                                                                           circleCenter,
                                                                                                           circleRadius);
        ASSERT_EQ(result, Gedim::GeometryUtilities::PointCirclePositionResult::Outside);
      }
      
      // check point on border
      {
        Eigen::Vector3d point(0.0, 4.0, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtilities.PointCirclePosition(point,
                                                                                                           circleCenter,
                                                                                                           circleRadius);
        ASSERT_EQ(result, Gedim::GeometryUtilities::PointCirclePositionResult::OnBorder);
      }
      
      // check point inside
      {
        Eigen::Vector3d point(0.1, 2.5, 0.0);
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;
        
        Gedim::GeometryUtilities::PointCirclePositionResult result = geometryUtilities.PointCirclePosition(point,
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
        
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> result = geometryUtilities.PointCirclePositions(points,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      const Gedim::GeometryUtilities::Polyhedron polyhedron =  geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0),
                                                                                                      1.0);
      const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
      const vector<Eigen::Vector3d> polyhedronFaceNormals = geometryUtilities.PolyhedronFaceNormals(polyhedronFace3DVertices);
      const vector<bool> polyhedronFaceNormalDirections = geometryUtilities.PolyhedronFaceNormalDirections(polyhedronFace3DVertices,
                                                                                                           polyhedronBarycenter,
                                                                                                           polyhedronFaceNormals);
      const vector<Eigen::Vector3d> polyhedronFaceTranslations = geometryUtilities.PolyhedronFaceTranslations(polyhedronFace3DVertices);
      const vector<Eigen::Matrix3d> polyhedronFaceRotationMatrices = geometryUtilities.PolyhedronFaceRotationMatrices(polyhedronFace3DVertices,
                                                                                                                      polyhedronFaceNormals,
                                                                                                                      polyhedronFaceTranslations);

      const vector<Eigen::MatrixXd> polyhedronFace2DVertices = geometryUtilities.PolyhedronFaceRotatedVertices(polyhedronFace3DVertices,
                                                                                                               polyhedronFaceTranslations,
                                                                                                               polyhedronFaceRotationMatrices);
      // check point outside
      {
        Gedim::GeometryUtilities::PointPolyhedronPositionResult result =
            geometryUtilities.PointPolyhedronPosition(Eigen::Vector3d(1.2, 1.7, -15.0),
                                                      polyhedron.Faces,
                                                      polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
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
            geometryUtilities.PointPolyhedronPosition(Eigen::Vector3d(0.5, 0.75, 0.25),
                                                      polyhedron.Faces,
                                                      polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
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
            geometryUtilities.PointPolyhedronPosition(Eigen::Vector3d(0.5, 0.75, 1.0),
                                                      polyhedron.Faces,
                                                      polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
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
            geometryUtilities.PointPolyhedronPosition(Eigen::Vector3d(0.0, 0.5, 1.0),
                                                      polyhedron.Faces,
                                                      polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
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
            geometryUtilities.PointPolyhedronPosition(Eigen::Vector3d(1.0, 1.0, 1.0),
                                                      polyhedron.Faces,
                                                      polyhedronFace3DVertices,
                                                      polyhedronFace2DVertices,
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

  TEST(TestGeometryUtilities, TestPointsBoundingBox)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-08;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      Eigen::MatrixXd points(3, 5);
      points.col(0)<< 0.50, 0.2, 0.9;
      points.col(1)<< 0.40, 0.5, 0.3;
      points.col(2)<< 0.80, 0.7, 0.1;
      points.col(3)<< 0.45, 0.4, 0.9;
      points.col(4)<< 0.59, 0.2, 0.3;

      Eigen::MatrixXd expected_result(3, 2);
      expected_result.col(0)<< 0.40, 0.2, 0.1;
      expected_result.col(1)<< 0.80, 0.7, 0.9;

      ASSERT_EQ(expected_result,
                geometryUtilities.PointsBoundingBox(points));
      ASSERT_TRUE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(0.40, 0.5, 0.3),
                                                         expected_result));
      ASSERT_TRUE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(0.41, 0.45, 0.25),
                                                         expected_result));
      ASSERT_FALSE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(1.39, 0.45, 0.25),
                                                          expected_result));
      ASSERT_FALSE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(0.39, 1.45, 0.25),
                                                          expected_result));
      ASSERT_FALSE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(0.39, 0.45, 1.25),
                                                          expected_result));
      ASSERT_FALSE(geometryUtilities.IsPointInBoundingBox(Eigen::Vector3d(1.39, 1.45, 1.25),
                                                          expected_result));
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_POINT_H
