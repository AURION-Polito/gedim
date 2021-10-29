#ifndef __TEST_GEOMETRY_POLYGON_H
#define __TEST_GEOMETRY_POLYGON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting {

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

  TEST(TestGeometryUtilities, TestPolygonIsConvex)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check convex polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        ASSERT_TRUE(geometryUtility.PolygonIsConvex(polygonVertices));
      }

      // check concave polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.25, 0.25, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        ASSERT_FALSE(geometryUtility.PolygonIsConvex(polygonVertices));
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

      // check rotation matrix of polygon 3D
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
      }

      // check rotation matrix of other polygon 3D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 1.0;
        polygonVertices.col(1)<< 0.0, 1.0, 0.0;
        polygonVertices.col(2)<< 1.0, 0.0, 0.0;

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
        ASSERT_DOUBLE_EQ(rotationMatrix(0, 0), -5.551115123125783e-17);
        ASSERT_DOUBLE_EQ(rotationMatrix(1, 1), -4.0824829046386296e-01);
        ASSERT_DOUBLE_EQ(rotationMatrix(2, 2), -5.7735026918962562e-01);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonTriangulation)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check triangle triangulation
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d internalPoint(0.25, 0.25, 0.0);

        ASSERT_EQ(geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices), vector<unsigned int>({ 0, 1, 2 }));
        ASSERT_EQ(geometryUtility.PolygonTriangulationByInternalPoint(polygonVertices,
                                                                      internalPoint), vector<unsigned int>({ 3, 0, 1, 3, 1, 2, 3, 2, 0 }));
      }

      // check square triangulation
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d internalPoint(0.25, 0.25, 0.0);

        ASSERT_EQ(geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices), vector<unsigned int>({ 0, 1, 2, 0, 2, 3 }));
        ASSERT_EQ(geometryUtility.PolygonTriangulationByInternalPoint(polygonVertices,
                                                                      internalPoint), vector<unsigned int>({ 4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonCirclePosition)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check circle outside and no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonOutsideCircleNoIntersection);
      }

      // check Polygon Inside Circle with center outside no intersections
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 10.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection);
      }

      // Polygon Inside Circle with center outside one intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(sqrt(2.0), sqrt(2.0), 0.0);
        double circleRadius = 2.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonInsideCircleOneVertexIntersection);
      }

      // check circle inside polygon no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
        double circleRadius = 0.125;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::CircleInsidePolygonNoIntersection);
      }

      // check Polygon Inside Circle with center inside no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.25, 0.25, 0.0);
        double circleRadius = 10.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonInsideCircleNoIntersection);
      }

      // check Polygon inside Circle intersects only vertices
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
        double circleRadius = sqrt(2) / 2;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonInsideCircleIntersectionOnlyOnVertices);
      }

      // check Circle Outside Polygon one intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 2.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionOnVertex);
      }

      // check Circle outside Polygon tangent to edge
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, -1.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::PolygonOutsideCircleOneIntersectionTangentOnEdge);
      }

      // check Circle inside Polygon tangent to edge
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.125, 0.0);
        double circleRadius = 0.125;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::CircleInsidePolygonOneIntersectionTangentOnEdge);
      }

      // check Circle Outside Polygon Intersects With Multiple SubPolygons
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.55, 0.0);
        double circleRadius = 0.5;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::CirclePolygonMultipleIntersections);
      }

      // check Circle Inside Polygon Intersects With Multiple SubPolygons
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.5, 0.5, 0.0);
        double circleRadius = 0.5;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::CirclePolygonMultipleIntersections);
      }

      // check generic intersections with quadrilateral
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 1.0, 3.0, 0.0;
        polygonVertices.col(1)<< 3.0, 3.0 - 2.0 / sqrt(3.0), 0.0;
        polygonVertices.col(2)<< 4.0 + 1.0 / 10.0, 3.0, 0.0;
        polygonVertices.col(3)<< 3.0, 4.0, 0.0;
        Eigen::Vector3d circleCenter(3.0, 3.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtility.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                         circleCenter,
                                                                                                                                         circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtility.PointCirclePositions(polygonVertices,
                                                                                                                           circleCenter,
                                                                                                                           circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtility.PolygonCirclePosition(polygonVertices,
                                                                                                              circleCenter,
                                                                                                              circleRadius,
                                                                                                              vertexPositions,
                                                                                                              polygonCircleIntersections);
        ASSERT_EQ(position, Gedim::GeometryUtilities::PolygonCirclePositionTypes::CirclePolygonMultipleIntersections);}
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_POLYGON_H
