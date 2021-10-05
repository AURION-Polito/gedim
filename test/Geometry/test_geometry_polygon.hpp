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

  TEST(TestGeometryUtilities, TestConvexHull)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check simple convex hull
      {
        Eigen::MatrixXd points(3, 4);
        points.col(0)<< 1.0, 1.0, 0.0;
        points.col(1)<< 0.0, 0.0, 0.0;
        points.col(2)<< 0.0, 1.0, 0.0;
        points.col(3)<< 1.0, 0.0, 0.0;

        vector<unsigned int> convexHull = geometryUtility.ConvexHull(points);
        ASSERT_EQ(convexHull, vector<unsigned int>({ 1, 3, 0, 2 }));
      }

      // check complex convex hull
      {
        Eigen::MatrixXd points(3, 8);
        points.col(0)<< 20.0, 0.0, 0.0;
        points.col(1)<< 30.0, 60.0, 0.0;
        points.col(2)<< 50.0, 40.0, 0.0;
        points.col(3)<< 70.0, 30.0, 0.0;
        points.col(4)<< 55.0, 20.0, 0.0;
        points.col(5)<< 50.0, 10.0, 0.0;
        points.col(6)<< 0.0, 30.0, 0.0;
        points.col(7)<< 15.0, 25.0, 0.0;

        vector<unsigned int> convexHull = geometryUtility.ConvexHull(points);
        ASSERT_EQ(convexHull, vector<unsigned int>({ 6, 0, 5, 3, 1 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_GEOMETRY_POLYGON_H
