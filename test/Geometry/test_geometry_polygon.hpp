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

        Eigen::VectorXd edgeLengths = geometryUtility.PolygonEdgeLengths(polygonVertices);
        ASSERT_EQ(edgeLengths, (Eigen::VectorXd(3) << 1.0,sqrt(2.0),1.0).finished());

        Eigen::MatrixXd edgeTangents = geometryUtility.PolygonEdgeTangents(polygonVertices);
        ASSERT_EQ(edgeTangents, (Eigen::MatrixXd(3, 3) << 1.0,-1.0,0.0, 0.0,1.0,-1.0, 0.0,0.0,0.0).finished());

        Eigen::MatrixXd edgeNormals = geometryUtility.PolygonEdgeNormals(polygonVertices);
        ASSERT_EQ(edgeNormals, (Eigen::MatrixXd(3, 3) << 0.0,1.0/sqrt(2),-1.0, -1.0,1.0/sqrt(2),0.0, 0.0,0.0,0.0).finished());
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

  TEST(TestGeometryUtilities, TestPolygonCentroid)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check centroid of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        double polygonArea = 1.0 / 2.0;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(barycenter[1], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroid[1], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

      // check area of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;
        double polygonArea = 1.0;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(barycenter[1], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroid[1], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

      // check area of generic triangle 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 3);
        polygonVertices.row(0) << -1.0, +5.0, +4.0;
        polygonVertices.row(1) << -2.0, -1.0, +5.0;
        double polygonArea = 1.850000000000000e+01;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 2.666666666666667e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 6.666666666666665e-01);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 2.666666666666667e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 6.666666666666665e-01);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

      // check area of generic quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 4);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;
        double polygonArea = 1.511000000000000e+01;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 3.100000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.850000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

      // check area of generic quadrilateral 2D with aligned points
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 3.000000000000000e+00, 5.700000000000000e+00, 5.000000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, 1.010638297872341e+00, -1.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;
        double polygonArea = 1.511000000000000e+01;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 3.400000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.401773049645390e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

      // check area of generic concave 2D polygon
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 2.500000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00, 2.000000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 3.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00, 3.000000000000000e+00;
        double polygonArea = 7.210000000000001e+00;

        Eigen::Vector3d barycenter = geometryUtility.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 2.816666666666666e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.900000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtility.PolygonCentroid(polygonVertices,
                                                                   polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 2.842903374942210e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.754923717059639e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);
      }

    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonDiameter)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check area of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        double diameter = geometryUtility.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, sqrt(2.0));
      }

      // check area of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        double diameter = geometryUtility.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, sqrt(2.0));
      }

      // check area of generic triangle 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 3);
        polygonVertices.row(0) << -1.0, +5.0, +4.0;
        polygonVertices.row(1) << -2.0, -1.0, +5.0;

        double diameter = geometryUtility.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 8.602325267042627e+00);
      }

      // check area of generic quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 4);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        double diameter = geometryUtility.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 7.300684899377592e+00);
      }

      // check area of generic quadrilateral 2D with aligned points
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 3.000000000000000e+00, 5.700000000000000e+00, 5.000000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, 1.010638297872341e+00, -1.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        double diameter = geometryUtility.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 7.300684899377592e+00);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonArea)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check area of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        vector<unsigned int> polygonTriangulation = geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtility.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 0.5);

        double areaWithTriangles = geometryUtility.PolygonArea(polygonVertices,
                                                               polygonTriangulation);
        ASSERT_DOUBLE_EQ(areaWithTriangles, 0.5);
      }

      // check area of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        vector<unsigned int> polygonTriangulation = geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtility.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.0);

        double areaWithTriangles = geometryUtility.PolygonArea(polygonVertices,
                                                               polygonTriangulation);
        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.0);
      }

      // check area of generic triangle 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 3);
        polygonVertices.row(0) << -1.0, +5.0, +4.0;
        polygonVertices.row(1) << -2.0, -1.0, +5.0;

        vector<unsigned int> polygonTriangulation = geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtility.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.850000000000000e+01);

        double areaWithTriangles = geometryUtility.PolygonArea(polygonVertices,
                                                               polygonTriangulation);
        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.850000000000000e+01);
      }

      // check area of generic quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 4);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        vector<unsigned int> polygonTriangulation = geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtility.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.511000000000000e+01);

        double areaWithTriangles = geometryUtility.PolygonArea(polygonVertices,
                                                               polygonTriangulation);
        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.511000000000000e+01);
      }

      // check area of generic quadrilateral 2D with aligned points
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 3.000000000000000e+00, 5.700000000000000e+00, 5.000000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, 1.010638297872341e+00, -1.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        vector<unsigned int> polygonTriangulation = geometryUtility.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtility.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.511000000000000e+01);

        double areaWithTriangles = geometryUtility.PolygonArea(polygonVertices,
                                                               polygonTriangulation);
        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.511000000000000e+01);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonType)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check triangle
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.PolygonType(polygonVertices),
                  Gedim::GeometryUtilities::PolygonTypes::Triangle);
      }

      // check quadrilateral polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.25, 0.25, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtility.PolygonType(polygonVertices),
                  Gedim::GeometryUtilities::PolygonTypes::Quadrilateral);
      }

      // check triangle with aligned edges polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.5, 0.5, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        vector<unsigned int> unalignedPoint = geometryUtility.UnalignedPoints(polygonVertices);

        Eigen::MatrixXd extraction = geometryUtility.ExtractPoints(polygonVertices,
                                                                   unalignedPoint);

        ASSERT_EQ(geometryUtility.PolygonType(extraction),
                  Gedim::GeometryUtilities::PolygonTypes::Triangle);
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

        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtility.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtility.PolygonRotationMatrix(polygonVertices,
                                                                               normal,
                                                                               translation);
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

        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtility.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtility.PolygonRotationMatrix(polygonVertices,
                                                                               normal,
                                                                               translation);


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

        Eigen::Vector3d normal = geometryUtility.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtility.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtility.PolygonRotationMatrix(polygonVertices,
                                                                               normal,
                                                                               translation);

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
