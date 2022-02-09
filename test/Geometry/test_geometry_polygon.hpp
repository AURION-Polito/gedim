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

        vector<unsigned int> polygonTriangulation = { 0, 1, 2 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 3 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<unsigned int> polygonTriangulation = { 0, 1, 2 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 2.666666666666667e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 6.666666666666665e-01);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 3 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<unsigned int> polygonTriangulation = { 0, 2, 3, 0, 3, 4, 0, 4, 5 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 5, 2, 3, 4, 2, 4, 5 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtility.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtility.PolygonCentroid(polygonTriangulationCentroids,
                                                                                polygonTriangulationAreas,
                                                                                polygonArea);

        ASSERT_DOUBLE_EQ(centroidWithTriangles[0], 2.842903374942210e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[1], 2.754923717059639e+00);
        ASSERT_DOUBLE_EQ(centroidWithTriangles[2], 0.0);
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

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

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

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

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

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

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

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

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

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtility.ExtractTriangulationPoints(polygonVertices,
                                                                                                        polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtility.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

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


  TEST(TestGeometryUtilities, TestPolygonDivisionByExternalPointAndEdge)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check triangle sub-division
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Matrix3d polygonEdgeTangents;
        polygonEdgeTangents.col(0)<<  1.0, 0.0, 0.0;
        polygonEdgeTangents.col(1)<<  -1.0, 1.0, 0.0;
        polygonEdgeTangents.col(2)<<  0.0, -1.0, 0.0;

        const Eigen::Vector3d circleCenter(0.5, -0.5, 0.0);
        const double circleRadius = sqrt(2.0) / 2.0;
        const unsigned int curvedEdgeIndex = 0;

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult result =
            geometryUtility.PolygonDivisionByCircle(polygonVertices,
                                                    polygonEdgeTangents,
                                                    circleCenter,
                                                    circleRadius,
                                                    curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 3)<< polygonVertices;
        expectedResult.Points.col(3)<< circleCenter;
        expectedResult.Points.col(4)<< Eigen::Vector3d(2.7639320225002101e-01,
                                                       1.7082039324993703e-01,
                                                       0.0000000000000000e+00);
        expectedResult.SubTriangles = { {3, 1, 2}, {3, 2, 0} };
        expectedResult.InternalTriangles = { {3, 1, 4}, {3, 4, 0} };
        expectedResult.SubPolygons = { {1, 2, 4}, {4, 2, 0} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }

      // check trapezioid sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.1, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 0.1, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.col(0)<<   9.0000000000000002e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<   0.0000000000000000e+00, -9.0000000000000002e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<<   1.0000000000000001e-01, -1.0000000000000001e-01, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.1;
        const unsigned int curvedEdgeIndex = 3;

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult result =
            geometryUtility.PolygonDivisionByCircle(polygonVertices,
                                                    polygonEdgeTangents,
                                                    circleCenter,
                                                    circleRadius,
                                                    curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.Points.col(4)<< circleCenter;

        expectedResult.SubTriangles = { {4, 1, 2} };
        expectedResult.InternalTriangles = { {4, 0, 3} };
        expectedResult.SubPolygons = { {0, 1, 2, 3} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }

      // check exagon sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 6);
        polygonVertices.col(0)<< 0.1, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 2.0 / 3.0, 1.0 / 3.0, 0.0;
        polygonVertices.col(3)<< 1.0 / 3.0, 2.0 / 3.0, 0.0;
        polygonVertices.col(4)<< 0.0, 1.0, 0.0;
        polygonVertices.col(5)<< 0.0, 0.1, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 6);
        polygonEdgeTangents.col(0)<<   9.0000000000000002e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  -3.3333333333333337e-01,  3.3333333333333331e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<  -3.3333333333333331e-01,  3.3333333333333331e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<<  -3.3333333333333331e-01,  3.3333333333333337e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(4)<<   0.0000000000000000e+00, -9.0000000000000002e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(5)<<   1.0000000000000001e-01, -1.0000000000000001e-01, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.1;
        const unsigned int curvedEdgeIndex = 5;

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult result =
            geometryUtility.PolygonDivisionByCircle(polygonVertices,
                                                    polygonEdgeTangents,
                                                    circleCenter,
                                                    circleRadius,
                                                    curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult expectedResult;
        expectedResult.Points.setZero(3, 9);
        expectedResult.Points.block(0, 0, 3, 6)<< polygonVertices;
        expectedResult.Points.col(6)<< circleCenter;
        expectedResult.Points.col(7)<< Eigen::Vector3d(8.9442719099991908e-02,
                                                       4.4721359549995954e-02,
                                                       0.0);
        expectedResult.Points.col(8)<< Eigen::Vector3d(4.4721359549995954e-02,
                                                       8.9442719099991908e-02,
                                                       0.0);

        expectedResult.SubTriangles = { {6, 1, 2}, {6, 2, 3}, {6, 3, 4} };
        expectedResult.InternalTriangles = { {6, 0, 7}, {6, 7, 8}, {6, 8, 5} };
        expectedResult.SubPolygons = { {0, 1, 2, 7}, {7, 2, 3, 8}, {8, 3, 4, 5} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonDivisionByAngleQuadrant)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      geometryUtilityConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check triangle no sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.1, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 0.1, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.col(0)<<   9.0000000000000002e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<   0.0000000000000000e+00, -9.0000000000000002e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<<   1.0000000000000001e-01, -1.0000000000000001e-01, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.01;
        const unsigned int curvedEdgeIndex = 3;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtility.PolygonDivisionByAngleQuadrant(polygonVertices,
                                                           polygonEdgeTangents,
                                                           circleCenter,
                                                           circleRadius,
                                                           curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 4);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.SubPolygons = { {3, 0, 1, 2} };
        expectedResult.SubPolygonTypes =
        { Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
      }

      // check triangle sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 6);
        polygonVertices.col(0)<< 0.1, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 0.1, 0.0;
        polygonVertices.col(4)<< 2.3542486889354099e-02, 7.6457513110645914e-02, 0.0;
        polygonVertices.col(5)<< 7.6457513110645914e-02, 2.3542486889354099e-02, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 6);
        polygonEdgeTangents.col(0)<<  9.0000000000000002e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<< -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<  0.0000000000000000e+00, -9.0000000000000002e-01, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<<  2.3542486889354099e-02, -2.3542486889354092e-02, 0.0000000000000000e+00;
        polygonEdgeTangents.col(4)<<  5.2915026221291815e-02, -5.2915026221291815e-02, 0.0000000000000000e+00;
        polygonEdgeTangents.col(5)<<  2.3542486889354092e-02, -2.3542486889354099e-02, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.08;
        const unsigned int curvedEdgeIndex = 4;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtility.PolygonDivisionByAngleQuadrant(polygonVertices,
                                                           polygonEdgeTangents,
                                                           circleCenter,
                                                           circleRadius,
                                                           curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 8);
        expectedResult.Points.block(0, 0, 3, 6)<< polygonVertices;
        expectedResult.Points.col(6)<< Eigen::Vector3d(2.3542486889354097e-01,
                                                       7.6457513110645903e-01,
                                                       0.0000000000000000e+00);
        expectedResult.Points.col(7)<< Eigen::Vector3d(7.6457513110645903e-01,
                                                       2.3542486889354097e-01,
                                                       0.0000000000000000e+00);
        expectedResult.SubPolygons = { {4, 6, 2, 3}, {5, 0, 1, 7}, {4, 5, 7, 6} };
        expectedResult.SubPolygonTypes =
        { Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalOrigin,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalEnd,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal
        };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
      }

      // check triangle other sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.col(0)<< -1.0, 0.0, 0.0;
        polygonVertices.col(1)<< -1.7320508075687857e-02, 0.0, 0.0;
        polygonVertices.col(2)<< 1.7320508075687746e-02, 0.0, 0.0;
        polygonVertices.col(3)<< 1.0, 0.0, 0.0;
        polygonVertices.col(4)<< 0.0, 1.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 5);
        polygonEdgeTangents.col(0)<<  9.8267949192431214e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  3.4641016151375603e-02,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<  9.8267949192431225e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<< -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(4)<< -1.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, -0.01, 0.0);
        const double circleRadius = 0.02;
        const unsigned int curvedEdgeIndex = 1;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtility.PolygonDivisionByAngleQuadrant(polygonVertices,
                                                           polygonEdgeTangents,
                                                           circleCenter,
                                                           circleRadius,
                                                           curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 7);
        expectedResult.Points.block(0, 0, 3, 5)<< polygonVertices;
        expectedResult.Points.col(5)<< Eigen::Vector3d(-6.4031434217770455e-01,
                                                       3.5968565782229545e-01,
                                                       0.0000000000000000e+00);
        expectedResult.Points.col(6)<< Eigen::Vector3d(6.4031434217770300e-01,
                                                       3.5968565782229700e-01,
                                                       0.0000000000000000e+00);
        expectedResult.SubPolygons = { {1, 5, 0}, {2, 3, 6}, {1, 2, 6, 4, 5} };
        expectedResult.SubPolygonTypes =
        { Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalOrigin,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalEnd,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal
        };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
      }

      // check triangle sub-division with only origin
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< -1.0, 0.0, 0.0;
        polygonVertices.col(1)<< -1.7320508075687857e-02, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 0.01, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.col(0)<<  9.8267949192431214e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  1.7320508075687857e-02,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<<  0.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<< -1.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, -0.01, 0.0);
        const double circleRadius = 0.02;
        const unsigned int curvedEdgeIndex = 1;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtility.PolygonDivisionByAngleQuadrant(polygonVertices,
                                                           polygonEdgeTangents,
                                                           circleCenter,
                                                           circleRadius,
                                                           curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.Points.col(4)<< Eigen::Vector3d(-6.4031434217770455e-01,
                                                       3.5968565782229545e-01,
                                                       0.0000000000000000e+00);
        expectedResult.SubPolygons = { {1, 4, 0}, {1, 2, 3, 4} };
        expectedResult.SubPolygonTypes =
        { Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalOrigin,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal,
        };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
      }

      // check triangle sub-division with only end
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.01, 0.0;
        polygonVertices.col(1)<< 1.7320508075687746e-02, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 0.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.col(0)<<  1.7320508075687746e-02,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<<  9.8267949192431225e-01,  0.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(2)<< -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        polygonEdgeTangents.col(3)<<  0.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, -0.01, 0.0);
        const double circleRadius = 0.02;
        const unsigned int curvedEdgeIndex = 0;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtility.PolygonDivisionByAngleQuadrant(polygonVertices,
                                                           polygonEdgeTangents,
                                                           circleCenter,
                                                           circleRadius,
                                                           curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.Points.col(4)<< Eigen::Vector3d(6.4031434217770300e-01,
                                                       3.5968565782229700e-01,
                                                       0.0000000000000000e+00);
        expectedResult.SubPolygons = { {1, 2, 4}, {0, 1, 4, 3} };
        expectedResult.SubPolygonTypes =
        {
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::ExternalEnd,
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal
        };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
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
