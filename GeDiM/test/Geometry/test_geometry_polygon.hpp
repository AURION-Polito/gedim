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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check normal of polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(polygonVertices);
        ASSERT_DOUBLE_EQ(normal[0], 0.0);
        ASSERT_DOUBLE_EQ(normal[1], 0.0);
        ASSERT_DOUBLE_EQ(normal[2], 1.0);

        Eigen::VectorXd edgeLengths = geometryUtilities.PolygonEdgeLengths(polygonVertices);
        ASSERT_EQ(edgeLengths, (Eigen::VectorXd(3) << 1.0,sqrt(2.0),1.0).finished());

        Eigen::MatrixXd edgeTangents = geometryUtilities.PolygonEdgeTangents(polygonVertices);
        ASSERT_EQ(edgeTangents, (Eigen::MatrixXd(3, 3) << 1.0,-1.0,0.0, 0.0,1.0,-1.0, 0.0,0.0,0.0).finished());

        Eigen::MatrixXd edgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);
        ASSERT_EQ(edgeNormals, (Eigen::MatrixXd(3, 3) << 0.0,1.0/sqrt(2),-1.0, -1.0,1.0/sqrt(2),0.0, 0.0,0.0,0.0).finished());
      }

      // check normal of polygon 3D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 1.0, 0.0, 0.0;
        polygonVertices.col(1)<< 0.0, 1.0, 0.0;
        polygonVertices.col(2)<< 0.0, 0.0, 1.0;

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(polygonVertices);
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check centroid of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        double polygonArea = 1.0 / 2.0;

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(barycenter[1], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroid[1], 1.0 / 3.0);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 1, 2 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(barycenter[1], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroid[1], 1.0 / 2.0);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 3 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 2.666666666666667e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 6.666666666666665e-01);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 2.666666666666667e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 6.666666666666665e-01);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 1, 2 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 3.100000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.850000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 3 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 3.400000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.401773049645390e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 3.338451356717406e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.617008603573792e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 2, 3, 0, 3, 4, 0, 4, 5 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

        Eigen::Vector3d barycenter = geometryUtilities.PolygonBarycenter(polygonVertices);
        ASSERT_DOUBLE_EQ(barycenter[0], 2.816666666666666e+00);
        ASSERT_DOUBLE_EQ(barycenter[1], 2.900000000000000e+00);
        ASSERT_DOUBLE_EQ(barycenter[2], 0.0);

        Eigen::Vector3d centroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                     polygonArea);
        ASSERT_DOUBLE_EQ(centroid[0], 2.842903374942210e+00);
        ASSERT_DOUBLE_EQ(centroid[1], 2.754923717059639e+00);
        ASSERT_DOUBLE_EQ(centroid[2], 0.0);

        vector<unsigned int> polygonTriangulation = { 0, 1, 2, 0, 2, 5, 2, 3, 4, 2, 4, 5 };

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::MatrixXd polygonTriangulationCentroids(3, polygonTriangulationPoints.size());
        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
        {
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);
          polygonTriangulationCentroids.col(t) = geometryUtilities.PolygonBarycenter(polygonTriangulationPoints[t]);
        }

        Eigen::Vector3d centroidWithTriangles = geometryUtilities.PolygonCentroid(polygonTriangulationCentroids,
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

  TEST(TestGeometryUtilities, TestPolygonInRadius)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check in radius of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        const double polygonArea = geometryUtilities.PolygonArea(polygonVertices);
        const Eigen::Vector3d polygonCentroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                                  polygonArea);
        const Eigen::MatrixXd polygonEdgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);

        ASSERT_DOUBLE_EQ(sqrt(1.0 / 18.0),
                         geometryUtilities.PolygonInRadius(polygonVertices,
                                                           polygonCentroid,
                                                           polygonEdgeNormals));
      }

      // check in radius of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;
        const double polygonArea = geometryUtilities.PolygonArea(polygonVertices);
        const Eigen::Vector3d polygonCentroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                                  polygonArea);
        const Eigen::MatrixXd polygonEdgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);

        ASSERT_DOUBLE_EQ(0.5,
                         geometryUtilities.PolygonInRadius(polygonVertices,
                                                           polygonCentroid,
                                                           polygonEdgeNormals));
      }

      // check in radius of reference triangle 2D with aligned edges
      {
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 0.5, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 0.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;
        polygonVertices.col(4)<< 0.0, 0.5, 0.0;
        const double polygonArea = 0.5;
        const Eigen::Vector3d polygonCentroid(1.0 / 3.0, 1.0 / 3.0, 0.0);
        const Eigen::MatrixXd polygonEdgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);

        ASSERT_DOUBLE_EQ(sqrt(1.0 / 18.0),
                         geometryUtilities.PolygonInRadius(polygonVertices,
                                                           polygonCentroid,
                                                           polygonEdgeNormals));
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check diameter of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        double diameter = geometryUtilities.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, sqrt(2.0));
      }

      // check diameter of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        double diameter = geometryUtilities.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, sqrt(2.0));
      }

      // check diameter of generic triangle 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 3);
        polygonVertices.row(0) << -1.0, +5.0, +4.0;
        polygonVertices.row(1) << -2.0, -1.0, +5.0;

        double diameter = geometryUtilities.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 8.602325267042627e+00);
      }

      // check diameter of generic quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 4);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        double diameter = geometryUtilities.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 7.300684899377592e+00);
      }

      // check diameter of generic quadrilateral 2D with aligned points
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 3.000000000000000e+00, 5.700000000000000e+00, 5.000000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, 1.010638297872341e+00, -1.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        double diameter = geometryUtilities.PolygonDiameter(polygonVertices);
        ASSERT_DOUBLE_EQ(diameter, 7.300684899377592e+00);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolygonAspectRatio)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check aspect ratio of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        const double polygonArea = geometryUtilities.PolygonArea(polygonVertices);
        const Eigen::Vector3d polygonCentroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                                  polygonArea);
        const Eigen::MatrixXd polygonEdgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);
        const double polygonInRadius = geometryUtilities.PolygonInRadius(polygonVertices,
                                                                         polygonCentroid,
                                                                         polygonEdgeNormals);
        const double polygonDiameter = geometryUtilities.PolygonDiameter(polygonVertices);

        ASSERT_DOUBLE_EQ(sqrt(2.0) / (2.0 * sqrt(1.0 / 18.0)),
                         geometryUtilities.PolygonAspectRatio(polygonDiameter,
                                                              polygonInRadius));
      }

      // check in radius of reference quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 1.0, 1.0, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;
        const double polygonArea = geometryUtilities.PolygonArea(polygonVertices);
        const Eigen::Vector3d polygonCentroid = geometryUtilities.PolygonCentroid(polygonVertices,
                                                                                  polygonArea);
        const Eigen::MatrixXd polygonEdgeNormals = geometryUtilities.PolygonEdgeNormals(polygonVertices);
        const double polygonInRadius = geometryUtilities.PolygonInRadius(polygonVertices,
                                                                         polygonCentroid,
                                                                         polygonEdgeNormals);
        const double polygonDiameter = geometryUtilities.PolygonDiameter(polygonVertices);

        ASSERT_DOUBLE_EQ(sqrt(2.0),
                         geometryUtilities.PolygonAspectRatio(polygonDiameter,
                                                              polygonInRadius));
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check area of reference triangle 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        vector<unsigned int> polygonTriangulation = geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtilities.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 0.5);

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);

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

        vector<unsigned int> polygonTriangulation = geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtilities.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.0);

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.0);
      }

      // check area of generic triangle 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 3);
        polygonVertices.row(0) << -1.0, +5.0, +4.0;
        polygonVertices.row(1) << -2.0, -1.0, +5.0;

        vector<unsigned int> polygonTriangulation = geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtilities.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.850000000000000e+01);

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.850000000000000e+01);
      }

      // check area of generic quadrilateral 2D
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 4);
        polygonVertices.row(0) << 1.000000000000000e+00, 5.700000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, -1.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        vector<unsigned int> polygonTriangulation = geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtilities.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.511000000000000e+01);

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);

        double areaWithTriangles = polygonTriangulationAreas.sum();

        ASSERT_DOUBLE_EQ(areaWithTriangles, 1.511000000000000e+01);
      }

      // check area of generic quadrilateral 2D with aligned points
      {
        Eigen::MatrixXd polygonVertices;
        polygonVertices.setZero(3, 6);
        polygonVertices.row(0) << 1.000000000000000e+00, 3.000000000000000e+00, 5.700000000000000e+00, 5.000000000000000e+00, 4.300000000000000e+00, 1.400000000000000e+00;
        polygonVertices.row(1) << 2.500000000000000e+00, 1.010638297872341e+00, -1.000000000000000e+00, 2.000000000000000e+00, 5.000000000000000e+00, 4.900000000000000e+00;

        vector<unsigned int> polygonTriangulation = geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices);

        double area = geometryUtilities.PolygonArea(polygonVertices);
        ASSERT_DOUBLE_EQ(area, 1.511000000000000e+01);

        vector<Eigen::Matrix3d> polygonTriangulationPoints = geometryUtilities.ExtractTriangulationPoints(polygonVertices,
                                                                                                          polygonTriangulation);

        Eigen::VectorXd polygonTriangulationAreas(polygonTriangulationPoints.size());
        for (unsigned int t = 0; t < polygonTriangulationPoints.size(); t++)
          polygonTriangulationAreas[t] = geometryUtilities.PolygonArea(polygonTriangulationPoints[t]);

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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check triangle
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtilities.PolygonType(polygonVertices),
                  Gedim::GeometryUtilities::PolygonTypes::Triangle);
      }

      // check quadrilateral polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.25, 0.25, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        ASSERT_EQ(geometryUtilities.PolygonType(polygonVertices),
                  Gedim::GeometryUtilities::PolygonTypes::Quadrilateral);
      }

      // check triangle with aligned edges polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.5, 0.5, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        vector<unsigned int> unalignedPoint = geometryUtilities.UnalignedPoints(polygonVertices);

        Eigen::MatrixXd extraction = geometryUtilities.ExtractPoints(polygonVertices,
                                                                     unalignedPoint);

        ASSERT_EQ(geometryUtilities.PolygonType(extraction),
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check convex polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(polygonVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(polygonVertices,
                                                                                   convexHull);

        ASSERT_TRUE(geometryUtilities.PolygonIsConvex(polygonVertices,
                                                      convexHullVertices));
      }

      // check concave polygon 2D
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.25, 0.25, 0.0;
        polygonVertices.col(3)<< 0.0, 1.0, 0.0;

        const vector<unsigned int> convexHull = geometryUtilities.ConvexHull(polygonVertices);
        const Eigen::MatrixXd convexHullVertices = geometryUtilities.ExtractPoints(polygonVertices,
                                                                                   convexHull);

        ASSERT_FALSE(geometryUtilities.PolygonIsConvex(polygonVertices,
                                                       convexHullVertices));
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check rotation matrix of polygon 2D
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtilities.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtilities.PolygonRotationMatrix(polygonVertices,
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

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtilities.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtilities.PolygonRotationMatrix(polygonVertices,
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

        Eigen::Vector3d normal = geometryUtilities.PolygonNormal(polygonVertices);
        Eigen::Vector3d translation = geometryUtilities.PolygonTranslation(polygonVertices);
        Eigen::Matrix3d rotationMatrix = geometryUtilities.PolygonRotationMatrix(polygonVertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check triangle triangulation
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;

        Eigen::Vector3d internalPoint(0.25, 0.25, 0.0);

        ASSERT_EQ(geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices), vector<unsigned int>({ 0, 1, 2 }));
        ASSERT_EQ(geometryUtilities.PolygonTriangulationByInternalPoint(polygonVertices,
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

        ASSERT_EQ(geometryUtilities.PolygonTriangulationByFirstVertex(polygonVertices), vector<unsigned int>({ 0, 1, 2, 0, 2, 3 }));
        ASSERT_EQ(geometryUtilities.PolygonTriangulationByInternalPoint(polygonVertices,
                                                                        internalPoint), vector<unsigned int>({ 4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestCircleDivisionByPolygon)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

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

        const Eigen::Vector3d circleCenter(-0.5, -0.5, 0.0);
        const double circleRadius = sqrt(10.0) / 2.0;
        const unsigned int curvedEdgeIndex = 1;

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult result =
            geometryUtilities.CircleDivisionByPolygon(polygonVertices,
                                                      polygonEdgeTangents,
                                                      circleCenter,
                                                      circleRadius,
                                                      curvedEdgeIndex);

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 3)<< polygonVertices;
        expectedResult.Points.col(3)<< circleCenter;
        expectedResult.Points.col(4)<< Eigen::Vector3d(6.1803398874989490e-01, 6.1803398874989490e-01, 0.0);

        expectedResult.SubTriangles = { {3, 4, 2}, {3, 1, 4} };
        expectedResult.InternalTriangles = { {3, 0, 2}, {3, 1, 0} };
        expectedResult.SubPolygons = { {2, 0, 4}, {4, 0, 1} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }

      // check trapezioid sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.col(0)<< 0.0, 1.0, 0.0;
        polygonVertices.col(1)<< -1.0 / 3.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, -1.0 / 3.0, 0.0;
        polygonVertices.col(3)<< 1.0, 0.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.row(0)<< -3.3333333333333331e-01,  3.3333333333333331e-01,  1.0000000000000000e+00, -1.0000000000000000e+00;
        polygonEdgeTangents.row(1)<< -1.0000000000000000e+00, -3.3333333333333331e-01,  3.3333333333333331e-01,  1.0000000000000000e+00;
        polygonEdgeTangents.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(-0.5, -0.5, 0.0);
        const double circleRadius = sqrt(10.0) / 2.0;
        const unsigned int curvedEdgeIndex = 3;

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult result =
            geometryUtilities.CircleDivisionByPolygon(polygonVertices,
                                                      polygonEdgeTangents,
                                                      circleCenter,
                                                      circleRadius,
                                                      curvedEdgeIndex);

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.Points.col(4)<< circleCenter;

        expectedResult.SubTriangles = { {4, 3, 0} };
        expectedResult.InternalTriangles = { {4, 2, 1} };
        expectedResult.SubPolygons = { {0, 1, 2, 3} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }

      // check pentagon sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 5);
        polygonVertices.col(0)<< 0.0, 1.0, 0.0;
        polygonVertices.col(1)<< -2.0 / 9.0, -1.0 / 9.0, 0.0;
        polygonVertices.col(2)<< -1.0 / 9.0, -2.0 / 9.0, 0.0;
        polygonVertices.col(3)<< 0.0, -1.0 / 3.0, 0.0;
        polygonVertices.col(4)<< 1.0, 0.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 5);
        polygonEdgeTangents.row(0)<< -2.2222222222222221e-01,  1.1111111111111110e-01,  1.1111111111111110e-01, 1.0000000000000000e+00, -1.0000000000000000e+00;
        polygonEdgeTangents.row(1)<< -1.1111111111111112e+00, -1.1111111111111110e-01, -1.1111111111111110e-01, 3.3333333333333331e-01,  1.0000000000000000e+00;
        polygonEdgeTangents.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(-0.5, -0.5, 0.0);
        const double circleRadius = sqrt(10.0) / 2.0;
        const unsigned int curvedEdgeIndex = 4;

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult result =
            geometryUtilities.CircleDivisionByPolygon(polygonVertices,
                                                      polygonEdgeTangents,
                                                      circleCenter,
                                                      circleRadius,
                                                      curvedEdgeIndex);

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult expectedResult;
        expectedResult.Points.setZero(3, 8);
        expectedResult.Points.block(0, 0, 3, 5)<< polygonVertices;
        expectedResult.Points.col(5)<< circleCenter;
        expectedResult.Points.col(6)<< Eigen::Vector3d(4.1901827761725985e-01, 7.8662558866416377e-01, 0.0);
        expectedResult.Points.col(7)<< Eigen::Vector3d(7.8662558866416377e-01, 4.1901827761725985e-01, 0.0);

        expectedResult.SubTriangles = { {5, 6, 0}, {5, 7, 6}, {5, 4, 7} };
        expectedResult.InternalTriangles = { {5, 1, 0}, {5, 2, 1}, {5, 3, 2} };
        expectedResult.SubPolygons = { {0, 1, 6}, {6, 1, 2, 7}, {7, 2, 3, 4} };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubTriangles, expectedResult.SubTriangles);
        ASSERT_EQ(result.InternalTriangles, expectedResult.InternalTriangles);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
      }

      // check exagon sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 6);
        polygonVertices.col(0)<< 0.0, 1.0, 0.0;
        polygonVertices.col(1)<< -1.0 / 3.0, 0.0, 0.0;
        polygonVertices.col(2)<< -2.0 / 9.0, -1.0 / 9.0, 0.0;
        polygonVertices.col(3)<< -1.0 / 9.0, -2.0 / 9.0, 0.0;
        polygonVertices.col(4)<< 0.0, -1.0 / 3.0, 0.0;
        polygonVertices.col(5)<< 1.0, 0.0, 0.0;

        Eigen::MatrixXd polygonEdgeTangents(3, 6);
        polygonEdgeTangents.row(0)<< -3.3333333333333331e-01,  1.1111111111111110e-01,  1.1111111111111110e-01,  1.1111111111111110e-01, 1.0000000000000000e+00, -1.0000000000000000e+00;
        polygonEdgeTangents.row(1)<< -1.0000000000000000e+00, -1.1111111111111110e-01, -1.1111111111111110e-01, -1.1111111111111110e-01, 3.3333333333333331e-01,  1.0000000000000000e+00;
        polygonEdgeTangents.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(-0.5, -0.5, 0.0);
        const double circleRadius = sqrt(10.0) / 2.0;
        const unsigned int curvedEdgeIndex = 5;

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult result =
            geometryUtilities.CircleDivisionByPolygon(polygonVertices,
                                                      polygonEdgeTangents,
                                                      circleCenter,
                                                      circleRadius,
                                                      curvedEdgeIndex);

        Gedim::GeometryUtilities::CircleDivisionByPolygonResult expectedResult;
        expectedResult.Points.setZero(3, 9);
        expectedResult.Points.block(0, 0, 3, 6)<< polygonVertices;
        expectedResult.Points.col(6)<< circleCenter;
        expectedResult.Points.col(7)<< Eigen::Vector3d(4.1901827761725985e-01, 7.8662558866416377e-01, 0.0);
        expectedResult.Points.col(8)<< Eigen::Vector3d(7.8662558866416377e-01, 4.1901827761725985e-01, 0.0);

        expectedResult.SubTriangles = { {6, 7, 0}, {6, 8, 7}, {6, 5, 8} };
        expectedResult.InternalTriangles = { {6, 2, 1}, {6, 3, 2}, {6, 4, 3} };
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

  TEST(TestGeometryUtilities, TestPolygonDivisionByCircle)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

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
            geometryUtilities.PolygonDivisionByCircle(polygonVertices,
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
            geometryUtilities.PolygonDivisionByCircle(polygonVertices,
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
            geometryUtilities.PolygonDivisionByCircle(polygonVertices,
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

      // check single aligned triangle
      {
        Eigen::MatrixXd polygonVertices(3, 3);
        polygonVertices.row(0)<< 0.0000000000000000e+00, 1.6808892342696363e-03, 0.0000000000000000e+00;
        polygonVertices.row(1)<< 1.0000000000000001e-01, 9.9985872058917000e-02, 1.0166676129318664e-01;
        polygonVertices.row(2)<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Eigen::MatrixXd polygonEdgeTangents(3, 3);
        polygonEdgeTangents.col(0)<<  1.6808892342696363e-03, -1.6808892342696363e-03,  0.0000000000000000e+00;
        polygonEdgeTangents.col(1)<< -1.4127941083005857e-05,  1.6808892342696363e-03, -1.6667612931866305e-03;
        polygonEdgeTangents.col(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.1;
        const unsigned int curvedEdgeIndex = 0;

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult result =
            geometryUtilities.PolygonDivisionByCircle(polygonVertices,
                                                      polygonEdgeTangents,
                                                      circleCenter,
                                                      circleRadius,
                                                      curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByCircleResult expectedResult;
        expectedResult.Points.setZero(3, 4);
        expectedResult.Points.block(0, 0, 3, 3)<< polygonVertices;
        expectedResult.Points.col(3)<< circleCenter;

        expectedResult.SubTriangles = { {3, 1, 2} };
        expectedResult.InternalTriangles = { {3, 1, 0} };
        expectedResult.SubPolygons = { {1, 2, 0} };

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

  TEST(TestGeometryUtilities, TestPolygonInsideCircleDivisionByAngleQuadrant)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check square sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.row(0)<<  0.0000000000000000e+00,  1.1616496581557671e-02,  1.9133698036719535e-02,  3.7692030705120110e-02;
        polygonVertices.row(1)<<  0.0000000000000000e+00, -9.9322993345804175e-02, -9.8152440618864065e-02,  5.0722730895791608e-02;
        polygonVertices.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.row(0)<<  1.1616496581557671e-02,  7.5172014551618642e-03,  1.8558332668400575e-02, -3.7692030705120110e-02;
        polygonEdgeTangents.row(1)<< -9.9322993345804175e-02,  1.1705527269401106e-03,  1.4887517151465568e-01, -5.0722730895791608e-02;
        polygonEdgeTangents.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 0.1;
        const unsigned int curvedEdgeIndex = 1;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtilities.PolygonInsideCircleDivisionByAngleQuadrant(polygonVertices,
                                                                         polygonEdgeTangents,
                                                                         circleCenter,
                                                                         circleRadius,
                                                                         curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 5);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.Points.col(4)<< 6.9388939039072284e-18, 1.3877787807814457e-17, 0.0;
        expectedResult.SubPolygons = { {2, 3, 4}, { 1, 2, 4, 0 } };
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

  TEST(TestGeometryUtilities, TestPolygonDivisionByAngleQuadrant)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      geometryUtilitiesConfig.Tolerance = 1.0e-12;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check square sub-division
      {
        Eigen::MatrixXd polygonVertices(3, 4);
        polygonVertices.row(0)<<  6.3425787364898834e-01,  8.6113631159405246e-01,  1.0000000000000000e+00,  7.7312156205493587e-01 ;
        polygonVertices.row(1)<< -7.7312156205493587e-01, -1.0000000000000000e+00, -8.6113631159405246e-01, -6.3425787364898834e-01;
        polygonVertices.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        Eigen::MatrixXd polygonEdgeTangents(3, 4);
        polygonEdgeTangents.row(0)<<  2.2687843794506413e-01,  1.3886368840594754e-01, -2.2687843794506413e-01, -1.3886368840594754e-01 ;
        polygonEdgeTangents.row(1)<< -2.2687843794506413e-01,  1.3886368840594754e-01,  2.2687843794506413e-01, -1.3886368840594754e-01;
        polygonEdgeTangents.row(2)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;

        const Eigen::Vector3d circleCenter(0.0, 0.0, 0.0);
        const double circleRadius = 1.0;
        const unsigned int curvedEdgeIndex = 3;

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult result =
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
                                                                          polygonEdgeTangents,
                                                                          circleCenter,
                                                                          circleRadius,
                                                                          curvedEdgeIndex);

        Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult expectedResult;
        expectedResult.Points.setZero(3, 4);
        expectedResult.Points.block(0, 0, 3, 4)<< polygonVertices;
        expectedResult.SubPolygons = { {3, 0, 1, 2} };
        expectedResult.SubPolygonTypes =
        {
          Gedim::GeometryUtilities::PolygonDivisionByAngleQuadrantResult::Types::Internal
        };

        ASSERT_EQ(result.Points, expectedResult.Points);
        ASSERT_EQ(result.SubPolygons, expectedResult.SubPolygons);
        ASSERT_EQ(result.SubPolygonTypes, expectedResult.SubPolygonTypes);
      }

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
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
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
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
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
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
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
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
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
            geometryUtilities.PolygonOutsideCircleDivisionByAngleQuadrant(polygonVertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check circle outside and no intersection
      {
        Eigen::Matrix3d polygonVertices;
        polygonVertices.col(0)<< 0.0, 0.0, 0.0;
        polygonVertices.col(1)<< 1.0, 0.0, 0.0;
        polygonVertices.col(2)<< 0.0, 1.0, 0.0;
        Eigen::Vector3d circleCenter(0.0, 3.0, 0.0);
        double circleRadius = 1.0;

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);
        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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

        Gedim::GeometryUtilities::IntersectionPolygonCircleResult polygonCircleIntersections = geometryUtilities.IntersectionPolygonCircle(polygonVertices,
                                                                                                                                           circleCenter,
                                                                                                                                           circleRadius);

        vector<Gedim::GeometryUtilities::PointCirclePositionResult> vertexPositions = geometryUtilities.PointCirclePositions(polygonVertices,
                                                                                                                             circleCenter,
                                                                                                                             circleRadius);
        Gedim::GeometryUtilities::PolygonCirclePositionTypes position = geometryUtilities.PolygonCirclePosition(polygonVertices,
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
