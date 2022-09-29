#ifndef __TEST_GEOMETRY_POLYHEDRON_H
#define __TEST_GEOMETRY_POLYHEDRON_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronBarycenter)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(geometryUtilities.PolyhedronBarycenter(cube.Vertices),
                                                          Eigen::Vector3d(0.5, 0.5, 0.5)));
      }

      // check tetrahedron barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices),
                                                          Eigen::Vector3d(0.25, 0.25, 0.25)));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronEdgeTangents)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        Eigen::MatrixXd expectedEdgeTangents(3, 12);
        expectedEdgeTangents.col(0)<<   1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(1)<<   0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(2)<<  -1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(3)<<   0.0000000000000000e+00, -1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(4)<<   1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(5)<<   0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(6)<<  -1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(7)<<   0.0000000000000000e+00, -1.0000000000000000e+00,  0.0000000000000000e+00;
        expectedEdgeTangents.col(8)<<   0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(9)<<   0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(10)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;
        expectedEdgeTangents.col(11)<<  0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00;

        ASSERT_EQ(geometryUtilities.PolyhedronEdgeTangents(cube.Vertices,
                                                           cube.Edges),
                  expectedEdgeTangents);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        Eigen::MatrixXd expectedEdgeTangents(3, 6);
        expectedEdgeTangents.col(0)<<  1.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(1)<<  0.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(2)<< -1.0000000000000000e+00,  1.0000000000000000e+00, 0.0000000000000000e+00;
        expectedEdgeTangents.col(3)<<  0.0000000000000000e+00,  0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedEdgeTangents.col(4)<< -1.0000000000000000e+00,  0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedEdgeTangents.col(5)<<  0.0000000000000000e+00, -1.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtilities.PolyhedronEdgeTangents(tetrahedron.Vertices,
                                                           tetrahedron.Edges),
                  expectedEdgeTangents);;

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceVertices)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        vector<Eigen::MatrixXd> expectedFaceVertices(6, Eigen::MatrixXd(3, 4));
        expectedFaceVertices[0]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedFaceVertices[1]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[2]<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[3]<< 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[4]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[5]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                           cube.Faces),
                  expectedFaceVertices);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        vector<Eigen::MatrixXd> expectedFaceVertices(4, Eigen::MatrixXd(3, 3));

        expectedFaceVertices[0]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        expectedFaceVertices[1]<< 0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[2]<< 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceVertices[3]<< 1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00;

        ASSERT_EQ(geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                           tetrahedron.Faces),
                  expectedFaceVertices);

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceNormals)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face normals
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::Vector3d barycenter(0.5, 0.5, 0.5);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(6);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[4]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[5]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(geometryUtilities.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices,
                                                                   barycenter,
                                                                   expectedFaceNormals),
                  vector<bool>({ false, true, false, true, true, false }));

      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::Vector3d barycenter(0.25, 0.25, 0.25);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(4);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +5.7735026918962584e-01, +5.7735026918962584e-01, +5.7735026918962584e-01;

        ASSERT_EQ(geometryUtilities.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtilities.PolyhedronFaceNormalDirections(faceVertices,
                                                                   barycenter,
                                                                   expectedFaceNormals),
                  vector<bool>({ false, true, false, true }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeDirections)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face edge directions
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        vector<vector<bool>> expectedFaceEdgeDirections(6);
        expectedFaceEdgeDirections[0] = vector<bool>({ true, true, true, true });
        expectedFaceEdgeDirections[1] = vector<bool>({ true, true, true, true });
        expectedFaceEdgeDirections[2] = vector<bool>({ false, true, true, false });
        expectedFaceEdgeDirections[3] = vector<bool>({ true, true, false, false });
        expectedFaceEdgeDirections[4] = vector<bool>({ true, true, false, false });
        expectedFaceEdgeDirections[5] = vector<bool>({ false, true, true, false });

        ASSERT_EQ(expectedFaceEdgeDirections,
                  geometryUtilities.PolyhedronFaceEdgeDirections(cube.Vertices,
                                                                 cube.Edges,
                                                                 cube.Faces));
      }

      // check tetrahedron face edge directions
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        vector<vector<bool>> expectedFaceEdgeDirections(4);
        expectedFaceEdgeDirections[0] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[1] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[2] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[3] = vector<bool>({ true, true, false });

        ASSERT_EQ(expectedFaceEdgeDirections,
                  geometryUtilities.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
                                                                 tetrahedron.Edges,
                                                                 tetrahedron.Faces));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceEdgeTangents)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::MatrixXd edgeTangents = geometryUtilities.PolyhedronEdgeTangents(cube.Vertices,
                                                                                      cube.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(cube.Vertices,
                                                                                                       cube.Edges,
                                                                                                       cube.Faces);

        vector<Eigen::MatrixXd> expectedFaceEdgeTangents(6);
        expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1,  0,
                                       0,  1,  0, -1,
                                       0,  0,  0,  0).finished();
        expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1,  0,
                                       0,  1,  0, -1,
                                       0,  0,  0,  0).finished();
        expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 4)<<
                                       -0,  0,  0, -0,
                                       1,  0, -1, -0,
                                       -0,  1,  0, -1).finished();
        expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 4)<<
                                       0,  0, -0, -0,
                                       1,  0, -1, -0,
                                       0,  1, -0, -1).finished();
        expectedFaceEdgeTangents[4] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1, -0,
                                       0,  0, -0, -0,
                                       0,  1, -0, -1).finished();
        expectedFaceEdgeTangents[5] = (Eigen::MatrixXd(3, 4)<<
                                       1,  0, -1, -0,
                                       -0,  0,  0, -0,
                                       -0,  1,  0, -1).finished();

        ASSERT_EQ(geometryUtilities.PolyhedronFaceEdgeTangents(cube.Vertices,
                                                               cube.Edges,
                                                               cube.Faces,
                                                               faceEdgeDirections,
                                                               edgeTangents),
                  expectedFaceEdgeTangents);

      }

      // check tetrahedron face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::MatrixXd edgeTangents = geometryUtilities.PolyhedronEdgeTangents(tetrahedron.Vertices,
                                                                                      tetrahedron.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtilities.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
                                                                                                       tetrahedron.Edges,
                                                                                                       tetrahedron.Faces);


        vector<Eigen::MatrixXd> expectedFaceEdgeTangents(4);
        expectedFaceEdgeTangents[0] = (Eigen::MatrixXd(3, 3)<<
                                       1, -1, -0,
                                       0,  1, -1,
                                       0,  0, -0).finished();
        expectedFaceEdgeTangents[1] = (Eigen::MatrixXd(3, 3)<<
                                       1, -1, -0,
                                       0,  0, -0,
                                       0,  1, -1).finished();
        expectedFaceEdgeTangents[2] = (Eigen::MatrixXd(3, 3)<<
                                       0,  0, -0,
                                       1, -1, -0,
                                       0,  1, -1).finished();
        expectedFaceEdgeTangents[3] = (Eigen::MatrixXd(3, 3)<<
                                       -1,  0,  1,
                                       1, -1, -0,
                                       0,  1, -1).finished();

        ASSERT_EQ(geometryUtilities.PolyhedronFaceEdgeTangents(tetrahedron.Vertices,
                                                               tetrahedron.Edges,
                                                               tetrahedron.Faces,
                                                               faceEdgeDirections,
                                                               edgeTangents),
                  expectedFaceEdgeTangents);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceTranslations)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face translations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(6);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[4]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[5]<< +0.0000000000000000e+00, +1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtilities.PolyhedronFaceTranslations(faceVertices));

      }

      // check tetrahedron face translations
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(4);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtilities.PolyhedronFaceTranslations(faceVertices));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceRotationMatrices)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);
        const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(cube.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);

        vector<Eigen::Matrix3d> expectedFaceRotationMatrices(6);
        expectedFaceRotationMatrices[0] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 9.9999999999999978e-01, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[1] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 9.9999999999999978e-01, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[2] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00,
                                           9.9999999999999989e-01, -1.1102230246251565e-16,  0.0000000000000000e+00,
                                           1.1102230246251565e-16,  9.9999999999999989e-01,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[3] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00,  0.0000000000000000e+00,  1.0000000000000000e+00,
                                           9.9999999999999989e-01, -1.1102230246251565e-16,  0.0000000000000000e+00,
                                           1.1102230246251565e-16,  9.9999999999999989e-01,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[4] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01,  0.0000000000000000e+00,  4.9650683064945600e-17,
                                           4.9650683064945576e-17,  2.4825341532472850e-17, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  9.9999999999999978e-01,  2.4825341532472825e-17).finished();
        expectedFaceRotationMatrices[5] = (Eigen::Matrix3d()<<
                                           9.9999999999999978e-01,  0.0000000000000000e+00,  4.9650683064945600e-17,
                                           4.9650683064945576e-17,  2.4825341532472850e-17, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  9.9999999999999978e-01,  2.4825341532472825e-17).finished();

        ASSERT_EQ(expectedFaceRotationMatrices,
                  geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices,
                                                                   faceNormals,
                                                                   faceTranslations));

      }

      // check tetrahedron face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);

        const Eigen::Vector3d barycenter = geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtilities.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtilities.PolyhedronFaceTranslations(faceVertices);

        vector<Eigen::Matrix3d> expectedFaceRotationMatrices(4);
        expectedFaceRotationMatrices[0] = (Eigen::Matrix3d()<<
                                           1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[1] = (Eigen::Matrix3d()<<
                                           1.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
                                           0.0000000000000000e+00,  0.0000000000000000e+00, -1.0000000000000000e+00,
                                           0.0000000000000000e+00,  1.0000000000000000e+00,  0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[2] = (Eigen::Matrix3d()<<
                                           0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00,
                                           1.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                                           0.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00).finished();
        expectedFaceRotationMatrices[3] = (Eigen::Matrix3d()<<
                                           -7.0710678118654724e-01, -4.0824829046386296e-01,  5.7735026918962595e-01,
                                           7.0710678118654724e-01, -4.0824829046386313e-01,  5.7735026918962573e-01,
                                           0.0000000000000000e+00,  8.1649658092772581e-01,  5.7735026918962651e-01).finished();

        ASSERT_EQ(expectedFaceRotationMatrices,
                  geometryUtilities.PolyhedronFaceRotationMatrices(faceVertices,
                                                                   faceNormals,
                                                                   faceTranslations));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronCoordinateSystem)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube coordinate system
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_EQ(geometryUtilities.PolyhedronCoordinateSystem(cube.Vertices,
                                                               cube.Edges),
                  vector<unsigned int>({ 0, 1, 3, 4 }));
      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_EQ(geometryUtilities.PolyhedronCoordinateSystem(tetrahedron.Vertices,
                                                               tetrahedron.Edges),
                  vector<unsigned int>({ 0, 1, 2, 3 }));
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceBarycenters)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face baycenters
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);
        ASSERT_EQ(geometryUtilities.PolyhedronFaceBarycenter(faceVertices),
                  vector<Eigen::Vector3d>({
                                            Eigen::Vector3d(0.5, 0.5, 0.0),
                                            Eigen::Vector3d(0.5, 0.5, 1.0),
                                            Eigen::Vector3d(0.0, 0.5, 0.5),
                                            Eigen::Vector3d(1.0, 0.5, 0.5),
                                            Eigen::Vector3d(0.5, 0.0, 0.5),
                                            Eigen::Vector3d(0.5, 1.0, 0.5)
                                          }));
      }

      // check tetrahedron face baycenters
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);

        ASSERT_EQ(geometryUtilities.PolyhedronFaceBarycenter(faceVertices),
                  vector<Eigen::Vector3d>({
                                            Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 0.0),
                                            Eigen::Vector3d(1.0 / 3.0, 0.0, 1.0 / 3.0),
                                            Eigen::Vector3d(0.0, 1.0 / 3.0, 1.0 / 3.0),
                                            Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
                                          }));

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronFaceTriangulations)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // check cube face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);

        const vector<vector<unsigned int>> faceTriangulationsByFirstVertex = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces,
                                                                                                                                         faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices,
                                                                                                                                             cube.Faces,
                                                                                                                                             faceVertices,
                                                                                                                                             faceBarycenters);

        ASSERT_EQ(faceTriangulationsByFirstVertex,
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 0,1,2,0,2,3 })
                                               }));
        ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(faceVertices,
                                                                                   faceTriangulationsByFirstVertex),
                  vector<vector<Eigen::Matrix3d>>({
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(2), cube.Vertices.col(3)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(4), cube.Vertices.col(5), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(4), cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(3), cube.Vertices.col(7)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(1), cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(1), cube.Vertices.col(6), cube.Vertices.col(5)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(1), cube.Vertices.col(5)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(0), cube.Vertices.col(5), cube.Vertices.col(4)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(3), cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< cube.Vertices.col(3), cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                                                    })
                                                  })
                  );


        ASSERT_EQ(faceTriangulationsByInternalPoint,
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 }),
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 }),
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 }),
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 }),
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 }),
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,3,4,3,0 })
                                               }));

        ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationPointsByInternalPoint(faceVertices,
                                                                                     faceBarycenters,
                                                                                     faceTriangulationsByInternalPoint),
                  vector<vector<Eigen::Matrix3d>>({
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[0], cube.Vertices.col(0), cube.Vertices.col(1)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[0], cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[0], cube.Vertices.col(2), cube.Vertices.col(3)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[0], cube.Vertices.col(3), cube.Vertices.col(0)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[1], cube.Vertices.col(4), cube.Vertices.col(5)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[1], cube.Vertices.col(5), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[1], cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[1], cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[2], cube.Vertices.col(0), cube.Vertices.col(3)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[2], cube.Vertices.col(3), cube.Vertices.col(7)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[2], cube.Vertices.col(7), cube.Vertices.col(4)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[2], cube.Vertices.col(4), cube.Vertices.col(0)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[3], cube.Vertices.col(1), cube.Vertices.col(2)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[3], cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[3], cube.Vertices.col(6), cube.Vertices.col(5)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[3], cube.Vertices.col(5), cube.Vertices.col(1)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[4], cube.Vertices.col(0), cube.Vertices.col(1)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[4], cube.Vertices.col(1), cube.Vertices.col(5)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[4], cube.Vertices.col(5), cube.Vertices.col(4)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[4], cube.Vertices.col(4), cube.Vertices.col(0)).finished(),
                                                    }),
                                                    vector<Eigen::Matrix3d>({
                                                      (Eigen::Matrix3d()<< faceBarycenters[5], cube.Vertices.col(3), cube.Vertices.col(2)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[5], cube.Vertices.col(2), cube.Vertices.col(6)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[5], cube.Vertices.col(6), cube.Vertices.col(7)).finished(),
                                                      (Eigen::Matrix3d()<< faceBarycenters[5], cube.Vertices.col(7), cube.Vertices.col(3)).finished(),
                                                    })
                                                  })
                  );
      }

      // check tetrahedron face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);

        ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces,
                                                                              faceVertices),
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 0,1,2 }),
                                                 vector<unsigned int>({ 0,1,2 }),
                                                 vector<unsigned int>({ 0,1,2 }),
                                                 vector<unsigned int>({ 0,1,2 })
                                               }));
        ASSERT_EQ(geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices,
                                                                                tetrahedron.Faces,
                                                                                faceVertices,
                                                                                faceBarycenters),
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 3,0,1,3,1,2,3,2,0 }),
                                                 vector<unsigned int>({ 3,0,1,3,1,2,3,2,0 }),
                                                 vector<unsigned int>({ 3,0,1,3,1,2,3,2,0 }),
                                                 vector<unsigned int>({ 3,0,1,3,1,2,3,2,0 })
                                               }));

      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronTetrahedrons)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      std::string exportFolder = "./Export/TestPolyhedronTetrahedrons";
      Gedim::Output::CreateFolder(exportFolder);

      // check cube face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtilities.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                           Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                           Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(cube.Vertices);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(cube.Vertices,
                                                                                              cube.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
        const vector<vector<unsigned int>> faceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces,
                                                                                                                            faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices,
                                                                                                                                             cube.Faces,
                                                                                                                                             faceVertices,
                                                                                                                                             faceBarycenters);
        const vector<unsigned int> tetrahedronList = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices,
                                                                                                                  cube.Faces,
                                                                                                                  faceTriangulations,
                                                                                                                  polyhedronBarycenter);
        const vector<unsigned int> tetrahedronByInternalPointsList = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices,
                                                                                                                                  cube.Faces,
                                                                                                                                  faceTriangulationsByInternalPoint,
                                                                                                                                  faceBarycenters,
                                                                                                                                  polyhedronBarycenter);
        ASSERT_EQ(tetrahedronList,
                  vector<unsigned int>({ 0,1,2,8,0,2,3,8,
                                         4,5,6,8,4,6,7,8,
                                         0,3,7,8,0,7,4,8,
                                         1,2,6,8,1,6,5,8,
                                         0,1,5,8,0,5,4,8,
                                         3,2,6,8,3,6,7,8 }));
        ASSERT_EQ(tetrahedronByInternalPointsList,
                  vector<unsigned int>({
                                         8 ,0,1,14,8 ,1,2,14,8 ,2,3,14,8 ,3,0,14,
                                         9 ,4,5,14,9 ,5,6,14,9 ,6,7,14,9 ,7,4,14,
                                         10,0,3,14,10,3,7,14,10,7,4,14,10,4,0,14,
                                         11,1,2,14,11,2,6,14,11,6,5,14,11,5,1,14,
                                         12,0,1,14,12,1,5,14,12,5,4,14,12,4,0,14,
                                         13,3,2,14,13,2,6,14,13,6,7,14,13,7,3,14
                                       }));
        // Export tetrahedrons
        {
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtilities.ExtractTetrahedronPoints(cube.Vertices,
                                                                                            polyhedronBarycenter,
                                                                                            tetrahedronList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                                                            tetrahedrons[t].col(1),
                                                                                                            tetrahedrons[t].col(2),
                                                                                                            tetrahedrons[t].col(3));
            vector<double> id(subTetra.Faces.size(), t);

            vtkExperter.AddPolyhedron(subTetra.Vertices,
                                      subTetra.Edges,
                                      subTetra.Faces,
                                      {
                                        {
                                          "Id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(id.size()),
                                          id.data()
                                        }
                                      });
          }

          vtkExperter.Export(exportFolder + "/Cube_Tetra_1.vtu",
                             Gedim::VTKUtilities::Ascii);
        }

        {
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtilities.ExtractTetrahedronPoints(cube.Vertices,
                                                                                            polyhedronBarycenter,
                                                                                            faceBarycenters,
                                                                                            tetrahedronByInternalPointsList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                                                            tetrahedrons[t].col(1),
                                                                                                            tetrahedrons[t].col(2),
                                                                                                            tetrahedrons[t].col(3));
            vector<double> id(subTetra.Faces.size(), t);

            vtkExperter.AddPolyhedron(subTetra.Vertices,
                                      subTetra.Edges,
                                      subTetra.Faces,
                                      {
                                        {
                                          "Id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(id.size()),
                                          id.data()
                                        }
                                      });
          }

          vtkExperter.Export(exportFolder + "/Cube_Tetra_2.vtu",
                             Gedim::VTKUtilities::Ascii);
        }
      }

      // check tetrahedron face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                               Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                               Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(tetrahedron.Vertices);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtilities.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                              tetrahedron.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtilities.PolyhedronFaceBarycenter(faceVertices);
        const vector<vector<unsigned int>> faceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces,
                                                                                                                            faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtilities.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices,
                                                                                                                                             tetrahedron.Faces,
                                                                                                                                             faceVertices,
                                                                                                                                             faceBarycenters);

        const vector<unsigned int> tetrahedronList = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
                                                                                                                  tetrahedron.Faces,
                                                                                                                  faceTriangulations,
                                                                                                                  polyhedronBarycenter);
        const vector<unsigned int> tetrahedronByInternalPointsList = geometryUtilities.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
                                                                                                                                  tetrahedron.Faces,
                                                                                                                                  faceTriangulationsByInternalPoint,
                                                                                                                                  faceBarycenters,
                                                                                                                                  polyhedronBarycenter);
        ASSERT_EQ(tetrahedronList,
                  vector<unsigned int>({ 0,1,2,4,
                                         0,1,3,4,
                                         0,2,3,4,
                                         1,2,3,4 }));
        ASSERT_EQ(tetrahedronByInternalPointsList,
                  vector<unsigned int>({ 4,0,1,8,4,1,2,8,4,2,0,8,
                                         5,0,1,8,5,1,3,8,5,3,0,8,
                                         6,0,2,8,6,2,3,8,6,3,0,8,
                                         7,1,2,8,7,2,3,8,7,3,1,8 }));

        // Export tetrahedrons
        {
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtilities.ExtractTetrahedronPoints(tetrahedron.Vertices,
                                                                                            polyhedronBarycenter,
                                                                                            tetrahedronList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                                                            tetrahedrons[t].col(1),
                                                                                                            tetrahedrons[t].col(2),
                                                                                                            tetrahedrons[t].col(3));
            vector<double> id(subTetra.Faces.size(), t);

            vtkExperter.AddPolyhedron(subTetra.Vertices,
                                      subTetra.Edges,
                                      subTetra.Faces,
                                      {
                                        {
                                          "Id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(id.size()),
                                          id.data()
                                        }
                                      });
          }

          vtkExperter.Export(exportFolder + "/Tetrahedrons_Tetra_1.vtu",
                             Gedim::VTKUtilities::Ascii);
        }

        {
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtilities.ExtractTetrahedronPoints(tetrahedron.Vertices,
                                                                                            polyhedronBarycenter,
                                                                                            faceBarycenters,
                                                                                            tetrahedronByInternalPointsList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtilities.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
                                                                                                            tetrahedrons[t].col(1),
                                                                                                            tetrahedrons[t].col(2),
                                                                                                            tetrahedrons[t].col(3));
            vector<double> id(subTetra.Faces.size(), t);

            vtkExperter.AddPolyhedron(subTetra.Vertices,
                                      subTetra.Edges,
                                      subTetra.Faces,
                                      {
                                        {
                                          "Id",
                                          Gedim::VTPProperty::Formats::Cells,
                                          static_cast<unsigned int>(id.size()),
                                          id.data()
                                        }
                                      });
          }

          vtkExperter.Export(exportFolder + "/Tetrahedrons_Tetra_2.vtu",
                             Gedim::VTKUtilities::Ascii);
        }
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronVolume)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-15;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    // check cube volume
    {
      const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                     1.0);

      const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    polyhedronFace3DVertices);

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

      const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(polyhedronFace2DVertices,
                                                                                                                                                             polyhedronFaceTriangulations);

      const double polyhedronVolume = geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                                         polyhedronFaceNormals,
                                                                         polyhedronFaceNormalDirections,
                                                                         polyhedronFaceTranslations,
                                                                         polyhedronFaceRotationMatrices);

      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0, polyhedronVolume));
    }

    // check tetrahedron volume
    {
      const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                            Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                            Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                            Eigen::Vector3d(0.0,1.0,0.0));

      const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    polyhedronFace3DVertices);

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

      const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(polyhedronFace2DVertices,
                                                                                                                                                             polyhedronFaceTriangulations);

      const double polyhedronVolume = geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                                         polyhedronFaceNormals,
                                                                         polyhedronFaceNormalDirections,
                                                                         polyhedronFaceTranslations,
                                                                         polyhedronFaceRotationMatrices);

      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(1.0/6.0, polyhedronVolume));
    }
  }

  TEST(TestGeometryUtilities, TestPolyhedron_TestPolyhedronCentroid)
  {
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance = 1.0e-14;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    // check cube centroid
    {
      const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                     1.0);

      const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    polyhedronFace3DVertices);

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

      const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(polyhedronFace2DVertices,
                                                                                                                                                             polyhedronFaceTriangulations);

      const double polyhedronVolume = geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                                         polyhedronFaceNormals,
                                                                         polyhedronFaceNormalDirections,
                                                                         polyhedronFaceTranslations,
                                                                         polyhedronFaceRotationMatrices);

      const Eigen::Vector3d polyhedronCentroid = geometryUtilities.PolyhedronCentroid(polyhedronFace2DTriangulationPoints,
                                                                                      polyhedronFaceNormals,
                                                                                      polyhedronFaceNormalDirections,
                                                                                      polyhedronFaceTranslations,
                                                                                      polyhedronFaceRotationMatrices,
                                                                                      polyhedronVolume);

      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.x(), 0.5));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.y(), 0.5));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.z(), 0.5));
    }

    // check tetrahedron volume
    {
      const Gedim::GeometryUtilities::Polyhedron polyhedron = geometryUtilities.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                            Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                            Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                            Eigen::Vector3d(0.0,1.0,0.0));

      const Eigen::Vector3d polyhedronBarycenter = geometryUtilities.PolyhedronBarycenter(polyhedron.Vertices);
      const vector<Eigen::MatrixXd> polyhedronFace3DVertices = geometryUtilities.PolyhedronFaceVertices(polyhedron.Vertices,
                                                                                                        polyhedron.Faces);
      const vector<vector<unsigned int>> polyhedronFaceTriangulations = geometryUtilities.PolyhedronFaceTriangulationsByFirstVertex(polyhedron.Faces,
                                                                                                                                    polyhedronFace3DVertices);

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

      const std::vector<std::vector<Eigen::Matrix3d>> polyhedronFace2DTriangulationPoints = geometryUtilities.PolyhedronFaceTriangulationPointsByFirstVertex(polyhedronFace2DVertices,
                                                                                                                                                             polyhedronFaceTriangulations);

      const double polyhedronVolume = geometryUtilities.PolyhedronVolume(polyhedronFace2DTriangulationPoints,
                                                                         polyhedronFaceNormals,
                                                                         polyhedronFaceNormalDirections,
                                                                         polyhedronFaceTranslations,
                                                                         polyhedronFaceRotationMatrices);

      const Eigen::Vector3d polyhedronCentroid = geometryUtilities.PolyhedronCentroid(polyhedronFace2DTriangulationPoints,
                                                                                      polyhedronFaceNormals,
                                                                                      polyhedronFaceNormalDirections,
                                                                                      polyhedronFaceTranslations,
                                                                                      polyhedronFaceRotationMatrices,
                                                                                      polyhedronVolume);

      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.x(), 1.0 / 4.0));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.y(), 1.0 / 4.0));
      ASSERT_TRUE(geometryUtilities.Are1DValuesEqual(polyhedronCentroid.z(), 1.0 / 4.0));
    }
  }
}

#endif // __TEST_GEOMETRY_POLYHEDRON_H
