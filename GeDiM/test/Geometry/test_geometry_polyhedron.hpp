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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(cube.Vertices),
                                                        Eigen::Vector3d(0.5, 0.5, 0.5)));
      }

      // check tetrahedron barycenter
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_TRUE(geometryUtility.PointsAreCoincident(geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices),
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
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

        ASSERT_EQ(geometryUtility.PolyhedronEdgeTangents(cube.Vertices,
                                                         cube.Edges),
                  expectedEdgeTangents);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
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

        ASSERT_EQ(geometryUtility.PolyhedronEdgeTangents(tetrahedron.Vertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
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

        ASSERT_EQ(geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                         cube.Faces),
                  expectedFaceVertices);
      }

      // check tetrahedron face vertices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
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

        ASSERT_EQ(geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face normals
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::Vector3d barycenter(0.5, 0.5, 0.5);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(6);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[4]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[5]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtility.PolyhedronFaceNormalDirections(faceVertices,
                                                                 barycenter,
                                                                 expectedFaceNormals),
                  vector<bool>({ false, true, false, true, true, false }));

      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::Vector3d barycenter(0.25, 0.25, 0.25);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceNormals(4);
        expectedFaceNormals[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, 1.0000000000000000e+00;
        expectedFaceNormals[1]<< +0.0000000000000000e+00, -1.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[2]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceNormals[3]<< +5.7735026918962584e-01, +5.7735026918962584e-01, +5.7735026918962584e-01;

        ASSERT_EQ(geometryUtility.PolyhedronFaceNormals(faceVertices),
                  expectedFaceNormals);
        ASSERT_EQ(geometryUtility.PolyhedronFaceNormalDirections(faceVertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face edge directions
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
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
                  geometryUtility.PolyhedronFaceEdgeDirections(cube.Vertices,
                                                               cube.Edges,
                                                               cube.Faces));
      }

      // check tetrahedron face edge directions
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        vector<vector<bool>> expectedFaceEdgeDirections(4);
        expectedFaceEdgeDirections[0] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[1] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[2] = vector<bool>({ true, true, false });
        expectedFaceEdgeDirections[3] = vector<bool>({ true, true, false });

        ASSERT_EQ(expectedFaceEdgeDirections,
                  geometryUtility.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::MatrixXd edgeTangents = geometryUtility.PolyhedronEdgeTangents(cube.Vertices,
                                                                                    cube.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(cube.Vertices,
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

        ASSERT_EQ(geometryUtility.PolyhedronFaceEdgeTangents(cube.Vertices,
                                                             cube.Edges,
                                                             cube.Faces,
                                                             faceEdgeDirections,
                                                             edgeTangents),
                  expectedFaceEdgeTangents);

      }

      // check tetrahedron face edge tangents
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::MatrixXd edgeTangents = geometryUtility.PolyhedronEdgeTangents(tetrahedron.Vertices,
                                                                                    tetrahedron.Edges);
        const vector<vector<bool>> faceEdgeDirections = geometryUtility.PolyhedronFaceEdgeDirections(tetrahedron.Vertices,
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

        ASSERT_EQ(geometryUtility.PolyhedronFaceEdgeTangents(tetrahedron.Vertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face translations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(6);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +1.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[4]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[5]<< +0.0000000000000000e+00, +1.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtility.PolyhedronFaceTranslations(faceVertices));

      }

      // check tetrahedron face translations
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        vector<Eigen::Vector3d> expectedFaceTranslations(4);
        expectedFaceTranslations[0]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[1]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[2]<< +0.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;
        expectedFaceTranslations[3]<< +1.0000000000000000e+00, +0.0000000000000000e+00, +0.0000000000000000e+00;

        ASSERT_EQ(expectedFaceTranslations,
                  geometryUtility.PolyhedronFaceTranslations(faceVertices));
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);
        const Eigen::Vector3d barycenter = geometryUtility.PolyhedronBarycenter(cube.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtility.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtility.PolyhedronFaceTranslations(faceVertices);

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
                  geometryUtility.PolyhedronFaceRotationMatrices(faceVertices,
                                                                 faceNormals,
                                                                 faceTranslations));

      }

      // check tetrahedron face rotation matrices
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        const Eigen::Vector3d barycenter = geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices);
        const vector<Eigen::Vector3d> faceNormals = geometryUtility.PolyhedronFaceNormals(faceVertices);
        const vector<Eigen::Vector3d> faceTranslations = geometryUtility.PolyhedronFaceTranslations(faceVertices);

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
                  geometryUtility.PolyhedronFaceRotationMatrices(faceVertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube coordinate system
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));

        ASSERT_EQ(geometryUtility.PolyhedronCoordinateSystem(cube.Vertices,
                                                             cube.Edges),
                  vector<unsigned int>({ 0, 1, 3, 4 }));
      }

      // check tetrahedron face normals
      {
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));
        ASSERT_EQ(geometryUtility.PolyhedronCoordinateSystem(tetrahedron.Vertices,
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face baycenters
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);
        ASSERT_EQ(geometryUtility.PolyhedronFaceBarycenter(faceVertices),
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
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);

        ASSERT_EQ(geometryUtility.PolyhedronFaceBarycenter(faceVertices),
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      // check cube face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtility.PolyhedronFaceBarycenter(faceVertices);

        const vector<vector<unsigned int>> faceTriangulationsByFirstVertex = geometryUtility.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces,
                                                                                                                                       faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtility.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices,
                                                                                                                                           cube.Faces,
                                                                                                                                           faceVertices,
                                                                                                                                           faceBarycenters);

        ASSERT_EQ(faceTriangulationsByFirstVertex,
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 0,1,2,0,2,3 }),
                                                 vector<unsigned int>({ 4,5,6,4,6,7 }),
                                                 vector<unsigned int>({ 0,3,7,0,7,4 }),
                                                 vector<unsigned int>({ 1,2,6,1,6,5 }),
                                                 vector<unsigned int>({ 0,1,5,0,5,4 }),
                                                 vector<unsigned int>({ 3,2,6,3,6,7 })
                                               }));
        ASSERT_EQ(geometryUtility.PolyhedronFaceTriangulationPointsByFirstVertex(cube.Vertices,
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
                                                 vector<unsigned int>({ 8 ,0,1,8 ,1,2,8 ,2,3,8 ,3,0 }),
                                                 vector<unsigned int>({ 9 ,4,5,9 ,5,6,9 ,6,7,9 ,7,4 }),
                                                 vector<unsigned int>({ 10,0,3,10,3,7,10,7,4,10,4,0 }),
                                                 vector<unsigned int>({ 11,1,2,11,2,6,11,6,5,11,5,1 }),
                                                 vector<unsigned int>({ 12,0,1,12,1,5,12,5,4,12,4,0 }),
                                                 vector<unsigned int>({ 13,3,2,13,2,6,13,6,7,13,7,3 })
                                               }));

        ASSERT_EQ(geometryUtility.PolyhedronFaceTriangulationPointsByInternalPoint(cube.Vertices,
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
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtility.PolyhedronFaceBarycenter(faceVertices);

        ASSERT_EQ(geometryUtility.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces,
                                                                            faceVertices),
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 0,1,2 }),
                                                 vector<unsigned int>({ 0,1,3 }),
                                                 vector<unsigned int>({ 0,2,3 }),
                                                 vector<unsigned int>({ 1,2,3 })
                                               }));
        ASSERT_EQ(geometryUtility.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices,
                                                                              tetrahedron.Faces,
                                                                              faceVertices,
                                                                              faceBarycenters),
                  vector<vector<unsigned int>>({
                                                 vector<unsigned int>({ 4,0,1,4,1,2,4,2,0 }),
                                                 vector<unsigned int>({ 5,0,1,5,1,3,5,3,0 }),
                                                 vector<unsigned int>({ 6,0,2,6,2,3,6,3,0 }),
                                                 vector<unsigned int>({ 7,1,2,7,2,3,7,3,1 })
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
      Gedim::GeometryUtilitiesConfig geometryUtilityConfig;
      Gedim::GeometryUtilities geometryUtility(geometryUtilityConfig);

      std::string exportFolder = "./Export/TestPolyhedronTetrahedrons";
      Gedim::Output::CreateFolder(exportFolder);

      // check cube face triangulations
      {
        const Gedim::GeometryUtilities::Polyhedron cube = geometryUtility.CreateParallelepipedWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                         Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                         Eigen::Vector3d(0.0,1.0,0.0));
        const Eigen::Vector3d polyhedronBarycenter = geometryUtility.PolyhedronBarycenter(cube.Vertices);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(cube.Vertices,
                                                                                            cube.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtility.PolyhedronFaceBarycenter(faceVertices);
        const vector<vector<unsigned int>> faceTriangulations = geometryUtility.PolyhedronFaceTriangulationsByFirstVertex(cube.Faces,
                                                                                                                          faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtility.PolyhedronFaceTriangulationsByInternalPoint(cube.Vertices,
                                                                                                                                           cube.Faces,
                                                                                                                                           faceVertices,
                                                                                                                                           faceBarycenters);
        const vector<unsigned int> tetrahedronList = geometryUtility.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices,
                                                                                                                faceTriangulations,
                                                                                                                polyhedronBarycenter);
        const vector<unsigned int> tetrahedronByInternalPointsList = geometryUtility.PolyhedronTetrahedronsByFaceTriangulations(cube.Vertices,
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
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtility.ExtractTetrahedronPoints(cube.Vertices,
                                                                                          polyhedronBarycenter,
                                                                                          tetrahedronList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtility.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
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
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtility.ExtractTetrahedronPoints(cube.Vertices,
                                                                                          polyhedronBarycenter,
                                                                                          faceBarycenters,
                                                                                          tetrahedronByInternalPointsList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtility.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
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
        const Gedim::GeometryUtilities::Polyhedron tetrahedron = geometryUtility.CreateTetrahedronWithOrigin(Eigen::Vector3d(0.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(1.0,0.0,0.0),
                                                                                                             Eigen::Vector3d(0.0,0.0,1.0),
                                                                                                             Eigen::Vector3d(0.0,1.0,0.0));

        const Eigen::Vector3d polyhedronBarycenter = geometryUtility.PolyhedronBarycenter(tetrahedron.Vertices);
        const vector<Eigen::MatrixXd> faceVertices = geometryUtility.PolyhedronFaceVertices(tetrahedron.Vertices,
                                                                                            tetrahedron.Faces);
        const vector<Eigen::Vector3d> faceBarycenters = geometryUtility.PolyhedronFaceBarycenter(faceVertices);
        const vector<vector<unsigned int>> faceTriangulations = geometryUtility.PolyhedronFaceTriangulationsByFirstVertex(tetrahedron.Faces,
                                                                                                                          faceVertices);
        const vector<vector<unsigned int>> faceTriangulationsByInternalPoint = geometryUtility.PolyhedronFaceTriangulationsByInternalPoint(tetrahedron.Vertices,
                                                                                                                                           tetrahedron.Faces,
                                                                                                                                           faceVertices,
                                                                                                                                           faceBarycenters);

        const vector<unsigned int> tetrahedronList = geometryUtility.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
                                                                                                                faceTriangulations,
                                                                                                                polyhedronBarycenter);
        const vector<unsigned int> tetrahedronByInternalPointsList = geometryUtility.PolyhedronTetrahedronsByFaceTriangulations(tetrahedron.Vertices,
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
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtility.ExtractTetrahedronPoints(tetrahedron.Vertices,
                                                                                          polyhedronBarycenter,
                                                                                          tetrahedronList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtility.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
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
          vector<Eigen::MatrixXd> tetrahedrons = geometryUtility.ExtractTetrahedronPoints(tetrahedron.Vertices,
                                                                                          polyhedronBarycenter,
                                                                                          faceBarycenters,
                                                                                          tetrahedronByInternalPointsList);

          Gedim::VTKUtilities vtkExperter;
          for (unsigned int t = 0; t < tetrahedrons.size(); t++)
          {
            Gedim::GeometryUtilities::Polyhedron subTetra = geometryUtility.CreateTetrahedronWithVertices(tetrahedrons[t].col(0),
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
}

#endif // __TEST_GEOMETRY_POLYHEDRON_H
