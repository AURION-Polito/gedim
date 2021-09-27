#ifndef __TEST_UnionMeshSegment_H
#define __TEST_UnionMeshSegment_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "UnionMeshSegment.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

  TEST(TestUnionMeshSegment, TestUnionMeshSegment)
  {
    try
    {
      Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
      Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

      // unify simple two points mesh
      {
        Gedim::IntersectorMesh2DSegment::IntersectionMesh meshOne;
        Gedim::IntersectorMesh2DSegment::IntersectionMesh meshTwo;
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                1.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));

        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                1.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));

        vector<double> curvilinearCoordinatesMeshOne;
        vector<double> curvilinearCoordinatesMeshTwo;
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(meshOne, curvilinearCoordinatesMeshOne);
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(meshTwo, curvilinearCoordinatesMeshTwo);

        Gedim::UnionMeshSegment UnionMeshSegment(geometryUtilities);

        Gedim::UnionMeshSegment::UnionMesh result;
        ASSERT_NO_THROW(UnionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne,
                                                         curvilinearCoordinatesMeshTwo,
                                                         result));

        EXPECT_EQ(result.Points.size(), 2);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Both);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[0], 0);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[1], 0);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Both);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[0], 1);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[1], 1);
      }

      // unify four points mesh with four points mesh
      {
        Gedim::IntersectorMesh2DSegment::IntersectionMesh meshOne;
        Gedim::IntersectorMesh2DSegment::IntersectionMesh meshTwo;
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.25,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.5,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshOne.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                1.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));

        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.5,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                0.75,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));
        meshTwo.Points.insert(pair<double,
                              Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint>(
                                1.0,
                                Gedim::IntersectorMesh2DSegment::IntersectionMesh::IntersectionMeshPoint()));

        vector<double> curvilinearCoordinatesMeshOne;
        vector<double> curvilinearCoordinatesMeshTwo;
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(meshOne, curvilinearCoordinatesMeshOne);
        Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(meshTwo, curvilinearCoordinatesMeshTwo);

        Gedim::UnionMeshSegment UnionMeshSegment(geometryUtilities);

        Gedim::UnionMeshSegment::UnionMesh result;
        ASSERT_NO_THROW(UnionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne,
                                                         curvilinearCoordinatesMeshTwo,
                                                         result));

        EXPECT_EQ(result.Points.size(), 5);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Both);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[0], 0);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[1], 0);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::First);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].MeshIndices[0], 1);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Both);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices[0], 2);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices[1], 1);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Second);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].MeshIndices[1], 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::Both);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[0], 3);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[1], 3);
      }
    }
    catch (const exception& exception)
    {
      cerr<< exception.what()<< endl;
      FAIL();
    }
  }
}

#endif // __TEST_UnionMeshSegment_H
