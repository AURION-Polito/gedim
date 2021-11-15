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
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[0], 0);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[1], 0);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[0], 1);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[1], 1);
        EXPECT_EQ(result.Segments.size(), 1);
        EXPECT_EQ(result.Segments[0].Points, vector<double>({ 0.0000000000000000e+00, 1.0000000000000000e+00 }));
        EXPECT_EQ(result.Segments[0].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[0].MeshIndices, vector<unsigned int>({ 0, 0 }));
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
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[0], 0);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices[1], 0);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::First);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].MeshIndices[0], 1);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices[0], 2);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices[1], 1);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Second);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].MeshIndices[1], 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices.size(), 2);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[0], 3);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices[1], 3);
        EXPECT_EQ(result.Segments.size(), 4);
        EXPECT_EQ(result.Segments[0].Points, vector<double>({ 0.0000000000000000e+00, 2.5000000000000000e-01 }));
        EXPECT_EQ(result.Segments[0].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[0].MeshIndices, vector<unsigned int>({ 0, 0 }));
        EXPECT_EQ(result.Segments[1].Points, vector<double>({ 2.5000000000000000e-01, 5.0000000000000000e-01 }));
        EXPECT_EQ(result.Segments[1].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[1].MeshIndices, vector<unsigned int>({ 1, 0 }));
        EXPECT_EQ(result.Segments[2].Points, vector<double>({ 5.0000000000000000e-01, 7.5000000000000000e-01 }));
        EXPECT_EQ(result.Segments[2].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[2].MeshIndices, vector<unsigned int>({ 2, 1 }));
        EXPECT_EQ(result.Segments[3].Points, vector<double>({ 7.5000000000000000e-01, 1.0000000000000000e+00 }));
        EXPECT_EQ(result.Segments[3].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[3].MeshIndices, vector<unsigned int>({ 2, 2 }));
      }

      // unify three points mesh with three points mesh
      {
        vector<double> curvilinearCoordinatesMeshOne = { 0.25, 0.75, 1.0 };
        vector<double> curvilinearCoordinatesMeshTwo = { 0.0, 0.5, 0.75 };

        Gedim::UnionMeshSegment UnionMeshSegment(geometryUtilities);

        Gedim::UnionMeshSegment::UnionMesh result;
        ASSERT_NO_THROW(UnionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne,
                                                         curvilinearCoordinatesMeshTwo,
                                                         result));

        EXPECT_EQ(result.Points.size(), 5);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Second);
        EXPECT_EQ(result.Points[0.0000000000000000e+00].MeshIndices, vector<unsigned int>({ 0, 0 }));
        EXPECT_EQ(result.Points[2.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::First);
        EXPECT_EQ(result.Points[2.5000000000000000e-01].MeshIndices, vector<unsigned int>({ 0, 0 }));
        EXPECT_EQ(result.Points[5.0000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Second);
        EXPECT_EQ(result.Points[5.0000000000000000e-01].MeshIndices, vector<unsigned int>({ 0, 1 }));
        EXPECT_EQ(result.Points[7.5000000000000000e-01].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::Both);
        EXPECT_EQ(result.Points[7.5000000000000000e-01].MeshIndices, vector<unsigned int>({ 1, 2 }));
        EXPECT_EQ(result.Points[1.0000000000000000e+00].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshPoint::UnionMeshPoint::Types::First);
        EXPECT_EQ(result.Points[1.0000000000000000e+00].MeshIndices, vector<unsigned int>({ 2, 0 }));
        EXPECT_EQ(result.Segments.size(), 4);
        EXPECT_EQ(result.Segments[0].Points, vector<double>({ 0.0000000000000000e+00, 2.5000000000000000e-01 }));
        EXPECT_EQ(result.Segments[0].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Second);
        EXPECT_EQ(result.Segments[0].MeshIndices, vector<unsigned int>({ 0, 0 }));
        EXPECT_EQ(result.Segments[1].Points, vector<double>({ 2.5000000000000000e-01, 5.0000000000000000e-01 }));
        EXPECT_EQ(result.Segments[1].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[1].MeshIndices, vector<unsigned int>({ 0, 0 }));
        EXPECT_EQ(result.Segments[2].Points, vector<double>({ 5.0000000000000000e-01, 7.5000000000000000e-01 }));
        EXPECT_EQ(result.Segments[2].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::Both);
        EXPECT_EQ(result.Segments[2].MeshIndices, vector<unsigned int>({ 0, 1 }));
        EXPECT_EQ(result.Segments[3].Points, vector<double>({ 7.5000000000000000e-01, 1.0000000000000000e+00 }));
        EXPECT_EQ(result.Segments[3].Type, Gedim::UnionMeshSegment::UnionMesh::UnionMeshSegment::Types::First);
        EXPECT_EQ(result.Segments[3].MeshIndices, vector<unsigned int>({ 1, 0 }));
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
