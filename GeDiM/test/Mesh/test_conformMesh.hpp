#ifndef __TEST_CONFORMMESH_H
#define __TEST_CONFORMMESH_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ConformerMeshPolygon.hpp"
#include "ConformerMeshSegment.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "UnionMeshSegment.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{

TEST(TestConformMesh, TestSerialization)
{
    try
    {
        Gedim::ConformerMeshSegment::ConformMesh mesh;

        mesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(
            0.0,
            Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        mesh.Points.at(0.0).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Types::Inherited;
        mesh.Points.at(0.0).Vertex2DIds = {2, 3};
        mesh.Points.at(0.0).Edge2DIds = {5, 6};
        mesh.Points.at(0.0).Cell2DIds = {8, 9};
        mesh.Points.insert(pair<double, Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint>(
            1.0,
            Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint()));
        mesh.Points.at(1.0).Type = Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Types::Original;
        mesh.Points.at(1.0).Vertex2DIds = {2, 3};
        mesh.Points.at(1.0).Edge2DIds = {5, 6};
        mesh.Points.at(1.0).Cell2DIds = {8, 9};
        mesh.Segments.resize(1);
        mesh.Segments.at(0).Points = {0.0, 1.0};
        mesh.Segments.at(0).Edge2DIds = {5, 6};
        mesh.Segments.at(0).Cell2DIds = {8, 9};

        std::stringstream buffer;
        ASSERT_NO_THROW(Gedim::ConformerMeshSegment::Serialize(buffer, mesh));
        Gedim::ConformerMeshSegment::ConformMesh result;
        ASSERT_NO_THROW(Gedim::ConformerMeshSegment::Deserialize(buffer, result));

        ASSERT_EQ(mesh.Points.size(), result.Points.size());
        ASSERT_EQ(mesh.Points.at(0.0).Type, result.Points.at(0.0).Type);
        ASSERT_EQ(mesh.Points.at(0.0).Vertex2DIds, result.Points.at(0.0).Vertex2DIds);
        ASSERT_EQ(mesh.Points.at(0.0).Edge2DIds, result.Points.at(0.0).Edge2DIds);
        ASSERT_EQ(mesh.Points.at(0.0).Cell2DIds, result.Points.at(0.0).Cell2DIds);
        ASSERT_EQ(mesh.Points.at(1.0).Type, result.Points.at(1.0).Type);
        ASSERT_EQ(mesh.Points.at(1.0).Vertex2DIds, result.Points.at(1.0).Vertex2DIds);
        ASSERT_EQ(mesh.Points.at(1.0).Edge2DIds, result.Points.at(1.0).Edge2DIds);
        ASSERT_EQ(mesh.Points.at(1.0).Cell2DIds, result.Points.at(1.0).Cell2DIds);
        ASSERT_EQ(mesh.Segments.size(), result.Segments.size());
        ASSERT_EQ(mesh.Segments.at(0).Edge2DIds, result.Segments.at(0).Edge2DIds);
        ASSERT_EQ(mesh.Segments.at(0).Cell2DIds, result.Segments.at(0).Cell2DIds);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh1D)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // conform of simple two points mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO domainMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.25, 0.5, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.5, 0.25, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(domainMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.25, 0.5, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.5, 0.25, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            EXPECT_EQ(conformMeshOne.Points.size(), 2);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[0].Cell2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            EXPECT_EQ(conformMeshTwo.Points.size(), 2);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[0].Cell2DIds.size(), 1);
        }

        // conform of intricate two mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.25, 0.5, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.5, 0.25, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            Eigen::Vector3d segmentTwoOrigin(0.25, 0.5, 0.0);
            Eigen::Vector3d segmentTwoEnd(0.5, 0.25, 0.0);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            EXPECT_EQ(conformMeshOne.Points.size(), 5);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[1.2500000000000003e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Inherited);
            EXPECT_EQ(conformMeshOne.Points[1.2500000000000003e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.2500000000000003e-01].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.2500000000000003e-01].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[3.3333333333333337e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Inherited);
            EXPECT_EQ(conformMeshOne.Points[3.3333333333333337e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[3.3333333333333337e-01].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[3.3333333333333337e-01].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments.size(), 4);
            EXPECT_EQ(conformMeshOne.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[0].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[1].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[1].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[2].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[2].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[3].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[3].Cell2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            EXPECT_EQ(conformMeshTwo.Points.size(), 5);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.2500000000000003e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[1.2500000000000003e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[1.2500000000000003e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.2500000000000003e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(conformMeshTwo.Points[3.3333333333333337e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[3.3333333333333337e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[3.3333333333333337e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[3.3333333333333337e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Inherited);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Vertex2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Edge2DIds.size(), 5);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Cell2DIds.size(), 5);
            EXPECT_EQ(conformMeshTwo.Segments.size(), 4);
            EXPECT_EQ(conformMeshTwo.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[0].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[1].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[1].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[2].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[2].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[3].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[3].Cell2DIds.size(), 1);
        }

        // add points to conform mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.25, 0.25, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.75, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.25, 0.25, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.75, 0.75, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            ASSERT_NO_THROW(conformMeshSegmentOne.InsertExternalPoint(segmentOneOrigin, segmentOneEnd, fractureMeshOne, 0.25, conformMeshOne));

            EXPECT_EQ(conformMeshOne.Points.size(), 4);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[2.5000000000000000e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External);
            EXPECT_EQ(conformMeshOne.Points[2.5000000000000000e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[2.5000000000000000e-01].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[2.5000000000000000e-01].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Points[4.9999999999999994e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Points[1.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments.size(), 3);
            EXPECT_EQ(conformMeshOne.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[0].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[1].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[1].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshOne.Segments[2].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshOne.Segments[2].Cell2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            ASSERT_NO_THROW(conformMeshSegmentTwo.InsertExternalPoint(segmentTwoOrigin, segmentTwoEnd, fractureMeshTwo, 0.75, conformMeshTwo));

            EXPECT_EQ(conformMeshTwo.Points.size(), 4);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[0.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[4.9999999999999994e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(conformMeshTwo.Points[7.5000000000000000e-01].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::External);
            EXPECT_EQ(conformMeshTwo.Points[7.5000000000000000e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[7.5000000000000000e-01].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[7.5000000000000000e-01].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Type,
                      Gedim::ConformerMeshSegment::ConformMesh::ConformMeshPoint::Original);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Points[1.0000000000000000e+00].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments.size(), 3);
            EXPECT_EQ(conformMeshTwo.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[0].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[1].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[1].Cell2DIds.size(), 1);
            EXPECT_EQ(conformMeshTwo.Segments[2].Edge2DIds.size(), 0);
            EXPECT_EQ(conformMeshTwo.Segments[2].Cell2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DTwoPoints)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        // conform of simple 2 points mesh
        {
            std::string exportFolder = "./Export/TestConformMesh2D/TestConformMesh2DTwoPoints";
            Gedim::Output::CreateFolder(exportFolder);
            Gedim::MeshUtilities meshUtilities;

            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "DomainOne_Original");

            const Eigen::Vector3d segmentOneOrigin(0.75, 0.0, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.0, 0.75, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.75, 0.0, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.0, 0.75, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "DomainOne_Conformed");

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 6);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 10);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 4);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshTwo);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 6);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 10);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 4);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DFivePoints)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // conform of simple 5 points mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            std::string exportFolder = "./Export/TestConformMesh2D/TestFivePointsMesh";
            Gedim::Output::CreateFolder(exportFolder);

            Gedim::MeshUtilities meshUtilities;
            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "Original");

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "Conformed");

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 8);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 16);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 9);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[1], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[3], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[6], false);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshTwo);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 8);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 16);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 9);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[1], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[3], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[6], false);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DTwentysixPoints)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        // conform of 26 points mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            std::string exportFolder = "./Export/TestConformMesh2D/TestTwentysixPointsMesh";
            Gedim::Output::CreateFolder(exportFolder);

            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "Original");

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            meshUtilities.ExportMeshToVTU(fractureMeshOne, exportFolder, "Conformed");

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 31);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 76);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 53);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[10], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[19], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[3], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[17], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[14], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[5], false);
            EXPECT_EQ(mockMeshOne.Mesh.ActiveCell2D[20], false);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshTwo);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 31);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 76);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 53);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[10], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[19], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[3], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[17], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[14], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[5], false);
            EXPECT_EQ(mockMeshTwo.Mesh.ActiveCell2D[20], false);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DTwoPointsFourPoints)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        // conform of simple 2 points and 4 points mesh
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 7);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 12);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 6);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 8);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 16);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 9);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DTwentysixPointsTwentysixPoints)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        // conform of simple 26 points and 26 points mesh with different traces
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.0, 0.75, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.0, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(1.0, 0.25, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.25, 1.0, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 33);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 80);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 57);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 32);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 76);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 53);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DTwoPointsWithNonPassingTrace)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-14;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        // conform of simple 2 points with trace non-passing by domain edges
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.25, 0.25, 0.0);
            const Eigen::Vector3d segmentOneEnd(0.75, 0.75, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(0.75, 0.75, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.25, 0.25, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            Gedim::MeshUtilities::CheckMesh2DConfiguration config;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshOne);

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 7);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 13);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 11);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            config.Cell0D_CheckCoordinates2D = true;
            config.Cell0D_CheckDuplications = true;
            config.Cell1D_CheckDuplications = true;
            config.Cell1D_CheckNeighbours = true;
            config.Cell2D_CheckEdges = true;
            config.Cell2D_CheckDuplications = true;
            meshUtilities.CheckMesh2D(config, geometryUtilities, fractureMeshTwo);

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 7);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 13);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 11);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_GE(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMesh, TestConformMesh2DOnlyOnEdges)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // conform of simple 2 points mesh only on edges
        {
            // Create mesh fracture one and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMeshOne;
            Gedim::MeshMatricesDAO fractureMeshOne(mockMeshOne.Mesh);

            const Eigen::Vector3d segmentOneOrigin(0.0, 1.0, 0.0);
            const Eigen::Vector3d segmentOneEnd(1.0, 0.0, 0.0);

            const Eigen::Vector3d segmentOneTangent = geometryUtilities.SegmentTangent(segmentOneOrigin, segmentOneEnd);
            const Eigen::Vector3d segmentOneBarycenter = geometryUtilities.SegmentBarycenter(segmentOneOrigin, segmentOneEnd);
            const double segmentOneLength = geometryUtilities.SegmentLength(segmentOneOrigin, segmentOneEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentOne(fractureMeshOne, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshOne;
            intersectorMeshSegmentOne.CreateIntersectionMesh(segmentOneOrigin,
                                                             segmentOneEnd,
                                                             segmentOneTangent,
                                                             segmentOneBarycenter,
                                                             segmentOneLength,
                                                             intersectionMeshOne);

            // Create mesh fracture two and corresponding intersection meshes
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMeshTwo;
            Gedim::MeshMatricesDAO fractureMeshTwo(mockMeshTwo.Mesh);

            const Eigen::Vector3d segmentTwoOrigin(1.0, 0.0, 0.0);
            const Eigen::Vector3d segmentTwoEnd(0.0, 1.0, 0.0);

            const Eigen::Vector3d segmentTwoTangent = geometryUtilities.SegmentTangent(segmentTwoOrigin, segmentTwoEnd);
            const Eigen::Vector3d segmentTwoBarycenter = geometryUtilities.SegmentBarycenter(segmentTwoOrigin, segmentTwoEnd);
            const double segmentTwoLength = geometryUtilities.SegmentLength(segmentTwoOrigin, segmentTwoEnd);

            Gedim::IntersectorMesh2DSegment intersectorMeshSegmentTwo(fractureMeshTwo, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh intersectionMeshTwo;
            intersectorMeshSegmentTwo.CreateIntersectionMesh(segmentTwoOrigin,
                                                             segmentTwoEnd,
                                                             segmentTwoTangent,
                                                             segmentTwoBarycenter,
                                                             segmentTwoLength,
                                                             intersectionMeshTwo);

            // Create mesh union
            vector<double> curvilinearCoordinatesMeshOne;
            vector<double> curvilinearCoordinatesMeshTwo;
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshOne, curvilinearCoordinatesMeshOne);
            Gedim::IntersectorMesh2DSegment::ToCurvilinearCoordinates(intersectionMeshTwo, curvilinearCoordinatesMeshTwo);
            Gedim::UnionMeshSegment unionMeshSegment(geometryUtilities);

            Gedim::UnionMeshSegment::UnionMesh unionMesh;
            unionMeshSegment.CreateUnionMesh(curvilinearCoordinatesMeshOne, curvilinearCoordinatesMeshTwo, unionMesh);

            // Create first conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentOne(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshOne;
            ASSERT_NO_THROW(conformMeshSegmentOne.CreateConformMesh(intersectionMeshOne, unionMesh, 0, conformMeshOne));

            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration conformerMeshPolygonOneConfig;
            conformerMeshPolygonOneConfig.Type = Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges;
            Gedim::ConformerMeshPolygon conformerMeshPolygonOne(geometryUtilities, conformerMeshPolygonOneConfig);

            ASSERT_NO_THROW(conformerMeshPolygonOne.CreateConformMesh(segmentOneOrigin, segmentOneEnd, segmentOneTangent, conformMeshOne, fractureMeshOne));

            EXPECT_EQ(mockMeshOne.Mesh.NumberCell0D, 5);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell1D, 7);
            EXPECT_EQ(mockMeshOne.Mesh.NumberCell2D, 4);
            for (const auto &mesh1Dpoint : conformMeshOne.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshOne.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 2);

            // Create second conform mesh
            Gedim::ConformerMeshSegment conformMeshSegmentTwo(geometryUtilities);

            Gedim::ConformerMeshSegment::ConformMesh conformMeshTwo;
            ASSERT_NO_THROW(conformMeshSegmentTwo.CreateConformMesh(intersectionMeshTwo, unionMesh, 1, conformMeshTwo));

            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration conformerMeshPolygonTwoConfig;
            conformerMeshPolygonTwoConfig.Type = Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges;
            Gedim::ConformerMeshPolygon conformerMeshPolygonTwo(geometryUtilities, conformerMeshPolygonTwoConfig);

            ASSERT_NO_THROW(conformerMeshPolygonTwo.CreateConformMesh(segmentTwoOrigin, segmentTwoEnd, segmentTwoTangent, conformMeshTwo, fractureMeshTwo));

            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell0D, 5);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell1D, 8);
            EXPECT_EQ(mockMeshTwo.Mesh.NumberCell2D, 4);
            for (const auto &mesh1Dpoint : conformMeshTwo.Points)
                EXPECT_EQ(mesh1Dpoint.second.Vertex2DIds.size(), 1);
            for (const auto &mesh1Dsegment : conformMeshTwo.Segments)
                EXPECT_EQ(mesh1Dsegment.Edge2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_CONFORMMESH_H
