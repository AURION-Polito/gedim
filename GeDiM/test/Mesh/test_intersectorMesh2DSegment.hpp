#ifndef __TEST_IntersectorMesh2DSegment_H
#define __TEST_IntersectorMesh2DSegment_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "IntersectorMesh2DSegment.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"
#include "MeshMatrices_2D_2Cells_Mock.hpp"
#include "MeshMatrices_2D_4Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"

using namespace testing;

namespace GedimUnitTesting
{
TEST(TestIntersectorMesh2DSegment, TestIntersectMesh)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1e-15;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // segment inside mesh2D
        {
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMesh;
            Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

            const Eigen::Vector3d segmentOrigin(0.25, 0.5, 0.0);
            const Eigen::Vector3d segmentEnd(0.5, 0.25, 0.0);
            const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);
            const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin, segmentEnd);
            const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin, segmentEnd);

            Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
            ASSERT_NO_THROW(
                IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin, segmentEnd, segmentTangent, segmentBarycenter, segmentLength, result));

            EXPECT_EQ(result.Points.size(), 2);
            EXPECT_EQ(result.Points[0.0].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0].Cell2DIds, result.Points[1.0].Cell2DIds);
            EXPECT_EQ(result.Segments.size(), 1);
            EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 1);
        }

        // segmentalong single edge
        {
            GedimUnitTesting::MeshMatrices_2D_2Cells_Mock mockMesh;
            Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

            const Eigen::Vector3d segmentOrigin(0.25, 0.75, 0.0);
            const Eigen::Vector3d segmentEnd(0.75, 0.25, 0.0);

            const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);
            const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin, segmentEnd);
            const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin, segmentEnd);

            Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
            ASSERT_NO_THROW(
                IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin, segmentEnd, segmentTangent, segmentBarycenter, segmentLength, result));

            EXPECT_EQ(result.Points.size(), 2);
            EXPECT_EQ(result.Points[0.0].Vertex2DIds.size(), 0);
            EXPECT_EQ(result.Points[0.0].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Points[0.0].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Points[1.0].Vertex2DIds.size(), 0);
            EXPECT_EQ(result.Points[1.0].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Points[0.0].Edge2DIds, result.Points[1.0].Edge2DIds);
            EXPECT_EQ(result.Points[0.0].Cell2DIds, result.Points[1.0].Cell2DIds);
            EXPECT_EQ(result.Segments.size(), 1);
            EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 2);
        }

        // segmentalong two edges
        {
            GedimUnitTesting::MeshMatrices_2D_4Cells_Mock mockMesh;
            Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

            const Eigen::Vector3d segmentOrigin(0.25, 0.75, 0.0);
            const Eigen::Vector3d segmentEnd(0.75, 0.25, 0.0);

            const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);
            const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin, segmentEnd);
            const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin, segmentEnd);

            Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
            ASSERT_NO_THROW(
                IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin, segmentEnd, segmentTangent, segmentBarycenter, segmentLength, result));

            EXPECT_EQ(result.Points.size(), 3);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Points[4.9999999999999994e-01].Vertex2DIds.size(), 1);
            EXPECT_EQ(result.Points[4.9999999999999994e-01].Edge2DIds.size(), 4);
            EXPECT_EQ(result.Points[4.9999999999999994e-01].Cell2DIds.size(), 4);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Vertex2DIds.size(), 0);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Segments.size(), 2);
            EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Segments[1].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Segments[1].Cell2DIds.size(), 2);
        }

        // segmentalong generic mesh
        {
            GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mockMesh;
            Gedim::MeshMatricesDAO fractureMesh(mockMesh.Mesh);

            const Eigen::Vector3d segmentOrigin(0.25, 0.75, 0.0);
            const Eigen::Vector3d segmentEnd(0.75, 0.25, 0.0);

            const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);
            const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin, segmentEnd);
            const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin, segmentEnd);

            Gedim::IntersectorMesh2DSegment IntersectorMesh2DSegment(fractureMesh, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
            ASSERT_NO_THROW(
                IntersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin, segmentEnd, segmentTangent, segmentBarycenter, segmentLength, result));

            EXPECT_EQ(result.Points.size(), 4);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Vertex2DIds.size(), 1);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Edge2DIds.size(), 7);
            EXPECT_EQ(result.Points[0.0000000000000000e+00].Cell2DIds.size(), 7);
            EXPECT_EQ(result.Points[2.5000000000000006e-01].Vertex2DIds.size(), 0);
            EXPECT_EQ(result.Points[2.5000000000000006e-01].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Points[2.5000000000000006e-01].Cell2DIds.size(), 2);
            EXPECT_EQ(result.Points[5.0000000000000000e-01].Vertex2DIds.size(), 1);
            EXPECT_EQ(result.Points[5.0000000000000000e-01].Edge2DIds.size(), 7);
            EXPECT_EQ(result.Points[5.0000000000000000e-01].Cell2DIds.size(), 7);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Vertex2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Edge2DIds.size(), 8);
            EXPECT_EQ(result.Points[1.0000000000000000e+00].Cell2DIds.size(), 8);
            EXPECT_EQ(result.Segments.size(), 3);
            EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Segments[1].Edge2DIds.size(), 0);
            EXPECT_EQ(result.Segments[1].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Segments[2].Edge2DIds.size(), 1);
            EXPECT_EQ(result.Segments[2].Cell2DIds.size(), 2);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestIntersectorMesh2DSegment, TestIntersectPolygonMesh)
{
    try
    {
        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-08;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        // segment inside single cell
        {
            Eigen::MatrixXd domainVertices(3, 18);
            domainVertices.row(0) << 0.0000000000000000e+00, 4.7870399680095710e-01, 9.2996690625192069e-01,
                1.3020525955672522e+00, 1.5523023989734441e+00, 1.6520254004587949e+00, 1.5971586348660534e+00,
                1.4070017796939815e+00, 1.2295922514998097e+00, 1.0816227865207027e+00, 8.9978717716903311e-01,
                6.6435084777128406e-01, 3.1127319189753433e-01, -1.2835226332822985e-01, -5.6896726300581724e-01,
                -4.2672544725436312e-01, -2.8448363150290878e-01, -1.4224181575145439e-01;
            domainVertices.row(1) << 0.0000000000000000e+00, -1.3877787807814457e-17, 1.5974763121761879e-01,
                4.6092823971113134e-01, 8.6901224456914061e-01, 1.3372032857881035e+00, 1.7872895129852631e+00,
                2.1904991819033635e+00, 2.4024538316062416e+00, 2.5288913940860285e+00, 2.4850252236959895e+00,
                2.4282283719781645e+00, 2.3430516309168499e+00, 2.2369960283014390e+00, 2.1307017071771104e+00,
                1.5980262803828331e+00, 1.0653508535885552e+00, 5.3267542679427760e-01;
            domainVertices.row(2) << 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
                0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

            Gedim::MeshMatrices meshData;
            Gedim::MeshMatricesDAO domainMesh(meshData);

            const unsigned int numDomainVertices = domainVertices.cols();
            std::vector<unsigned int> vertexMarkers(numDomainVertices);
            std::vector<unsigned int> edgeMarkers(numDomainVertices);
            for (unsigned int v = 0; v < numDomainVertices; v++)
            {
                vertexMarkers[v] = v + 1;
                edgeMarkers[v] = numDomainVertices + v + 1;
            }

            Gedim::MeshUtilities meshUtilities;

            meshUtilities.Mesh2DFromPolygon(domainVertices, vertexMarkers, edgeMarkers, domainMesh);
            meshUtilities.ComputeCell1DCell2DNeighbours(domainMesh);

            std::string exportFolder = "./Export/TestIntersectorMesh2DSegment/TestIntersectPolygonMesh";
            Gedim::Output::CreateFolder(exportFolder);

            {
                meshUtilities.ExportMeshToVTU(domainMesh, exportFolder, "Mesh");
            }

            const Eigen::Vector3d segmentOrigin(4.3603822504511447e-01, 1.1740861406549512e+00, 0.0);
            const Eigen::Vector3d segmentEnd(1.5971586348660536e+00, 1.7872895129852631e+00, 0.0);
            const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin, segmentEnd);
            const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin, segmentEnd);
            const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin, segmentEnd);

            {
                Gedim::VTKUtilities exporter;
                exporter.AddSegment(segmentOrigin, segmentEnd);
                exporter.Export(exportFolder + "/Segment.vtu");
            }

            Gedim::IntersectorMesh2DSegment intersectorMesh2DSegment(domainMesh, geometryUtilities);

            Gedim::IntersectorMesh2DSegment::IntersectionMesh result;
            ASSERT_NO_THROW(
                intersectorMesh2DSegment.CreateIntersectionMesh(segmentOrigin, segmentEnd, segmentTangent, segmentBarycenter, segmentLength, result));

            EXPECT_EQ(result.Points.size(), 2);
            EXPECT_EQ(result.Points[0.0].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0].Cell2DIds.size(), 1);
            EXPECT_EQ(result.Points[1.0].Cell2DIds, result.Points[1.0].Cell2DIds);
            EXPECT_EQ(result.Segments.size(), 1);
            EXPECT_EQ(result.Segments[0].Edge2DIds.size(), 0);
            EXPECT_EQ(result.Segments[0].Cell2DIds.size(), 1);
        }
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_IntersectorMesh2DSegment_H
