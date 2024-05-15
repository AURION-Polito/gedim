#ifndef __TEST_IntersectorMesh3DSegment_H
#define __TEST_IntersectorMesh3DSegment_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "GeometryUtilities.hpp"
#include "MeshMatrices_3D_68Cells_Mock.hpp"
#include "MeshUtilities.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_3D_1Cells_Mock.hpp"
#include "MeshMatrices_3D_22Cells_Mock.hpp"
#include "MeshMatrices_3D_329Cells_Mock.hpp"
#include "IntersectorMesh3DSegment.hpp"
#include "VTKUtilities.hpp"

using namespace testing;

namespace GedimUnitTesting
{
  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentFullInsideOneCell)
  {
    std::string exportFolder = "./Export/TestIntersectorMesh3DSegment/TestIntersectMesh_SegmentInside";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    // segment inside mesh3D
    GedimUnitTesting::MeshMatrices_3D_1Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh3D(mockMesh.Mesh);

    const Gedim::MeshUtilities::MeshGeometricData3D mesh3D_geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                                 mesh3D);

    const Eigen::Vector3d segmentOrigin(0.25, 0.5, 0.75);
    const Eigen::Vector3d segmentEnd(0.5, 0.25, 0.25);
    const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                            segmentEnd);
    const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                  segmentEnd);
    const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                 segmentEnd);

    {
      meshUtilities.ExportMeshToVTU(mesh3D,
                                    exportFolder,
                                    "mesh3D");

      {
        Gedim::VTKUtilities exporter;
        exporter.AddSegment(segmentOrigin,
                            segmentEnd);

        exporter.Export(exportFolder + "/segment.vtu");
      }
    }

    const Gedim::IntersectorMesh3DSegment intersectorMesh3DSegment(geometryUtilities,
                                                                   meshUtilities);

    const Gedim::IntersectorMesh3DSegment::IntersectionMesh intersections =
        intersectorMesh3DSegment.CreateIntersectionMesh(segmentOrigin,
                                                        segmentEnd,
                                                        segmentTangent,
                                                        segmentBarycenter,
                                                        segmentLength,
                                                        mesh3D,
                                                        mesh3D_geometricData);


    EXPECT_EQ(intersections.Points.size(), 2);
    EXPECT_EQ(intersections.Points[0].CurvilinearCoordinate, 0.0);
    EXPECT_EQ(intersections.Points[1].CurvilinearCoordinate, 1.0);
    EXPECT_EQ(intersections.Points[0].Cell3DIds.size(), 1);
    EXPECT_EQ(intersections.Points[1].Cell3DIds.size(), 1);
    EXPECT_EQ(intersections.Points[0].Cell3DIds[0], 0);
    EXPECT_EQ(intersections.Points[1].Cell3DIds[0], 0);
    EXPECT_EQ(intersections.Segments.size(), 1);
    EXPECT_EQ(intersections.Segments[0].Cell3DIds.size(), 1);
    EXPECT_EQ(intersections.Segments[0].Cell3DIds[0], 0);

    Gedim::MeshMatrices mesh_1D_data;
    Gedim::MeshMatricesDAO mesh_1D(mesh_1D_data);

    const std::vector<double> coordinates = Gedim::IntersectorMesh3DSegment::ToCurvilinearCoordinates(intersections);

    meshUtilities.FillMesh1D(geometryUtilities,
                             segmentOrigin,
                             segmentTangent,
                             coordinates,
                             mesh_1D);
    {
      meshUtilities.ExportMeshToVTU(mesh_1D,
                                    exportFolder,
                                    "mesh1D");
    }
  }

  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentOnFace)
  {
    std::string exportFolder = "./Export/TestIntersectorMesh3DSegment/TestIntersectMesh_SegmentOnFace";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    // segment inside mesh3D
    GedimUnitTesting::MeshMatrices_3D_22Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh3D(mockMesh.Mesh);
    meshUtilities.ComputeCell2DCell3DNeighbours(mesh3D);

    const Gedim::MeshUtilities::MeshGeometricData3D mesh3D_geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                                 mesh3D);

    const Eigen::Vector3d segmentOrigin(0.25, 0.5, 0.75);
    const Eigen::Vector3d segmentEnd(0.5, 0.25, 0.25);
    const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                            segmentEnd);
    const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                  segmentEnd);
    const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                 segmentEnd);

    {
      meshUtilities.ExportMeshToVTU(mesh3D,
                                    exportFolder,
                                    "mesh3D");

      {
        Gedim::VTKUtilities exporter;
        exporter.AddSegment(segmentOrigin,
                            segmentEnd);

        exporter.Export(exportFolder + "/segment.vtu");
      }
    }

    const Gedim::IntersectorMesh3DSegment intersectorMesh3DSegment(geometryUtilities,
                                                                   meshUtilities);

    const Gedim::IntersectorMesh3DSegment::IntersectionMesh intersections =
        intersectorMesh3DSegment.CreateIntersectionMesh(segmentOrigin,
                                                        segmentEnd,
                                                        segmentTangent,
                                                        segmentBarycenter,
                                                        segmentLength,
                                                        mesh3D,
                                                        mesh3D_geometricData);

    EXPECT_EQ(intersections.Points.size(), 2);
    EXPECT_EQ(intersections.Points[0].CurvilinearCoordinate, 0.0);
    EXPECT_EQ(intersections.Points[1].CurvilinearCoordinate, 1.0);
    EXPECT_EQ(intersections.Points[0].Cell3DIds, std::vector<unsigned int>({ 6, 11 }));
    EXPECT_EQ(intersections.Points[1].Cell3DIds, std::vector<unsigned int>({ 0,2,6,9,11,12 }));
    EXPECT_EQ(intersections.Segments.size(), 1);
    EXPECT_EQ(intersections.Segments[0].Cell3DIds, std::vector<unsigned int>({ 6, 11 }));

    Gedim::MeshMatrices mesh_1D_data;
    Gedim::MeshMatricesDAO mesh_1D(mesh_1D_data);

    const std::vector<double> coordinates = Gedim::IntersectorMesh3DSegment::ToCurvilinearCoordinates(intersections);

    meshUtilities.FillMesh1D(geometryUtilities,
                             segmentOrigin,
                             segmentTangent,
                             coordinates,
                             mesh_1D);
    {
      meshUtilities.ExportMeshToVTU(mesh_1D,
                                    exportFolder,
                                    "mesh1D");
    }
  }

  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentInsideOneCell)
  {
    std::string exportFolder = "./Export/TestIntersectorMesh3DSegment/TestIntersectMesh_SegmentInsideOneCell";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    // segment inside mesh3D
    GedimUnitTesting::MeshMatrices_3D_68Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh3D(mockMesh.Mesh);
    meshUtilities.ComputeCell2DCell3DNeighbours(mesh3D);

    const Gedim::MeshUtilities::MeshGeometricData3D mesh3D_geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                                 mesh3D);

    const Eigen::Vector3d segmentOrigin(0.25, 0.5, 0.75);
    const Eigen::Vector3d segmentEnd(0.5, 0.25, 0.25);
    const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                            segmentEnd);
    const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                  segmentEnd);
    const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                 segmentEnd);

    {
      meshUtilities.ExportMeshToVTU(mesh3D,
                                    exportFolder,
                                    "mesh3D");

      {
        Gedim::VTKUtilities exporter;
        exporter.AddSegment(segmentOrigin,
                            segmentEnd);

        exporter.Export(exportFolder + "/segment.vtu");
      }
    }

    const Gedim::IntersectorMesh3DSegment intersectorMesh3DSegment(geometryUtilities,
                                                                   meshUtilities);

    const Gedim::IntersectorMesh3DSegment::IntersectionMesh intersections =
        intersectorMesh3DSegment.CreateIntersectionMesh(segmentOrigin,
                                                        segmentEnd,
                                                        segmentTangent,
                                                        segmentBarycenter,
                                                        segmentLength,
                                                        mesh3D,
                                                        mesh3D_geometricData);
    EXPECT_EQ(intersections.Points.size(), 2);
    EXPECT_EQ(intersections.Points[0].CurvilinearCoordinate, 0.0);
    EXPECT_EQ(intersections.Points[1].CurvilinearCoordinate, 1.0);
    EXPECT_EQ(intersections.Points[0].Cell3DIds, std::vector<unsigned int>({ 0,18,21,25,40,60,62 }));
    EXPECT_EQ(intersections.Points[1].Cell3DIds, std::vector<unsigned int>({ 1,17,27,50,52,62 }));
    EXPECT_EQ(intersections.Segments.size(), 1);
    EXPECT_EQ(intersections.Segments[0].Cell3DIds, std::vector<unsigned int>({ 62 }));

    Gedim::MeshMatrices mesh_1D_data;
    Gedim::MeshMatricesDAO mesh_1D(mesh_1D_data);

    const std::vector<double> coordinates = Gedim::IntersectorMesh3DSegment::ToCurvilinearCoordinates(intersections);

    meshUtilities.FillMesh1D(geometryUtilities,
                             segmentOrigin,
                             segmentTangent,
                             coordinates,
                             mesh_1D);
    {
      meshUtilities.ExportMeshToVTU(mesh_1D,
                                    exportFolder,
                                    "mesh1D");
    }
  }

  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentOnEdge)
  {
    std::string exportFolder = "./Export/TestIntersectorMesh3DSegment/TestIntersectMesh_SegmentOnEdge";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-8;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    // segment inside mesh3D
    GedimUnitTesting::MeshMatrices_3D_329Cells_Mock mockMesh;
    Gedim::MeshMatricesDAO mesh3D(mockMesh.Mesh);
    meshUtilities.ComputeCell2DCell3DNeighbours(mesh3D);

    const Gedim::MeshUtilities::MeshGeometricData3D mesh3D_geometricData = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                                                                 mesh3D);

    const Eigen::Vector3d segment_test_Origin(3.1628787878787878e-01,
                                              3.1628787878787878e-01,
                                              3.1628787878787878e-01);
    const Eigen::Vector3d segment_test_End(2.9578836930455638e-01,
                                           6.7768285371702630e-01,
                                           3.2231714628297359e-01);
    const Eigen::Vector3d segment_test_Tangent = geometryUtilities.SegmentTangent(segment_test_Origin,
                                                                                  segment_test_End);

    const Eigen::Vector3d segmentOrigin = segment_test_Origin - 0.5 * segment_test_Tangent;
    const Eigen::Vector3d segmentEnd = segment_test_Origin + 1.1 * segment_test_Tangent;
    const Eigen::Vector3d segmentTangent = geometryUtilities.SegmentTangent(segmentOrigin,
                                                                            segmentEnd);
    const Eigen::Vector3d segmentBarycenter = geometryUtilities.SegmentBarycenter(segmentOrigin,
                                                                                  segmentEnd);
    const double segmentLength = geometryUtilities.SegmentLength(segmentOrigin,
                                                                 segmentEnd);

    {
      meshUtilities.ExportMeshToVTU(mesh3D,
                                    exportFolder,
                                    "mesh3D");

      {
        Gedim::VTKUtilities exporter;
        exporter.AddSegment(segmentOrigin,
                            segmentEnd);

        exporter.Export(exportFolder + "/segment.vtu");
      }
    }

    const Gedim::IntersectorMesh3DSegment intersectorMesh3DSegment(geometryUtilities,
                                                                   meshUtilities);

    const Gedim::IntersectorMesh3DSegment::IntersectionMesh intersections =
        intersectorMesh3DSegment.CreateIntersectionMesh(segmentOrigin,
                                                        segmentEnd,
                                                        segmentTangent,
                                                        segmentBarycenter,
                                                        segmentLength,
                                                        mesh3D,
                                                        mesh3D_geometricData);

    std::cout<< Gedim::IntersectorMesh3DSegment::ToString(intersections)<< std::endl;

    EXPECT_EQ(intersections.Points.size(), 4);
    EXPECT_EQ(intersections.Points[0].CurvilinearCoordinate, 0.0000000000000000e+00);
    EXPECT_EQ(intersections.Points[1].CurvilinearCoordinate, 3.1249999999999967e-01);
    EXPECT_EQ(intersections.Points[2].CurvilinearCoordinate, 9.3749999999999989e-01);
    EXPECT_EQ(intersections.Points[3].CurvilinearCoordinate, 1.0000000000000000e+00);
    EXPECT_EQ(intersections.Points[0].Cell3DIds, std::vector<unsigned int>({ 81 }));
    EXPECT_EQ(intersections.Points[1].Cell3DIds, std::vector<unsigned int>({ 12,59,69,81,111,116,148,152,156,168,196,225,246,250,262,263,270,276,279,280,281,288,298,305,311,312,327,328 }));
    EXPECT_EQ(intersections.Points[2].Cell3DIds, std::vector<unsigned int>({ 1,127,141,154,173,181,188,191,200,231,246,248,250,272,273,277,278,280,281,282,298,300,303,306,314,315,317,323 }));
    EXPECT_EQ(intersections.Points[3].Cell3DIds, std::vector<unsigned int>({ 323 }));
    EXPECT_EQ(intersections.Segments.size(), 3);
    EXPECT_EQ(intersections.Segments[0].Cell3DIds, std::vector<unsigned int>({ 81 }));
    EXPECT_EQ(intersections.Segments[1].Cell3DIds, std::vector<unsigned int>({ 246,250,280,281,298 }));
    EXPECT_EQ(intersections.Segments[2].Cell3DIds, std::vector<unsigned int>({ 323 }));

    Gedim::MeshMatrices mesh_1D_data;
    Gedim::MeshMatricesDAO mesh_1D(mesh_1D_data);

    const std::vector<double> coordinates = Gedim::IntersectorMesh3DSegment::ToCurvilinearCoordinates(intersections);

    meshUtilities.FillMesh1D(geometryUtilities,
                             segmentOrigin,
                             segmentTangent,
                             coordinates,
                             mesh_1D);
    {
      meshUtilities.ExportMeshToVTU(mesh_1D,
                                    exportFolder,
                                    "mesh1D");
    }
  }
}

#endif // __TEST_IntersectorMesh3DSegment_H
