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
#include "IntersectorMesh3DSegment.hpp"
#include "VTKUtilities.hpp"

using namespace testing;

namespace GedimUnitTesting
{
  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentInside)
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

  TEST(TestIntersectorMesh3DSegment, TestIntersectMesh_SegmentOnVertex)
  {
    std::string exportFolder = "./Export/TestIntersectorMesh3DSegment/TestIntersectMesh_SegmentOnVertex";
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

    std::cout<< Gedim::IntersectorMesh3DSegment::ToString(intersections)<< std::endl;

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
}

#endif // __TEST_IntersectorMesh3DSegment_H
