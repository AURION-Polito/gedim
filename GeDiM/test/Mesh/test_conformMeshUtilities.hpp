#ifndef __TEST_ConformMeshUtilities_H
#define __TEST_ConformMeshUtilities_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "GeometryUtilities.hpp"
#include "MeshDAOExporterToCsv.hpp"
#include "MeshUtilities.hpp"

#include "MeshMatrices.hpp"
#include "MeshMatricesDAO.hpp"
#include "MeshMatrices_2D_26Cells_Mock.hpp"

#include "VTKUtilities.hpp"

#include "ConformMeshUtilities.hpp"

using namespace testing;
using namespace std;

namespace GedimUnitTesting
{
TEST(TestConformMeshUtilities, TestConformMeshTwoSegments)
{
    try
    {
        string exportFolder = "./Export";
        Gedim::Output::CreateFolder(exportFolder);
        exportFolder = "./Export/TestConformMeshUtilities/TestConformMeshTwoSegments";
        Gedim::Output::CreateFolder(exportFolder);

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
        Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

        meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ConformedMesh");

        Gedim::ConformMeshUtilities conformMeshUtilities(geometryUtilities, meshUtilities);

        const unsigned int numSegments = 2;

        std::vector<Eigen::MatrixXd> segmentsVertices(numSegments, Eigen::MatrixXd(3, 2));
        std::vector<Eigen::Vector3d> segmentsTangent(numSegments);
        std::vector<Eigen::Vector3d> segmentsBarycenter(numSegments);
        std::vector<double> segmentsLength(numSegments);
        std::vector<double> segmentsSquaredLength(numSegments);

        segmentsVertices[0].col(0) << 0.1, 0.1, 0.0;
        segmentsVertices[0].col(1) << 0.9, 0.7, 0.0;
        segmentsVertices[1].col(0) << 0.75, 0.0, 0.0;
        segmentsVertices[1].col(1) << 0.4, 0.8, 0.0;

        // export segments
        {
            Gedim::VTKUtilities exporter;
            for (unsigned int s = 0; s < numSegments; s++)
                exporter.AddSegment(segmentsVertices[s]);
            exporter.Export(exportFolder + "/Segments.vtu");
        }

        // compute segment geometric properties
        for (unsigned int s = 0; s < numSegments; s++)
        {
            segmentsTangent[s] = geometryUtilities.SegmentTangent(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsBarycenter[s] =
                geometryUtilities.SegmentBarycenter(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsLength[s] = geometryUtilities.SegmentLength(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsSquaredLength[s] = segmentsLength[s] * segmentsLength[s];
        }

        // compute segment intersections
        const vector<list<double>> segmentsAdditionalPoints =
            geometryUtilities.IntersectionsBetweenSegments(segmentsVertices, segmentsTangent, segmentsBarycenter, segmentsLength);

        std::vector<Gedim::IntersectorMesh2DSegment::IntersectionMesh> segmentsIntersectionMesh(numSegments);
        std::vector<std::vector<double>> segmentsCurvilinearCoordinatesMesh(numSegments);
        std::vector<Gedim::UnionMeshSegment::UnionMesh> segmentsUnionMesh(numSegments);
        std::vector<Gedim::ConformerMeshSegment::ConformMesh> segmentsConformMesh(numSegments);
        Gedim::ConformMeshUtilities::ComputeDomainConformedMeshOptions options;
        options.PrintStatus = false;
        options.VtkExportFolder = exportFolder;

        conformMeshUtilities.ComputeConformedMeshWithSegments(
            segmentsAdditionalPoints,
            segmentsVertices,
            segmentsTangent,
            segmentsBarycenter,
            segmentsLength,
            segmentsSquaredLength,
            meshDAO,
            segmentsIntersectionMesh,
            segmentsCurvilinearCoordinatesMesh,
            segmentsUnionMesh,
            segmentsConformMesh,
            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::Generalized,
            options);

        segmentsUnionMesh.clear();
        segmentsIntersectionMesh.clear();
        segmentsCurvilinearCoordinatesMesh.clear();
        segmentsConformMesh.clear();

        segmentsUnionMesh.resize(numSegments);
        segmentsIntersectionMesh.resize(numSegments);
        segmentsCurvilinearCoordinatesMesh.resize(numSegments);
        segmentsConformMesh.resize(numSegments);

        conformMeshUtilities.ComputeConformedMeshWithSegments(
            vector<list<double>>(numSegments),
            segmentsVertices,
            segmentsTangent,
            segmentsBarycenter,
            segmentsLength,
            segmentsSquaredLength,
            meshDAO,
            segmentsIntersectionMesh,
            segmentsCurvilinearCoordinatesMesh,
            segmentsUnionMesh,
            segmentsConformMesh,
            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges,
            options);

        // check meshes 1D
        for (unsigned int i = 0; i < numSegments; i++)
        {
            for (const auto &point : segmentsConformMesh[i].Points)
                Gedim::Output::Assert(point.second.Vertex2DIds.size() == 1);

            for (const auto &segment : segmentsConformMesh[i].Segments)
                Gedim::Output::Assert(segment.Edge2DIds.size() == 1);
        }

        // check meshes 2D
        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        meshUtilities.CheckMesh2D(config, geometryUtilities, meshDAO);

        conformMeshUtilities.AddConformedMeshProperties(segmentsConformMesh, meshDAO);

        ASSERT_EQ(50, mesh.Mesh.NumberCell0D);
        ASSERT_EQ(103, mesh.Mesh.NumberCell1D);
        ASSERT_EQ(54, mesh.Mesh.NumberCell2D);

        meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ConformedMesh");
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}

TEST(TestConformMeshUtilities, TestConformMeshRandomSegments)
{
    try
    {
        string exportFolder = "./Export";
        Gedim::Output::CreateFolder(exportFolder);
        exportFolder = "./Export/TestConformMeshUtilities/TestConformMeshRandomSegments";
        Gedim::Output::CreateFolder(exportFolder);

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = 1.0e-12;
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);
        Gedim::MeshUtilities meshUtilities;

        GedimUnitTesting::MeshMatrices_2D_26Cells_Mock mesh;
        Gedim::MeshMatricesDAO meshDAO(mesh.Mesh);

        meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ConformedMesh");

        Gedim::ConformMeshUtilities conformMeshUtilities(geometryUtilities, meshUtilities);

        const unsigned int numSegments = 10;

        std::vector<Eigen::MatrixXd> segmentsVertices(numSegments, Eigen::MatrixXd(3, 2));
        std::vector<Eigen::Vector3d> segmentsTangent(numSegments);
        std::vector<Eigen::Vector3d> segmentsBarycenter(numSegments);
        std::vector<double> segmentsLength(numSegments);
        std::vector<double> segmentsSquaredLength(numSegments);

        // random segments
        srand(112233);
        for (unsigned int s = 0; s < numSegments; s++)
        {
            segmentsVertices[s].col(0) << rand() / double(RAND_MAX), rand() / double(RAND_MAX), 0.0;
            segmentsVertices[s].col(1) << rand() / double(RAND_MAX), rand() / double(RAND_MAX), 0.0;
        }

        // export segments
        {
            Gedim::VTKUtilities exporter;
            for (unsigned int s = 0; s < numSegments; s++)
            {
                vector<double> id(1, s);
                exporter.AddSegment(segmentsVertices[s],
                                    {{"Id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(id.size()), id.data()}});
            }
            exporter.Export(exportFolder + "/Segments.vtu");
        }

        // compute segment geometric properties
        for (unsigned int s = 0; s < numSegments; s++)
        {
            segmentsTangent[s] = geometryUtilities.SegmentTangent(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsBarycenter[s] =
                geometryUtilities.SegmentBarycenter(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsLength[s] = geometryUtilities.SegmentLength(segmentsVertices.at(s).col(0), segmentsVertices.at(s).col(1));
            segmentsSquaredLength[s] = segmentsLength[s] * segmentsLength[s];
        }

        // compute segment intersections
        const vector<list<double>> segmentsAdditionalPoints =
            geometryUtilities.IntersectionsBetweenSegments(segmentsVertices, segmentsTangent, segmentsBarycenter, segmentsLength);

        std::vector<Gedim::IntersectorMesh2DSegment::IntersectionMesh> segmentsIntersectionMesh(numSegments);
        std::vector<std::vector<double>> segmentsCurvilinearCoordinatesMesh(numSegments);
        std::vector<Gedim::UnionMeshSegment::UnionMesh> segmentsUnionMesh(numSegments);
        std::vector<Gedim::ConformerMeshSegment::ConformMesh> segmentsConformMesh(numSegments);
        Gedim::ConformMeshUtilities::ComputeDomainConformedMeshOptions options;
        options.PrintStatus = false;
        options.VtkExportFolder = exportFolder;

        conformMeshUtilities.ComputeConformedMeshWithSegments(
            segmentsAdditionalPoints,
            segmentsVertices,
            segmentsTangent,
            segmentsBarycenter,
            segmentsLength,
            segmentsSquaredLength,
            meshDAO,
            segmentsIntersectionMesh,
            segmentsCurvilinearCoordinatesMesh,
            segmentsUnionMesh,
            segmentsConformMesh,
            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::Generalized,
            options);

        segmentsUnionMesh.clear();
        segmentsIntersectionMesh.clear();
        segmentsCurvilinearCoordinatesMesh.clear();
        segmentsConformMesh.clear();

        segmentsUnionMesh.resize(numSegments);
        segmentsIntersectionMesh.resize(numSegments);
        segmentsCurvilinearCoordinatesMesh.resize(numSegments);
        segmentsConformMesh.resize(numSegments);

        conformMeshUtilities.ComputeConformedMeshWithSegments(
            vector<list<double>>(numSegments),
            segmentsVertices,
            segmentsTangent,
            segmentsBarycenter,
            segmentsLength,
            segmentsSquaredLength,
            meshDAO,
            segmentsIntersectionMesh,
            segmentsCurvilinearCoordinatesMesh,
            segmentsUnionMesh,
            segmentsConformMesh,
            Gedim::ConformerMeshPolygon::ConformerMeshPolygonConfiguration::Types::OnlyOnEdges,
            options);

        // check meshes 1D
        for (unsigned int i = 0; i < numSegments; i++)
        {
            for (const auto &point : segmentsConformMesh[i].Points)
                Gedim::Output::Assert(point.second.Vertex2DIds.size() == 1);

            for (const auto &segment : segmentsConformMesh[i].Segments)
                Gedim::Output::Assert(segment.Edge2DIds.size() == 1);
        }

        // check meshes 2D
        Gedim::MeshUtilities::CheckMesh2DConfiguration config;
        meshUtilities.CheckMesh2D(config, geometryUtilities, meshDAO);

        conformMeshUtilities.AddConformedMeshProperties(segmentsConformMesh, meshDAO);

        meshUtilities.ExportMeshToVTU(meshDAO, exportFolder, "ConformedMesh");

        // export mesh
        //      Gedim::MeshFromCsvUtilities meshFromCsvUtilities;
        //      Gedim::MeshFromCsvUtilities::Configuration exportConfiguration;
        //      exportConfiguration.Folder = exportFolder;
        //      Gedim::MeshDAOExporterToCsv exporter(meshFromCsvUtilities);
        //      exporter.Export(exportConfiguration,
        //                      meshDAO);
    }
    catch (const exception &exception)
    {
        cerr << exception.what() << endl;
        FAIL();
    }
}
} // namespace GedimUnitTesting

#endif // __TEST_ConformMeshUtilities_H
