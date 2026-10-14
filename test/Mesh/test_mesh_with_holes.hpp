// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __test_mesh_with_holes_H
#define __test_mesh_with_holes_H

#include "CreateMeshWithHoles.hpp"
#include "GraphUtilities.hpp"
#include "ImportExportUtilities.hpp"
#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VoroInterface.hpp"
#include "VtkMeshExporter.hpp"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace GedimUnitTesting
{
TEST(Test_mesh_with_holes, Test_no_mark)
{
    GTEST_SKIP(); // NOT WORKING YET

    std::string exportFolder = "./Export/Test_mesh_with_holes/Test_no_mark";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::GraphUtilities graphUtilities;

    // Create domain 3D
    const double domain_3D_length = 1.0;
    const double domain_3D_height = 1.0;
    const double domain_3D_width = 1.0;
    const unsigned int domain_3D_length_num_cells = 1;
    const unsigned int domain_3D_height_num_cells = 2;
    const unsigned int domain_3D_width_num_cells = 2;
    const Eigen::Vector3d domain_3D_origin(0.0, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_length_direction(domain_3D_length, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_height_direction(0.0, 0.0, domain_3D_height);
    const Eigen::Vector3d domain_3D_width_direction(0.0, domain_3D_width, 0.0);

    const auto domain_3D = geometryUtilities.CreateParallelepipedWithOrigin(domain_3D_origin,
                                                                            domain_3D_length_direction,
                                                                            domain_3D_height_direction,
                                                                            domain_3D_width_direction);

    {
        Gedim::VTKUtilities exporter;
        exporter.AddPolyhedron(domain_3D.Vertices, domain_3D.Edges, domain_3D.Faces);
        exporter.Export(exportFolder + "/Domain_3D.vtu");
    }

    // Create mesh 3D
    const std::vector<double> domain_3D_mesh_length_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_length_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_height_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_height_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_width_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_width_num_cells + 1, 0.0, 1.0, 1);

    Gedim::MeshMatrices mesh_3D_data;
    Gedim::MeshMatricesDAO mesh_3D(mesh_3D_data);

    unsigned int starting_mesh_type = 0;

    switch (starting_mesh_type)
    {
    case 0: {
        const Eigen::Vector3d code_domain_3D_origin = domain_3D.Vertices.col(0);
        const Eigen::Vector3d code_domain_3D_length_direction = domain_3D.Vertices.col(1) - domain_3D_origin;
        const Eigen::Vector3d code_domain_3D_height_direction = domain_3D.Vertices.col(3) - domain_3D_origin;
        const Eigen::Vector3d code_domain_3D_width_direction = domain_3D.Vertices.col(4) - domain_3D_origin;

        meshUtilities.CreateParallelepipedMesh(code_domain_3D_origin,
                                               code_domain_3D_length_direction,
                                               code_domain_3D_height_direction,
                                               code_domain_3D_width_direction,
                                               domain_3D_mesh_length_coordinates,
                                               domain_3D_mesh_height_coordinates,
                                               domain_3D_mesh_width_coordinates,
                                               mesh_3D);
    }
    break;
    case 1: {
        Gedim::External::MeshUtilities::CreatePolyhedralMesh(geometryUtilities,
                                                             domain_3D.Vertices,
                                                             domain_3D.Edges,
                                                             domain_3D.Faces,
                                                             domain_3D_length_num_cells * domain_3D_width_num_cells * domain_3D_height_num_cells,
                                                             10,
                                                             mesh_3D,
                                                             256);
    }
    break;
    default:
        throw std::runtime_error("Unknown mesh type");
    }

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
    }

    auto empty_function = [&](const Eigen::MatrixXd &points) {
        Eigen::VectorXi marks(points.cols());

        auto surface_one = [domain_3D_width, domain_3D_height](const Eigen::Vector3d &p) {
            return Eigen::Vector3d(p.x(), p.y(), domain_3D_height * 0.251);
        };

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            const Eigen::Vector3d point = points.col(p);

            if (geometryUtilities.IsValueGreaterOrEqual(point.z(), surface_one(point).z(), geometryUtilities.Tolerance1D()))
            {
                marks[p] = 1;
                continue;
            }

            marks[p] = 0;
        }

        return marks;
    };

    //    auto empty_function = [&](const Eigen::MatrixXd& points)
    //    {
    //      return Eigen::VectorXi::Zero(points.cols());
    //    };

    const auto agglomeration_info =
        Gedim::External::MeshUtilities::AgglomerateMeshByFunction(geometryUtilities, graphUtilities, mesh_3D_data, empty_function, 0.0);

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_agglomerated";

        {
            Gedim::Output::CreateFolder(mesh_3D_export_folder);
            Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
        }

        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(agglomeration_info.Cell3Ds_mark.begin(), agglomeration_info.Cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                mesh_3D.Cell0DsCoordinates(),
                mesh_3D.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark.vtu");
        }

        {
            std::unordered_set<unsigned int> original_cell3Ds_index;
            Gedim::VTKUtilities exporter;

            for (unsigned int c_3D = 0; c_3D < mesh_3D.Cell3DTotalNumber(); c_3D++)
            {
                if (!mesh_3D.Cell3DIsActive(c_3D))
                    continue;

                const auto &cell3D_tetra_vertices = agglomeration_info.Cell3Ds_tetra_vertices.at(c_3D);

                for (const auto &tetra_vertices : cell3D_tetra_vertices)
                {
                    std::vector<double> cell3D(1, c_3D);

                    const auto tetra_polyhedron = geometryUtilities.CreateTetrahedronWithVertices(tetra_vertices.col(0),
                                                                                                  tetra_vertices.col(1),
                                                                                                  tetra_vertices.col(2),
                                                                                                  tetra_vertices.col(3));
                    exporter.AddPolyhedron(
                        tetra_polyhedron.Vertices,
                        tetra_polyhedron.Edges,
                        tetra_polyhedron.Faces,
                        {{"cell3D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell3D.size()), cell3D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell3Ds_tetra.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c_2D = 0; c_2D < mesh_3D.Cell2DTotalNumber(); c_2D++)
            {
                if (!mesh_3D.Cell2DIsActive(c_2D))
                    continue;

                const auto &cell2D_triangle_vertices = agglomeration_info.Cell2Ds_triangles_vertices.at(c_2D);

                for (const auto &triangle_vertices : cell2D_triangle_vertices)
                {
                    std::vector<double> cell2D(1, c_2D);

                    exporter.AddPolygon(
                        triangle_vertices,
                        {{"cell2D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell2D.size()), cell2D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell2Ds_triangles.vtu");
        }
    }

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_Data";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        meshUtilities.ExportMeshToCsv(mesh_3D, ';', mesh_3D_export_folder);

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell3Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell3Ds_mark << std::endl;
            file << std::scientific << agglomeration_info.Cell3Ds_tetra_vertices << std::endl;

            file.close();
        }

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell2Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell2Ds_triangles_vertices << std::endl;

            file.close();
        }
    }

    Gedim::MeshUtilities::CheckMesh3DConfiguration checkMeshDataConfig;
    checkMeshDataConfig.Cell2D_CheckConvexity = false;
    checkMeshDataConfig.Cell3D_CheckConvexity = false;
    checkMeshDataConfig.Cell2D_CheckMeasure = false;
    checkMeshDataConfig.Cell3D_CheckMeasure = false;
    ASSERT_NO_THROW(meshUtilities.CheckMesh3D(checkMeshDataConfig, geometryUtilities, mesh_3D));

    const auto mesh_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                           mesh_3D,
                                                                           agglomeration_info.Cell3Ds_tetra_vertices,
                                                                           agglomeration_info.Cell2Ds_triangles_vertices);

    {
        const std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_agglomerated_cell3Ds";

        for (unsigned int c_3D = 0; c_3D < mesh_3D.Cell3DTotalNumber(); c_3D++)
        {
            if (!mesh_3D.Cell3DIsActive(c_3D))
                continue;

            const std::string mesh_3D_export_folder_c_3D = mesh_3D_export_folder + "/c_" + std::to_string(c_3D);
            Gedim::Output::CreateFolder(mesh_3D_export_folder_c_3D);

            Gedim::External::MeshUtilities::ExportCell3DToVTU(geometryUtilities,
                                                              mesh_3D,
                                                              c_3D,
                                                              mesh_geometric_data.Cell3DsVertices.at(c_3D),
                                                              mesh_geometric_data.Cell3DsTetrahedronPoints.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFaces3DTriangulations.at(c_3D),
                                                              mesh_geometric_data.Cell3DsVolumes.at(c_3D),
                                                              mesh_geometric_data.Cell3DsCentroids.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesTranslations.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesRotationMatrices.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesAreas.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFaces2DVertices.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFaces3DVertices.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesEdgeLengths.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesEdgeDirections.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesEdge2DTangents.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesEdge2DNormals.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesNormals.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFacesNormalDirections.at(c_3D),
                                                              mesh_geometric_data.Cell3DsFaces2DCentroids.at(c_3D),
                                                              mesh_3D_export_folder_c_3D);
        }
    }

    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig, geometryUtilities, mesh_3D, mesh_geometric_data));

    const bool is_agglomerated_convex = false;

    if (!is_agglomerated_convex)
        return;

    const auto mesh_geometric_data_convex = meshUtilities.FillMesh3DGeometricData(geometryUtilities, mesh_3D);

    {
        using namespace Gedim;

        for (unsigned int c = 0; c < mesh_3D.Cell3DTotalNumber(); c++)
        {
            for (unsigned int f = 0; f < mesh_3D.Cell3DNumberFaces(c); f++)
            {
                ASSERT_TRUE(geometryUtilities.AreValuesEqual(mesh_geometric_data_convex.Cell3DsFacesAreas.at(c).at(f),
                                                             mesh_geometric_data.Cell3DsFacesAreas.at(c).at(f),
                                                             geometryUtilities.Tolerance2D()));

                ASSERT_TRUE(
                    geometryUtilities.PointsAreCoincident(mesh_geometric_data_convex.Cell3DsFaces2DCentroids.at(c).at(f),
                                                          mesh_geometric_data.Cell3DsFaces2DCentroids.at(c).at(f)));

                const auto &new_triangles_3D = mesh_geometric_data.Cell3DsFaces3DTriangulations.at(c).at(f);
                const auto &original_triangles_3D = mesh_geometric_data_convex.Cell3DsFaces3DTriangulations.at(c).at(f);

                const auto &new_triangles_2D = mesh_geometric_data.Cell3DsFaces2DTriangulations.at(c).at(f);
                const auto &original_triangles_2D = mesh_geometric_data_convex.Cell3DsFaces2DTriangulations.at(c).at(f);

                ASSERT_EQ(original_triangles_3D.size(), new_triangles_3D.size());
                ASSERT_EQ(original_triangles_2D.size(), new_triangles_2D.size());

                for (unsigned int t = 0; t < original_triangles_3D.size(); t++)
                {
                    for (unsigned int v = 0; v < 3; v++)
                    {
                        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(original_triangles_3D.at(t).col(v),
                                                                          new_triangles_3D.at(t).col(v)));
                    }
                }

                for (unsigned int t = 0; t < original_triangles_2D.size(); t++)
                {
                    for (unsigned int v = 0; v < 3; v++)
                    {
                        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(original_triangles_2D.at(t).col(v),
                                                                          new_triangles_2D.at(t).col(v)));
                    }
                }
            }
        }
    }

    ASSERT_EQ(mesh_geometric_data.Cell3DsBoundingBox, mesh_geometric_data_convex.Cell3DsBoundingBox);

    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsVertices, mesh_geometric_data.Cell3DsVertices);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsEdges, mesh_geometric_data.Cell3DsEdges);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFaces, mesh_geometric_data.Cell3DsFaces);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsDiameters, mesh_geometric_data.Cell3DsDiameters);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsEdgesCentroid, mesh_geometric_data.Cell3DsEdgesCentroid);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsEdgeLengths, mesh_geometric_data.Cell3DsEdgeLengths);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsEdgeTangents, mesh_geometric_data.Cell3DsEdgeTangents);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsEdgeDirections, mesh_geometric_data.Cell3DsEdgeDirections);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsTetrahedronPoints, mesh_geometric_data.Cell3DsTetrahedronPoints);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesTranslations, mesh_geometric_data.Cell3DsFacesTranslations);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesRotationMatrices, mesh_geometric_data.Cell3DsFacesRotationMatrices);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesNormals, mesh_geometric_data.Cell3DsFacesNormals);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesNormalDirections, mesh_geometric_data.Cell3DsFacesNormalDirections);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdgeDirections, mesh_geometric_data.Cell3DsFacesEdgeDirections);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFaces3DVertices, mesh_geometric_data.Cell3DsFaces3DVertices);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFaces2DVertices, mesh_geometric_data.Cell3DsFaces2DVertices);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesDiameters, mesh_geometric_data.Cell3DsFacesDiameters);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdgeLengths, mesh_geometric_data.Cell3DsFacesEdgeLengths);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdge3DTangents, mesh_geometric_data.Cell3DsFacesEdge3DTangents);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdges3DCentroid, mesh_geometric_data.Cell3DsFacesEdges3DCentroid);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdge2DTangents, mesh_geometric_data.Cell3DsFacesEdge2DTangents);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdges2DCentroid, mesh_geometric_data.Cell3DsFacesEdges2DCentroid);
    ASSERT_EQ(mesh_geometric_data_convex.Cell3DsFacesEdge2DNormals, mesh_geometric_data.Cell3DsFacesEdge2DNormals);

    for (unsigned int c = 0; c < mesh_3D.Cell3DTotalNumber(); c++)
    {
        ASSERT_TRUE(geometryUtilities.AreValuesEqual(mesh_geometric_data_convex.Cell3DsVolumes.at(c),
                                                     mesh_geometric_data.Cell3DsVolumes.at(c),
                                                     geometryUtilities.Tolerance3D()));
        ASSERT_TRUE(geometryUtilities.PointsAreCoincident(mesh_geometric_data_convex.Cell3DsCentroids.at(c),
                                                          mesh_geometric_data.Cell3DsCentroids.at(c)));
    }
}

TEST(Test_mesh_with_holes, Test_one_hole)
{
    GTEST_SKIP(); // NOT WORKING YET

    std::string exportFolder = "./Export/Test_mesh_with_holes/Test_one_hole";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::GraphUtilities graphUtilities;

    // Create domain 3D
    const double domain_3D_length = 1.0;
    const double domain_3D_height = 1.0;
    const double domain_3D_width = 1.0;
    const unsigned int domain_3D_length_num_cells = 5;
    const unsigned int domain_3D_height_num_cells = 5;
    const unsigned int domain_3D_width_num_cells = 5;
    const Eigen::Vector3d domain_3D_origin(0.0, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_length_direction(domain_3D_length, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_height_direction(0.0, 0.0, domain_3D_height);
    const Eigen::Vector3d domain_3D_width_direction(0.0, domain_3D_width, 0.0);

    {
        const auto domain_3D = geometryUtilities.CreateParallelepipedWithOrigin(domain_3D_origin,
                                                                                domain_3D_length_direction,
                                                                                domain_3D_height_direction,
                                                                                domain_3D_width_direction);

        Gedim::VTKUtilities exporter;
        exporter.AddPolyhedron(domain_3D.Vertices, domain_3D.Edges, domain_3D.Faces);
        exporter.Export(exportFolder + "/Domain_3D.vtu");
    }

    // Create mesh 3D
    const std::vector<double> domain_3D_mesh_length_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_length_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_height_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_height_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_width_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_width_num_cells + 1, 0.0, 1.0, 1);

    Gedim::MeshMatrices mesh_3D_data;
    Gedim::MeshMatricesDAO mesh_3D(mesh_3D_data);
    meshUtilities.CreateParallelepipedMesh(domain_3D_origin,
                                           domain_3D_length_direction,
                                           domain_3D_height_direction,
                                           domain_3D_width_direction,
                                           domain_3D_mesh_length_coordinates,
                                           domain_3D_mesh_height_coordinates,
                                           domain_3D_mesh_width_coordinates,
                                           mesh_3D);

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
    }

    auto rock_function = [&](const Eigen::MatrixXd &points) {
        Eigen::VectorXi marks(points.cols());

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            const Eigen::Vector3d point = points.col(p);

            if (!geometryUtilities.IsValueGreaterOrEqual(
                    (Eigen::Vector3d(0.75 * domain_3D_length, 0.75 * domain_3D_width, 0.75 * domain_3D_height) - point).norm(),
                    domain_3D_length * 0.25,
                    geometryUtilities.Tolerance1D()))
            {
                marks[p] = 1;
                continue;
            }

            if (!geometryUtilities.IsValueGreaterOrEqual(
                    (Eigen::Vector3d(0.25 * domain_3D_length, 0.25 * domain_3D_width, 0.25 * domain_3D_height) - point).norm(),
                    domain_3D_length * 0.25,
                    geometryUtilities.Tolerance1D()))
            {
                marks[p] = 1;
                continue;
            }

            marks[p] = 0;
        }

        return marks;
    };

    const auto agglomeration_info =
        Gedim::External::MeshUtilities::AgglomerateMeshByFunction(geometryUtilities, graphUtilities, mesh_3D_data, rock_function, 0.0);

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_agglomerated";

        {
            Gedim::Output::CreateFolder(mesh_3D_export_folder);
            Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
        }

        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(agglomeration_info.Cell3Ds_mark.begin(), agglomeration_info.Cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                mesh_3D.Cell0DsCoordinates(),
                mesh_3D.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark.vtu");
        }

        {
            std::unordered_set<unsigned int> original_cell3Ds_index;
            Gedim::VTKUtilities exporter;

            for (unsigned int c_3D = 0; c_3D < mesh_3D.Cell3DTotalNumber(); c_3D++)
            {
                if (!mesh_3D.Cell3DIsActive(c_3D))
                    continue;

                const auto &cell3D_tetra_vertices = agglomeration_info.Cell3Ds_tetra_vertices.at(c_3D);

                for (const auto &tetra_vertices : cell3D_tetra_vertices)
                {
                    std::vector<double> cell3D(1, c_3D);

                    const auto tetra_polyhedron = geometryUtilities.CreateTetrahedronWithVertices(tetra_vertices.col(0),
                                                                                                  tetra_vertices.col(1),
                                                                                                  tetra_vertices.col(2),
                                                                                                  tetra_vertices.col(3));
                    exporter.AddPolyhedron(
                        tetra_polyhedron.Vertices,
                        tetra_polyhedron.Edges,
                        tetra_polyhedron.Faces,
                        {{"cell3D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell3D.size()), cell3D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell3Ds_tetra.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c_2D = 0; c_2D < mesh_3D.Cell2DTotalNumber(); c_2D++)
            {
                if (!mesh_3D.Cell2DIsActive(c_2D))
                    continue;

                const auto &cell2D_triangle_vertices = agglomeration_info.Cell2Ds_triangles_vertices.at(c_2D);

                for (const auto &triangle_vertices : cell2D_triangle_vertices)
                {
                    std::vector<double> cell2D(1, c_2D);

                    exporter.AddPolygon(
                        triangle_vertices,
                        {{"cell2D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell2D.size()), cell2D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell2Ds_triangles.vtu");
        }
    }

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_Data";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        meshUtilities.ExportMeshToCsv(mesh_3D, ';', mesh_3D_export_folder);

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell3Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell3Ds_mark << std::endl;
            file << std::scientific << agglomeration_info.Cell3Ds_tetra_vertices << std::endl;

            file.close();
        }

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell2Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell2Ds_triangles_vertices << std::endl;

            file.close();
        }
    }

    Gedim::MeshUtilities::CheckMesh3DConfiguration checkMeshDataConfig;
    checkMeshDataConfig.Cell2D_CheckConvexity = false;
    checkMeshDataConfig.Cell3D_CheckConvexity = false;
    checkMeshDataConfig.Cell2D_CheckMeasure = false;
    checkMeshDataConfig.Cell3D_CheckMeasure = false;
    ASSERT_NO_THROW(meshUtilities.CheckMesh3D(checkMeshDataConfig, geometryUtilities, mesh_3D));

    const auto mesh_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                           mesh_3D,
                                                                           agglomeration_info.Cell3Ds_tetra_vertices,
                                                                           agglomeration_info.Cell2Ds_triangles_vertices);
    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig, geometryUtilities, mesh_3D, mesh_geometric_data));
}

TEST(Test_mesh_with_holes, Test_one_plane)
{
    GTEST_SKIP(); // NOT WORKING YET

    std::string exportFolder = "./Export/Test_mesh_with_holes/Test_one_plane";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.MinTolerance = 1.0e-14;
    geometryUtilitiesConfig.Tolerance1D = 1.0e-6;
    geometryUtilitiesConfig.Tolerance2D = 1.0e-12;
    geometryUtilitiesConfig.Tolerance3D = 1.0e-10;
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;
    Gedim::GraphUtilities graphUtilities;

    // Create domain 3D
    const double domain_3D_length = 1.0;
    const double domain_3D_height = 2.0;
    const double domain_3D_width = 3.0;
    const unsigned int domain_3D_length_num_cells = 5;
    const unsigned int domain_3D_height_num_cells = 5;
    const unsigned int domain_3D_width_num_cells = 5;
    const Eigen::Vector3d domain_3D_origin(0.0, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_length_direction(domain_3D_length, 0.0, 0.0);
    const Eigen::Vector3d domain_3D_height_direction(0.0, 0.0, domain_3D_height);
    const Eigen::Vector3d domain_3D_width_direction(0.0, domain_3D_width, 0.0);

    {
        const auto domain_3D = geometryUtilities.CreateParallelepipedWithOrigin(domain_3D_origin,
                                                                                domain_3D_length_direction,
                                                                                domain_3D_height_direction,
                                                                                domain_3D_width_direction);

        Gedim::VTKUtilities exporter;
        exporter.AddPolyhedron(domain_3D.Vertices, domain_3D.Edges, domain_3D.Faces);
        exporter.Export(exportFolder + "/Domain_3D.vtu");
    }

    // Create mesh 3D
    const std::vector<double> domain_3D_mesh_length_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_length_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_height_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_height_num_cells + 1, 0.0, 1.0, 1);
    const std::vector<double> domain_3D_mesh_width_coordinates =
        geometryUtilities.EquispaceCoordinates(domain_3D_width_num_cells + 1, 0.0, 1.0, 1);

    Gedim::MeshMatrices mesh_3D_data;
    Gedim::MeshMatricesDAO mesh_3D(mesh_3D_data);
    meshUtilities.CreateParallelepipedMesh(domain_3D_origin,
                                           domain_3D_length_direction,
                                           domain_3D_height_direction,
                                           domain_3D_width_direction,
                                           domain_3D_mesh_length_coordinates,
                                           domain_3D_mesh_height_coordinates,
                                           domain_3D_mesh_width_coordinates,
                                           mesh_3D);

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
    }

    auto surface_function = [&geometryUtilities, domain_3D_width, domain_3D_height](const Eigen::MatrixXd &points) {
        Eigen::VectorXi marks(points.cols());

        auto surface_one = [domain_3D_width, domain_3D_height](const Eigen::Vector3d &p) {
            return Eigen::Vector3d(p.x(), p.y(), domain_3D_height * (0.1 * (1.0 - std::cos(M_PI * p.y() / domain_3D_width)) + 0.45));
        };

        auto surface_two = [domain_3D_width, domain_3D_height](const Eigen::Vector3d &p) {
            return Eigen::Vector3d(p.x(), p.y(), domain_3D_height * (0.1 * (1.0 - std::cos(M_PI * p.y() / domain_3D_width)) + 0.15));
        };

        for (unsigned int p = 0; p < points.cols(); p++)
        {
            const Eigen::Vector3d point = points.col(p);

            if (geometryUtilities.IsValueGreaterOrEqual(point.z(), surface_one(point).z(), geometryUtilities.Tolerance1D()))
            {
                marks[p] = 1;
                continue;
            }

            if (!geometryUtilities.IsValueGreaterOrEqual(point.z(), surface_two(point).z(), geometryUtilities.Tolerance1D()))
            {
                marks[p] = 1;
                continue;
            }

            marks[p] = 0;
        }

        return marks;
    };

    const auto agglomeration_info =
        Gedim::External::MeshUtilities::AgglomerateMeshByFunction(geometryUtilities, graphUtilities, mesh_3D_data, surface_function, 0.0);

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_agglomerated";

        {
            Gedim::Output::CreateFolder(mesh_3D_export_folder);
            Gedim::External::MeshUtilities::ExportMeshToVTU(mesh_3D, mesh_3D_export_folder, "Mesh3D");
        }

        {
            Gedim::VTKUtilities exporter;
            std::vector<double> marks(agglomeration_info.Cell3Ds_mark.begin(), agglomeration_info.Cell3Ds_mark.end());

            exporter.AddPolyhedrons(
                mesh_3D.Cell0DsCoordinates(),
                mesh_3D.Cell3DsFacesVertices(),
                {{"Mark", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(marks.size()), marks.data()}});
            exporter.Export(mesh_3D_export_folder + "/cell3Ds_mark.vtu");
        }

        {
            std::unordered_set<unsigned int> original_cell3Ds_index;
            Gedim::VTKUtilities exporter;

            for (unsigned int c_3D = 0; c_3D < mesh_3D.Cell3DTotalNumber(); c_3D++)
            {
                if (!mesh_3D.Cell3DIsActive(c_3D))
                    continue;

                const auto &cell3D_tetra_vertices = agglomeration_info.Cell3Ds_tetra_vertices.at(c_3D);

                for (const auto &tetra_vertices : cell3D_tetra_vertices)
                {
                    std::vector<double> cell3D(1, c_3D);

                    const auto tetra_polyhedron = geometryUtilities.CreateTetrahedronWithVertices(tetra_vertices.col(0),
                                                                                                  tetra_vertices.col(1),
                                                                                                  tetra_vertices.col(2),
                                                                                                  tetra_vertices.col(3));
                    exporter.AddPolyhedron(
                        tetra_polyhedron.Vertices,
                        tetra_polyhedron.Edges,
                        tetra_polyhedron.Faces,
                        {{"cell3D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell3D.size()), cell3D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell3Ds_tetra.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c_2D = 0; c_2D < mesh_3D.Cell2DTotalNumber(); c_2D++)
            {
                if (!mesh_3D.Cell2DIsActive(c_2D))
                    continue;

                const auto &cell2D_triangle_vertices = agglomeration_info.Cell2Ds_triangles_vertices.at(c_2D);

                for (const auto &triangle_vertices : cell2D_triangle_vertices)
                {
                    std::vector<double> cell2D(1, c_2D);

                    exporter.AddPolygon(
                        triangle_vertices,
                        {{"cell2D", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(cell2D.size()), cell2D.data()}});
                }
            }

            exporter.Export(mesh_3D_export_folder + "/cell2Ds_triangles.vtu");
        }
    }

    {
        std::string mesh_3D_export_folder = exportFolder + "/Mesh3D_Data";
        Gedim::Output::CreateFolder(mesh_3D_export_folder);
        meshUtilities.ExportMeshToCsv(mesh_3D, ';', mesh_3D_export_folder);

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell3Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell3Ds_mark << std::endl;
            file << std::scientific << agglomeration_info.Cell3Ds_tetra_vertices << std::endl;

            file.close();
        }

        {
            using namespace Gedim_ImportExport_Utilities;

            const std::string file_path = mesh_3D_export_folder + "/Cell2Ds_Additional_Info.txt";

            std::ofstream file(file_path);

            if (!file.is_open())
                throw std::runtime_error("Export file '" + file_path + "' not opened");

            file.precision(16);

            file << std::scientific << agglomeration_info.Cell2Ds_triangles_vertices << std::endl;

            file.close();
        }
    }

    Gedim::MeshUtilities::CheckMesh3DConfiguration checkMeshDataConfig;
    checkMeshDataConfig.Cell2D_CheckConvexity = false;
    checkMeshDataConfig.Cell3D_CheckConvexity = false;
    checkMeshDataConfig.Cell2D_CheckMeasure = false;
    checkMeshDataConfig.Cell3D_CheckMeasure = false;
    ASSERT_NO_THROW(meshUtilities.CheckMesh3D(checkMeshDataConfig, geometryUtilities, mesh_3D));

    const auto mesh_geometric_data = meshUtilities.FillMesh3DGeometricData(geometryUtilities,
                                                                           mesh_3D,
                                                                           agglomeration_info.Cell3Ds_tetra_vertices,
                                                                           agglomeration_info.Cell2Ds_triangles_vertices);
    Gedim::MeshUtilities::CheckMeshGeometricData3DConfiguration checkMeshGeometricDataConfig;
    ASSERT_NO_THROW(meshUtilities.CheckMeshGeometricData3D(checkMeshGeometricDataConfig, geometryUtilities, mesh_3D, mesh_geometric_data));
}
} // namespace GedimUnitTesting

#endif // __test_3D1D_hexa_mesh_H
