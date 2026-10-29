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

#ifndef __TEST_meshUtilities3D_MarchingSquares_H
#define __TEST_meshUtilities3D_MarchingSquares_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace GedimUnitTesting
{

TEST(TestMeshUtilities, TestMeshUtilities_2D_MarchingSquares)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_2D_MarchingSquares";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
        const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d square_LengthTangent = square.col(1) - square.col(0);
        const Eigen::Vector3d square_HeightTangent = square.col(3) - square.col(0);

        const std::vector<double> lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4, 0.0, 1.0, 1);
        const std::vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(5, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateRectangleMesh(square.col(0),
                                                  square_LengthTangent,
                                                  square_HeightTangent,
                                                  lengthMeshCurvilinearCoordinates,
                                                  heightMeshCurvilinearCoordinates,
                                                  mesh);
    }

    {
        const std::string mesh_folder = exportFolder + "/original_mesh";
        Gedim::Output::CreateFolder(mesh_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    auto plane_x = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.4, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.0, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
    };

    auto generic_plane = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.6, 0.4, 0.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.5, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
    };

    auto circle = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = 1.0;
        const double radius = 0.45;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5);
        return x * x + y * y - radius * radius;
    };

    auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y();
        return std::sin(x) * std::cos(y) + std::sin(y) - c * scale;
    };

    auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
        return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
    };

    auto surface_implicit_equation = [&plane_x, &generic_plane, &circle, &gyroid, &gyroid_shifted](const Eigen::MatrixXd &points) {
        Eigen::VectorXd surface_values(points.cols());
        for (unsigned int p = 0; p < surface_values.size(); p++)
            surface_values[p] = circle(points.col(p), 1.0, 0.0);

        return surface_values;
    };

    const auto marching_squares_result =
        Gedim::MeshUtilities::MarchingSquares(geometry_utilities, surface_implicit_equation, mesh);

    {
        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
            {
                const auto &cell0D_surface_point = marching_squares_result.surface_intersections.cell0Ds_surface_point.at(c);

                double cell0D_sign = 0.0;
                switch (cell0D_surface_point.surface_position)
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    cell0D_sign = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    cell0D_sign = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    cell0D_sign = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPoint(
                    mesh.Cell0DCoordinates(c),
                    {{"sign", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_sign},
                     {"value", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_surface_point.surface_value}});
            }

            exporter.Export(exportFolder + "/cell0Ds_surface_position.vtu");
        }

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int p = 0; p < marching_squares_result.surface_intersections.surface_intersection_points.size(); ++p)
            {
                const auto &surface_point = marching_squares_result.surface_intersections.surface_intersection_points.at(p);
                const double surface_id = surface_point.mesh_id;

                double surface_type = 0.0;
                switch (surface_point.Type)
                {
                case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Vertex:
                    surface_type = 1.0;
                    break;
                case Gedim::MeshUtilities::IntersectMeshWithSurfaceResult::SurfacePointData::Types::Edge:
                    surface_type = 2.0;
                    break;
                default:
                    throw std::runtime_error("Surface Type not managed");
                }

                exporter.AddPoint(surface_point.coordinate,
                                  {{"type", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &surface_type},
                                   {"id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &surface_id}});
            }

            exporter.Export(exportFolder + "/marching_square_points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < marching_squares_result.surface_segments.surface_segments.size(); ++p)
            {
                const auto &segment_ids = marching_squares_result.surface_segments.surface_segments.at(p);
                Eigen::MatrixXd segment;
                segment.resize(3, 2);
                segment.col(0) = marching_squares_result.surface_intersections.surface_intersection_points
                                     .at(segment_ids.indices.at(0))
                                     .coordinate;
                segment.col(1) = marching_squares_result.surface_intersections.surface_intersection_points
                                     .at(segment_ids.indices.at(1))
                                     .coordinate;
                const double cell2D_id = segment_ids.cell2D_id;

                exporter.AddSegment(segment,
                                    {{"cell2D_id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell2D_id}});
            }

            exporter.Export(exportFolder + "/marching_squares_segments.vtu");
        }
    }

    const auto surface_intersection_graph =
        Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_squares_result.surface_intersections,
                                                                  marching_squares_result.surface_segments,
                                                                  {},
                                                                  mesh);

    {
        const std::string graph_folder = exportFolder + "/intersection_graph";
        Gedim::Output::CreateFolder(graph_folder);

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c = 0; c < surface_intersection_graph.vertices.size(); ++c)
            {
                const auto &vertex = surface_intersection_graph.vertices.at(c);

                const double vertex_original_surface = vertex.original_surface_point_index;

                double vertex_position = 0.0;
                switch (vertex.surface_position)
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    vertex_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    vertex_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    vertex_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPoint(
                    vertex.coordinate,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_position},
                     {"original_surface_index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_original_surface}});
            }

            exporter.Export(graph_folder + "/vertices.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < surface_intersection_graph.edges.size(); ++p)
            {
                const auto &edge = surface_intersection_graph.edges.at(p);

                Eigen::MatrixXd segment;
                segment.resize(3, 2);
                segment.col(0) = surface_intersection_graph.vertices.at(edge.extremes.at(0)).coordinate;
                segment.col(1) = surface_intersection_graph.vertices.at(edge.extremes.at(1)).coordinate;

                double edge_position = 0.0;
                switch (edge.surface_position)
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    edge_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    edge_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    edge_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddSegment(
                    segment,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &edge_position}});
            }

            exporter.Export(graph_folder + "/edges.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int f = 0; f < surface_intersection_graph.faces.size(); ++f)
            {
                const auto &face = surface_intersection_graph.faces.at(f);

                Eigen::MatrixXd face_vertices;
                face_vertices.resize(3, face.intersection_graph_vertices.size());
                for (unsigned int fv = 0; fv < face.intersection_graph_vertices.size(); ++fv)
                {
                    face_vertices.col(fv) =
                        surface_intersection_graph.vertices.at(face.intersection_graph_vertices.at(fv)).coordinate;
                }

                double face_position = 0.0;
                switch (face.surface_position)
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    face_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    face_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    face_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPolygon(
                    face_vertices,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &face_position}});
            }

            exporter.Export(graph_folder + "/faces.vtu");
        }
    }

    Gedim::MeshMatrices conformed_mesh_data;
    Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

    const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(mesh,
                                                                           surface_intersection_graph,
                                                                           conformed_mesh);

    {
        const std::string mesh_folder = exportFolder + "/conformed_mesh";
        Gedim::Output::CreateFolder(mesh_folder);

        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c = 0; c < conformed_mesh.Cell0DTotalNumber(); ++c)
            {
                const double vertex_original_cell0D = conform_result.new_cell0Ds_to_original_cell0D.at(c);
                const double vertex_original_surface = conform_result.new_cell0Ds_to_surface_vertex.at(c);

                double vertex_position = 0.0;
                switch (conform_result.new_cell0Ds_surface_position.at(c))
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    vertex_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    vertex_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    vertex_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPoint(
                    conformed_mesh.Cell0DCoordinates(c),
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_position},
                     {"original_surface_vertex", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_original_surface},
                     {"original_cell0D_index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &vertex_original_cell0D}});
            }

            exporter.Export(mesh_folder + "/conformed_mesh_cell0Ds.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < conformed_mesh.Cell1DTotalNumber(); ++p)
            {
                const double edge_original_cell1D = conform_result.new_cell1Ds_to_original_cell1D.at(p);
                const double edge_original_surface = conform_result.new_cell1Ds_to_surface_segment.at(p);

                Eigen::MatrixXd segment;
                segment.resize(3, 2);
                segment.col(0) = conformed_mesh.Cell1DOriginCoordinates(p);
                segment.col(1) = conformed_mesh.Cell1DEndCoordinates(p);

                double edge_position = 0.0;
                switch (conform_result.new_cell1Ds_surface_position.at(p))
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    edge_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    edge_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    edge_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddSegment(
                    segment,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &edge_position},
                     {"original_surface_edge", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &edge_original_surface},
                     {"original_cell1D_index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &edge_original_cell1D}});
            }

            exporter.Export(mesh_folder + "/conformed_mesh_cell1Ds.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int f = 0; f < conformed_mesh.Cell2DTotalNumber(); ++f)
            {
                const double face_original_cell2D = conform_result.new_cell2Ds_to_original_cell2D.at(f);
                const double face_original_surface = conform_result.new_cell2Ds_to_surface_triangle.at(f);
                const auto face_vertices = conformed_mesh.Cell2DVerticesCoordinates(f);

                double face_position = 0.0;
                switch (conform_result.new_cell2Ds_surface_position.at(f))
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    face_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    face_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    face_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPolygon(
                    face_vertices,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &face_position},
                     {"original_surface_triangle", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &face_original_surface},
                     {"original_cell2D_index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &face_original_cell2D}});
            }

            exporter.Export(mesh_folder + "/conformed_mesh_cell2Ds.vtu");
        }
    }

    Gedim::MeshUtilities::ComputeCell1DCell2DNeighbours(conformed_mesh);
    Gedim::MeshUtilities::CheckMesh2DConfiguration check_mesh_config;
    Gedim::MeshUtilities::CheckMesh2D(check_mesh_config, geometry_utilities, conformed_mesh);
}

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingSquares_H
