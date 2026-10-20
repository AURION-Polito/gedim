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

    Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, exportFolder, "original_mesh");

    auto circle = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = 1.0;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5);
        return x * x + y * y - std::sqrt(0.05);
    };

    auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y();
        return std::sin(x) * std::cos(y) + std::sin(y) - c;
    };

    auto surface_implicit_equation = [&circle, &gyroid](const Eigen::MatrixXd &points) {
        float L = 1.0f;
        float c = 0.0f;

        Eigen::VectorXd surface_values(points.cols());
        for (unsigned int p = 0; p < surface_values.size(); p++)
            surface_values[p] = circle(points.col(p), L, c);

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
                case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
                    cell0D_sign = 0.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Positive:
                    cell0D_sign = +1.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Negative:
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

    const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(marching_squares_result.surface_intersections,
                                                                             marching_squares_result.surface_segments,
                                                                           {},
                                                                             mesh);

    {
        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c = 0; c < conform_result.surface_intersection_graph.vertices.size(); ++c)
            {
                const auto &vertex = conform_result.surface_intersection_graph.vertices.at(c);

                const double vertex_original_surface = vertex.original_surface_point_index;

                double vertex_position = 0.0;
                switch (vertex.surface_type)
                {
                case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
                    vertex_position = 0.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Positive:
                    vertex_position = +1.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Negative:
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

            exporter.Export(exportFolder + "/conformed_mesh_cell0Ds.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < conform_result.surface_intersection_graph.edges.size(); ++p)
            {
                const auto &edge = conform_result.surface_intersection_graph.edges.at(p);

                Eigen::MatrixXd segment;
                segment.resize(3, 2);
                segment.col(0) = conform_result.surface_intersection_graph.vertices.at(edge.extremes.at(0)).coordinate;
                segment.col(1) = conform_result.surface_intersection_graph.vertices.at(edge.extremes.at(1)).coordinate;

                double edge_position = 0.0;
                switch (edge.surface_type)
                {
                case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
                    edge_position = 0.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Positive:
                    edge_position = +1.0;
                    break;
                case Gedim::MeshUtilities::PointSurfacePositions::Negative:
                    edge_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddSegment(
                    segment,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &edge_position}});
            }

            exporter.Export(exportFolder + "/conformed_mesh_cell1Ds.vtu");
        }
    }
}

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingSquares_H
