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

#ifndef __TEST_meshUtilities3D_MarchingCubes_H
#define __TEST_meshUtilities3D_MarchingCubes_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshUtilities.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace GedimUnitTesting
{

TEST(TestMeshUtilities, TestMeshUtilities_3D_MarchingCubes)
{
    std::string exportFolder = "./Export/TestMeshUtilities/TestMeshUtilities_3D_MarchingCubes";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
        const Gedim::GeometryUtilities::Polyhedron cube =
            geometry_utilities.CreateCubeWithOrigin(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d parallelepipedsLengthTangent = cube.Vertices.col(1) - cube.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsHeightTangent = cube.Vertices.col(3) - cube.Vertices.col(0);
        const Eigen::Vector3d parallelepipedsWidthTangent = cube.Vertices.col(4) - cube.Vertices.col(0);

        const std::vector<double> lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(4, 0.0, 1.0, 1);
        const std::vector<double> heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(5, 0.0, 1.0, 1);
        const std::vector<double> widthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(6, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateParallelepipedMesh(cube.Vertices.col(0),
                                                       parallelepipedsLengthTangent,
                                                       parallelepipedsHeightTangent,
                                                       parallelepipedsWidthTangent,
                                                       lengthMeshCurvilinearCoordinates,
                                                       heightMeshCurvilinearCoordinates,
                                                       widthMeshCurvilinearCoordinates,
                                                       mesh);
    }

    {
        const std::string mesh_folder = exportFolder + "/original_mesh";
        Gedim::Output::CreateFolder(mesh_folder);
        Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    auto plane_x = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.4, 0.4, 0.4);
        const Eigen::Vector3d n = Eigen::Vector3d(0.0, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
    };

    auto generic_plane = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const Eigen::Vector3d p0 = Eigen::Vector3d(0.6, 0.4, 1.0);
        const Eigen::Vector3d n = Eigen::Vector3d(0.5, 1.0, 0.0).normalized();
        return (point - p0).dot(n);
    };

    auto sphere = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = 1.0;
        const double radius = 0.45;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5), z = scale * (point.z() - c - 0.5);
        return x * x + y * y + z * z - radius * radius;
    };

    auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();

        return std::sin(x) * std::cos(y) + std::sin(y) * std::cos(z) + std::sin(z) * std::cos(x) - c * scale;
    };

    auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
        return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
    };

    auto diamond = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
        return std::sin(x) * std::sin(y) * std::sin(z) + std::sin(x) * std::cos(y) * std::cos(z) +
               std::cos(x) * std::sin(y) * std::cos(z) + std::cos(x) * std::cos(y) * std::sin(z) - c * scale;
    };

    auto diamond_shifted = [&diamond](const Eigen::Vector3d &point, const double &L, const double &c) {
        return diamond(point, L, c + 0.1) * diamond(point, L, c + 0.0);
    };

    auto split_p = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
        return 1.1 * (std::sin(2.0 * x) * std::sin(z) * std::cos(y) + std::sin(2.0 * y) * std::sin(x) * std::cos(z) +
                      std::sin(2.0 * z) * std::sin(y) * std::cos(x)) -
               0.2 * (std::cos(2.0 * x) * std::cos(2.0 * y) + std::cos(2.0 * y) * std::cos(2.0 * z) +
                      std::cos(2.0 * z) * std::cos(2.0 * x)) -
               0.4 * (std::cos(2.0 * x) * std::cos(2.0 * y) * std::cos(2.0 * z)) - c * scale;
    };

    auto lidinoid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
        return std::sin(2.0 * x) * std::cos(y) * std::sin(z) + std::sin(2.0 * y) * std::cos(z) * std::sin(x) +
               std::sin(2.0 * z) * std::cos(x) * std::sin(y) - std::cos(2.0 * x) * std::cos(2.0 * y) -
               std::cos(2.0 * y) * std::cos(2.0 * z) - std::cos(2.0 * z) * std::cos(2.0 * x) - c * scale;
    };

    auto neovius = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
        return 3.0 * (std::cos(x) + std::cos(y) + std::cos(z)) + 4.0 * std::cos(x) * std::cos(y) * std::cos(z) - c * scale;
    };

    auto surface_implicit_equation =
        [&plane_x, &generic_plane, &sphere, &gyroid, &diamond, &split_p, &lidinoid, &neovius, &gyroid_shifted, &diamond_shifted](
            const Eigen::MatrixXd &points) {
            float L = 1.0f;
            float c = 0.0f;

            Eigen::VectorXd surface_values(points.cols());
            for (unsigned int p = 0; p < surface_values.size(); p++)
                surface_values[p] = sphere(points.col(p), L, c);

            return surface_values;
        };

    const auto marching_cubes_result = Gedim::MeshUtilities::MarchingCubes(geometry_utilities, surface_implicit_equation, mesh);

    {
        {
            Gedim::VTKUtilities exporter;

            for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
            {
                const auto &cell0D_surface_point = marching_cubes_result.surface_intersections.cell0Ds_surface_point.at(c);

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

            for (unsigned int p = 0; p < marching_cubes_result.surface_intersections.surface_intersection_points.size(); ++p)
            {
                const auto &surface_point = marching_cubes_result.surface_intersections.surface_intersection_points.at(p);
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

            exporter.Export(exportFolder + "/marching_cubes_points.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < marching_cubes_result.surface_segments.surface_segments.size(); ++p)
            {
                const auto &segment_ids = marching_cubes_result.surface_segments.surface_segments.at(p);
                Eigen::MatrixXd segment;
                segment.resize(3, 2);
                segment.col(0) = marching_cubes_result.surface_intersections.surface_intersection_points
                                     .at(segment_ids.indices.at(0))
                                     .coordinate;
                segment.col(1) = marching_cubes_result.surface_intersections.surface_intersection_points
                                     .at(segment_ids.indices.at(1))
                                     .coordinate;
                const double cell2D_id = segment_ids.cell2D_id;

                exporter.AddSegment(segment,
                                    {{"cell2D_id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell2D_id}});
            }

            exporter.Export(exportFolder + "/marching_cubes_segments.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int p = 0; p < marching_cubes_result.surface_triangles.surface_triangles.size(); ++p)
            {
                const auto &triangle_ids = marching_cubes_result.surface_triangles.surface_triangles.at(p);
                Eigen::Matrix3d triangle;
                triangle.col(0) = marching_cubes_result.surface_intersections.surface_intersection_points
                                      .at(triangle_ids.indices.at(0))
                                      .coordinate;
                triangle.col(1) = marching_cubes_result.surface_intersections.surface_intersection_points
                                      .at(triangle_ids.indices.at(1))
                                      .coordinate;
                triangle.col(2) = marching_cubes_result.surface_intersections.surface_intersection_points
                                      .at(triangle_ids.indices.at(2))
                                      .coordinate;
                const double cell3D_id = triangle_ids.cell3D_id;

                exporter.AddPolygon(triangle,
                                    {{"cell3D_id", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell3D_id}});
            }

            exporter.Export(exportFolder + "/marching_cubes_triangles.vtu");
        }
    }

    const auto surface_intersection_graph =
        Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_cubes_result.surface_intersections,
                                                                  marching_cubes_result.surface_segments,
                                                                  marching_cubes_result.surface_triangles,
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

            exporter.Export(graph_folder + "/graph_vertices.vtu");
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

            exporter.Export(graph_folder + "/graph_edges.vtu");
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

            exporter.Export(graph_folder + "/graph_faces.vtu");
        }

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int v = 0; v < surface_intersection_graph.volumes.size(); ++v)
            {
                const auto &volume = surface_intersection_graph.volumes.at(v);

                Eigen::MatrixXd volume_vertices;
                volume_vertices.resize(3, volume.intersection_graph_vertices.size());
                std::map<unsigned int, unsigned int> vertices_to_idx;
                for (unsigned int vv = 0; vv < volume.intersection_graph_vertices.size(); ++vv)
                {
                    volume_vertices.col(vv) =
                        surface_intersection_graph.vertices.at(volume.intersection_graph_vertices.at(vv)).coordinate;
                    vertices_to_idx.insert(std::make_pair(volume.intersection_graph_vertices.at(vv), vv));
                }

                std::vector<std::vector<unsigned int>> volume_facets(volume.intersection_graph_faces.size());
                for (unsigned int vf = 0; vf < volume.intersection_graph_faces.size(); ++vf)
                {
                    const auto volume_face = surface_intersection_graph.faces.at(volume.intersection_graph_faces.at(vf));
                    volume_facets.at(vf).resize(volume_face.intersection_graph_vertices.size());
                    for (unsigned vfv = 0; vfv < volume_face.intersection_graph_vertices.size(); ++vfv)
                        volume_facets.at(vf).at(vfv) = vertices_to_idx.at(volume_face.intersection_graph_vertices.at(vfv));
                }

                double volume_position = 0.0;
                switch (volume.surface_position)
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    volume_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    volume_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    volume_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPolyhedron(
                    volume_vertices,
                    volume_facets,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &volume_position}});
            }

            exporter.Export(graph_folder + "/graph_volumes.vtu");
        }
    }

    Gedim::MeshMatrices conformed_mesh_data;
    Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

    const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(marching_cubes_result.surface_intersections,
                                                                           marching_cubes_result.surface_segments,
                                                                           marching_cubes_result.surface_triangles,
                                                                           mesh,
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

        {
            Gedim::VTKUtilities exporter;
            for (unsigned int v = 0; v < conformed_mesh.Cell3DTotalNumber(); ++v)
            {
                const double volume_original_cell3D = conform_result.new_cell3Ds_to_original_cell3D.at(v);

                const auto volume_vertices = conformed_mesh.Cell3DVerticesCoordinates(v);
                std::map<unsigned int, unsigned int> vertices_to_idx;
                for (unsigned int vv = 0; vv < conformed_mesh.Cell3DNumberVertices(v); ++vv)
                    vertices_to_idx.insert(std::make_pair(conformed_mesh.Cell3DVertex(v, vv), vv));

                std::vector<std::vector<unsigned int>> volume_facets(conformed_mesh.Cell3DNumberFaces(v));
                for (unsigned int vf = 0; vf < conformed_mesh.Cell3DNumberFaces(v); ++vf)
                {
                    const auto face_vertices = conformed_mesh.Cell2DVertices(conformed_mesh.Cell3DFace(v, vf));
                    volume_facets.at(vf).resize(face_vertices.size());
                    for (unsigned vfv = 0; vfv < face_vertices.size(); ++vfv)
                        volume_facets.at(vf).at(vfv) = vertices_to_idx.at(face_vertices.at(vfv));
                }

                double volume_position = 0.0;
                switch (conform_result.new_cell3Ds_surface_position.at(v))
                {
                case Gedim::MeshUtilities::SurfacePositions::OnSurface:
                    volume_position = 0.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Positive:
                    volume_position = +1.0;
                    break;
                case Gedim::MeshUtilities::SurfacePositions::Negative:
                    volume_position = -1.0;
                    break;
                default:
                    throw std::runtime_error("Position not managed");
                }

                exporter.AddPolyhedron(
                    volume_vertices,
                    volume_facets,
                    {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &volume_position},
                     {"original_cell3D_index", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &volume_original_cell3D}});
            }

            exporter.Export(mesh_folder + "/conformed_mesh_cell3Ds.vtu");
        }
    }

    Gedim::MeshUtilities::CheckMesh3DConfiguration check_mesh_config;
    Gedim::MeshUtilities::CheckMesh3D(check_mesh_config, geometry_utilities, conformed_mesh);
}

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingCubes_H
