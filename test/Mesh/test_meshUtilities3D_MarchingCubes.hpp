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

    Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, exportFolder, "original_mesh");

    auto sphere = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = 1.0;
      const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5), z = scale * (point.z() - c - 0.5);
      return x * x + y * y + z * z - std::sqrt(0.05);
    };

    auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = (2.0 * std::numbers::pi) / L;
      const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
      return std::sin(x) * std::sin(y) * std::sin(z) + std::sin(x) * std::cos(y) * std::cos(z) +
          std::cos(x) * std::sin(y) * std::cos(z) + std::cos(x) * std::cos(y) * std::sin(z) - c;
    };

    auto diamond = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = (2.0 * std::numbers::pi) / L;
      const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
      return std::sin(x) * std::sin(y) * std::sin(z) + std::sin(x) * std::cos(y) * std::cos(z) +
          std::cos(x) * std::sin(y) * std::cos(z) + std::cos(x) * std::cos(y) * std::sin(z) - c;
    };

    auto split_p = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = (2.0 * std::numbers::pi) / L;
      const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
      return 1.1 * (std::sin(2.0 * x) * std::sin(z) * std::cos(y) + std::sin(2.0 * y) * std::sin(x) * std::cos(z) +
                    std::sin(2.0 * z) * std::sin(y) * std::cos(x)) -
          0.2 * (std::cos(2.0 * x) * std::cos(2.0 * y) + std::cos(2.0 * y) * std::cos(2.0 * z) +
                 std::cos(2.0 * z) * std::cos(2.0 * x)) -
          0.4 * (std::cos(2.0 * x) * std::cos(2.0 * y) * std::cos(2.0 * z)) - c;
    };

    auto lidinoid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = (2.0 * std::numbers::pi) / L;
      const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
      return std::sin(2.0 * x) * std::cos(y) * std::sin(z) + std::sin(2.0 * y) * std::cos(z) * std::sin(x) +
          std::sin(2.0 * z) * std::cos(x) * std::sin(y) - std::cos(2.0 * x) * std::cos(2.0 * y) -
          std::cos(2.0 * y) * std::cos(2.0 * z) - std::cos(2.0 * z) * std::cos(2.0 * x) - c;
    };

    auto neovius = [](const Eigen::Vector3d &point, const double &L, const double &c) {
      const double scale = (2.0 * std::numbers::pi) / L;
      const double x = scale * point.x(), y = scale * point.y(), z = scale * point.z();
      return 3.0 * (std::cos(x) + std::cos(y) + std::cos(z)) + 4.0 * std::cos(x) * std::cos(y) * std::cos(z) - c;
    };

    auto surface_implicit_equation = [&sphere, &gyroid, &diamond, &split_p, &lidinoid, &neovius](const Eigen::MatrixXd &points) {
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

    const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(marching_cubes_result.surface_intersections,
                                                                           marching_cubes_result.surface_segments,
                                                                           marching_cubes_result.surface_triangles,
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

      {
        Gedim::VTKUtilities exporter;
        for (unsigned int f = 0; f < conform_result.surface_intersection_graph.faces.size(); ++f)
        {
          const auto &face = conform_result.surface_intersection_graph.faces.at(f);

          Eigen::MatrixXd face_vertices;
          face_vertices.resize(3,
                               face.intersection_graph_vertices.size());
          for (unsigned int fv = 0; fv < face.intersection_graph_vertices.size(); ++fv)
          {
            face_vertices.col(fv) =
                conform_result.surface_intersection_graph.vertices.at(face.intersection_graph_vertices.at(fv)).coordinate;
          }


          double face_position = 0.0;
          switch (face.surface_type)
          {
            case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
              face_position = 0.0;
              break;
            case Gedim::MeshUtilities::PointSurfacePositions::Positive:
              face_position = +1.0;
              break;
            case Gedim::MeshUtilities::PointSurfacePositions::Negative:
              face_position = -1.0;
              break;
            default:
              throw std::runtime_error("Position not managed");
          }

          exporter.AddPolygon(
                face_vertices,
                {{"surface_position", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &face_position}});
        }

        exporter.Export(exportFolder + "/conformed_mesh_cell2Ds.vtu");
      }

      {
        Gedim::VTKUtilities exporter;
        for (unsigned int v = 0; v < conform_result.surface_intersection_graph.volumes.size(); ++v)
        {
          const auto &volume = conform_result.surface_intersection_graph.volumes.at(v);

          Eigen::MatrixXd volume_vertices;
          volume_vertices.resize(3,
                                 volume.intersection_graph_vertices.size());
          std::map<unsigned int, unsigned int> vertices_to_idx;
          for (unsigned int vv = 0; vv < volume.intersection_graph_vertices.size(); ++vv)
          {
            volume_vertices.col(vv) =
                conform_result.surface_intersection_graph.vertices.at(volume.intersection_graph_vertices.at(vv)).coordinate;
            vertices_to_idx.insert(std::make_pair(volume.intersection_graph_vertices.at(vv), vv));
          }

          std::vector<std::vector<unsigned int>> volume_facets(volume.intersection_graph_faces.size());
          for (unsigned int vf = 0; vf < volume.intersection_graph_faces.size(); ++vf)
          {
            const auto volume_face = conform_result.surface_intersection_graph.faces.at(volume.intersection_graph_faces.at(vf));
            volume_facets.at(vf).resize(volume_face.intersection_graph_vertices.size());
            for (unsigned vfv = 0; vfv < volume_face.intersection_graph_vertices.size(); ++vfv)
              volume_facets.at(vf).at(vfv) = vertices_to_idx.at(volume_face.intersection_graph_vertices.at(vfv));
          }

          double volume_position = 0.0;
          switch (volume.surface_type)
          {
            case Gedim::MeshUtilities::PointSurfacePositions::OnSurface:
              volume_position = 0.0;
              break;
            case Gedim::MeshUtilities::PointSurfacePositions::Positive:
              volume_position = +1.0;
              break;
            case Gedim::MeshUtilities::PointSurfacePositions::Negative:
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

        exporter.Export(exportFolder + "/conformed_mesh_cell3Ds.vtu");
      }
    }
  }

} // namespace GedimUnitTesting

#endif // __TEST_meshUtilities3D_MarchingCubes_H
