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

#ifndef __TEST_FastMarchingMethod_H
#define __TEST_FastMarchingMethod_H

#include <gmock/gmock-matchers.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "FastMarchingMethod.hpp"
#include "MeshUtilities.hpp"
#include "TriangleInterface.hpp"
#include "VTKUtilities.hpp"
#include "VtkMeshExporter.hpp"

namespace UnitTesting
{
  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Simple)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Simple";
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

      Gedim::External::MeshUtilities::CreateTriangularMesh(square,
                                                           0.01,
                                                           mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 3 ? 4 : 0;
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method_2D(geometry_utilities,
                                                                                       mesh,
                                                                                       speed,
                                                                                       source);

    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          double cell0D_status = 0.0;
          switch (result.vertices_status.at(c))
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Far:
              cell0D_status = -1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Trial:
              cell0D_status = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Verified:
              cell0D_status = +1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

          exporter.AddPoint(
                mesh.Cell0DCoordinates(c),
                {{"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                 {"time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}});
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

    }
  }

  TEST(TestFastMarchingMethod, TestFastMarchingMethod_Curves)
  {
    std::string exportFolder = "./Export/TestFastMarchingMethod/TestFastMarchingMethod_Curves";
    Gedim::Output::CreateFolder(exportFolder);

    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = 1e-8;
    geometry_utilities_config.Tolerance2D = 1e-12;
    geometry_utilities_config.Tolerance3D = 1e-10;
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Gedim::MeshMatrices mesh_data;
    Gedim::MeshMatricesDAO mesh(mesh_data);

    {
      Gedim::MeshMatrices original_mesh_data;
      Gedim::MeshMatricesDAO original_mesh(original_mesh_data);

      {
        const double domain_area = 1.0;
        const double max_relative_area = 0.001;
        const double max_cell_edge = std::sqrt(domain_area * max_relative_area);
        const unsigned int num_cells_base = std::ceil(1.0 / max_cell_edge) + 1;
        const unsigned int num_cells_height = std::ceil(1.0 / max_cell_edge) + 1;

        const auto square = geometry_utilities.CreateSquare(Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
        const Eigen::Vector3d square_LengthTangent = square.col(1) - square.col(0);
        const Eigen::Vector3d square_HeightTangent = square.col(3) - square.col(0);

        const auto lengthMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_base, 0.0, 1.0, 1);
        const auto heightMeshCurvilinearCoordinates = geometry_utilities.EquispaceCoordinates(num_cells_height, 0.0, 1.0, 1);

        Gedim::MeshUtilities::CreateRectangleMesh(square.col(0),
                                                  square_LengthTangent,
                                                  square_HeightTangent,
                                                  lengthMeshCurvilinearCoordinates,
                                                  heightMeshCurvilinearCoordinates,
                                                  original_mesh);
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
        const double radius = 0.2;
        const double x = scale * (point.x() - c - 0.5), y = scale * (point.y() - c - 0.5);
        return x * x + y * y - radius * radius;
      };

      auto gyroid = [](const Eigen::Vector3d &point, const double &L, const double &c) {
        const double scale = (2.0 * std::numbers::pi) / L;
        const double x = scale * point.x(), y = scale * point.y();
        return std::sin(x) * std::cos(y) + std::sin(y) - c * scale;
      };

      auto circle_shifted = [&circle](const Eigen::Vector3d &point, const double &L, const double &c) {
        return circle(point, L, c) * circle(point, L, c + 0.1);
      };

      auto gyroid_shifted = [&gyroid](const Eigen::Vector3d &point, const double &L, const double &c) {
        return gyroid(point, L, c + 0.1) * gyroid(point, L, c + 0.0);
      };

      auto surface_implicit_equation = [&plane_x, &generic_plane, &circle, &gyroid, &gyroid_shifted, &circle_shifted](const Eigen::MatrixXd &points) {
        Eigen::VectorXd surface_values(points.cols());
        for (unsigned int p = 0; p < surface_values.size(); p++)
          surface_values[p] = gyroid_shifted(points.col(p), 1.0, 0.0);

        return surface_values;
      };

      const auto marching_squares_result =
          Gedim::MeshUtilities::MarchingSquares(geometry_utilities,
                                                surface_implicit_equation,
                                                original_mesh);

      const auto surface_intersection_graph =
          Gedim::MeshUtilities::ComputeMeshSurfaceIntersectionGraph(marching_squares_result.surface_intersections,
                                                                    marching_squares_result.surface_segments,
                                                                    {},
                                                                    original_mesh);

      Gedim::MeshMatrices conformed_mesh_data;
      Gedim::MeshMatricesDAO conformed_mesh(conformed_mesh_data);

      const auto conform_result = Gedim::MeshUtilities::ConformMeshToSurface(original_mesh,
                                                                             surface_intersection_graph,
                                                                             conformed_mesh);

      std::vector<unsigned int> cell2D_to_mantain;
      cell2D_to_mantain.reserve(conformed_mesh.Cell2DTotalNumber());
      for (unsigned int c = 0; c < conformed_mesh.Cell2DTotalNumber(); ++c)
      {
        switch (conform_result.new_cell2Ds_surface_position.at(c))
        {
          case Gedim::MeshUtilities::SurfacePositions::OnSurface:
          case Gedim::MeshUtilities::SurfacePositions::Positive:
            cell2D_to_mantain.push_back(c);
            break;
          case Gedim::MeshUtilities::SurfacePositions::Negative:
            break;
          default:
            throw std::runtime_error("Position not managed");
        }
      }
      cell2D_to_mantain.shrink_to_fit();

      const auto mesh_filter = Gedim::MeshUtilities::FilterMesh2D(cell2D_to_mantain,
                                                                  conformed_mesh);

      Gedim::MeshMatrices filtered_mesh_data;
      Gedim::MeshMatricesDAO filtered_mesh(filtered_mesh_data);

      Gedim::MeshUtilities::ExtractMesh2D(mesh_filter.Cell0Ds,
                                          mesh_filter.Cell1Ds,
                                          mesh_filter.Cell2Ds,
                                          conformed_mesh,
                                          filtered_mesh);

      std::vector<std::vector<unsigned int>> cell2Ds_triangles(filtered_mesh.Cell2DTotalNumber());

      for (unsigned int f = 0; f < filtered_mesh.Cell2DTotalNumber(); f++)
      {
        const auto face_vertices = filtered_mesh.Cell2DVerticesCoordinates(f);

        cell2Ds_triangles[f] = geometry_utilities.PolygonTriangulationByFirstVertex(face_vertices);
      }

      mesh_data = filtered_mesh_data;

      Gedim::MeshUtilities::MakeMeshTriangularFaces(cell2Ds_triangles, mesh);
      Gedim::MeshUtilities::ComputeCell0DCell1DNeighbours(mesh);
    }

    {
      const std::string mesh_folder = exportFolder + "/mesh";
      Gedim::Output::CreateFolder(mesh_folder);
      Gedim::External::MeshUtilities::ExportMeshToVTU(mesh, mesh_folder, "mesh");
    }

    const double speed = 1.0;
    const unsigned int source = mesh.Cell0DTotalNumber() > 170 ? 170 : 0;
    const auto result = Gedim::Algorithms::FastMarchingMethod::fast_marching_method_2D(geometry_utilities,
                                                                                       mesh,
                                                                                       speed,
                                                                                       source);

    {
      {
        Gedim::VTKUtilities exporter;

        for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); ++c)
        {
          double cell0D_status = 0.0;
          switch (result.vertices_status.at(c))
          {
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Far:
              cell0D_status = -1.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Trial:
              cell0D_status = 0.0;
              break;
            case Gedim::Algorithms::FastMarchingMethod::FastMarchingMethodResult::StatusTypes::Verified:
              cell0D_status = +1.0;
              break;
            default:
              throw std::runtime_error("status not managed");
          }

          exporter.AddPoint(
                mesh.Cell0DCoordinates(c),
                {{"status", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &cell0D_status},
                 {"time", Gedim::VTPProperty::Formats::Cells, static_cast<unsigned int>(1), &result.vertices_time.at(c)}});
        }

        exporter.Export(exportFolder + "/cell0Ds_time.vtu");
      }

      {
        Gedim::VTKUtilities exporter;

        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"time", Gedim::VTPProperty::Formats::Points, static_cast<unsigned int>(result.vertices_time.size()), result.vertices_time.data()}});

        exporter.Export(exportFolder + "/cell2Ds_time.vtu");
      }

    }
  }

} // namespace UnitTesting

#endif // __TEST_FastMarchingMethod_H
